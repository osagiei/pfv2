package bio.igm.utils.discovery;

import bio.igm.utils.init.Logging;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Turns STAR output into candidate junction coordinates.
 *
 * Reads {@code star_Chimeric.out.junction} to call putative backsplices and
 * {@code star_SJ.out.tab} to call the flanking canonical junctions used for
 * normalisation, writing anchor coordinates for each.
 *
 * @author osagie izuogu
 */
public class ProcessShuffledCoordinates {

    /** Default lower bound on the genomic span of a backsplice, in bp. */
    public static final int DEFAULT_MIN_SPAN = 50;

    /**
     * Highest STAR intron motif code accepted for a canonical junction. Codes 1-4 are
     * GT/AG, CT/AC, GC/AG and CT/GC; code 0 is non-canonical.
     */
    private static final int MAX_MOTIF = 4;

    private static final String CHIMERIC_INPUT = "star_Chimeric.out.junction";
    private static final String CANONICAL_INPUT = "star_SJ.out.tab";
    private static final String PUTATIVE_OUTPUT = "putative_structures.txt";
    private static final String CANONICAL_OUTPUT = "canonical_structures.txt";

    private final File path;
    private final int maxSpan;
    private final int minSpan;
    private final int segmentSize;
    private final boolean legacy;
    private final Logger log;

    /**
     * @param _path        working directory holding the STAR output
     * @param _max_span    maximum genomic distance between the two anchors, in bp
     * @param _segment_size length of each construct arm, in bp
     * @throws IOException if an input is missing or an output cannot be written
     */
    public ProcessShuffledCoordinates(String _path, int _max_span, int _segment_size) throws IOException {
        this(_path, _max_span, DEFAULT_MIN_SPAN, _segment_size);
    }

    /**
     * @param _path        working directory holding the STAR output
     * @param _max_span    maximum genomic distance between the two anchors, in bp
     * @param _min_span    minimum genomic distance between the two anchors, in bp
     * @param _segment_size length of each construct arm, in bp
     * @throws IOException if an input is missing or an output cannot be written
     */
    public ProcessShuffledCoordinates(String _path, int _max_span, int _min_span, int _segment_size)
            throws IOException {
        this(_path, _max_span, _min_span, _segment_size, false);
    }

    /**
     * @param _path        working directory holding the STAR output
     * @param _max_span    maximum genomic distance between the two anchors, in bp
     * @param _min_span    minimum genomic distance between the two anchors, in bp
     * @param _segment_size length of each construct arm, in bp
     * @throws IOException if an input is missing or an output cannot be written
     */
    public ProcessShuffledCoordinates(String _path, int _max_span, int _min_span, int _segment_size,
            boolean _legacy) throws IOException {
        this.legacy = _legacy;
        this.path = new File(_path);
        this.log = Logging.forWorkingDir(_path, ProcessShuffledCoordinates.class);

        if (_segment_size < 1) {
            throw new IllegalArgumentException("segment size must be positive, got " + _segment_size);
        }
        if (_max_span < _min_span) {
            throw new IllegalArgumentException(
                    "maximum span (" + _max_span + ") is below minimum span (" + _min_span + ")");
        }
        this.maxSpan = _max_span;
        this.minSpan = _min_span;
        this.segmentSize = _segment_size;

        requireReadable(CHIMERIC_INPUT);
        requireReadable(CANONICAL_INPUT);

        log.info("Screening STAR output in " + path.getAbsolutePath()
                + " (span " + minSpan + "-" + maxSpan + " bp, arm " + segmentSize + " bp)");

        // The two inputs are independent files, so read them concurrently; any
        // failure on the worker is rethrown here rather than silently logged.
        final AtomicReference<Throwable> chimericFailure = new AtomicReference<>();
        Thread worker = new Thread(new Runnable() {
            @Override
            public void run() {
                try {
                    read_star_chimeric_junction();
                } catch (Throwable t) {
                    chimericFailure.set(t);
                }
            }
        }, "chimeric-junctions");
        worker.start();

        read_star_canonical_junctions();

        try {
            worker.join();
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while reading " + CHIMERIC_INPUT, ex);
        }

        Throwable failure = chimericFailure.get();
        if (failure != null) {
            throw new IOException("Failed to process " + CHIMERIC_INPUT + ": " + failure.getMessage(), failure);
        }
    }

    private void requireReadable(String name) throws IOException {
        File f = new File(path, name);
        if (!f.isFile() || !f.canRead()) {
            throw new IOException("Required STAR output is missing or unreadable: " + f.getAbsolutePath());
        }
    }

    private void read_star_chimeric_junction() throws IOException {
        File in = new File(path, CHIMERIC_INPUT);
        File out = new File(path, PUTATIVE_OUTPUT);

        long read = 0;
        long written = 0;
        long malformed = 0;

        log.info("Reading STAR chimeric junction file " + in.getName());
        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter bw = new BufferedWriter(new FileWriter(out))) {

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#' || line.startsWith("chr_donorA")) {
                    continue; // comment or the header emitted by newer STAR releases
                }
                read++;

                String[] f = line.split("\t");
                if (f.length < 9) {
                    malformed++;
                    continue;
                }

                String to_print;
                try {
                    to_print = chimericRecordToAnchors(f);
                } catch (NumberFormatException e) {
                    malformed++;
                    continue;
                }

                if (to_print != null) {
                    bw.write(to_print);
                    bw.write('\n');
                    written++;
                }
            }
        }

        if (malformed > 0) {
            log.warning("Skipped " + malformed + " malformed line(s) in " + in.getName());
        }
        log.info("Finished " + in.getName() + ": " + read + " chimeric alignment(s) read, "
                + written + " putative backsplice(s) retained");
        if (written == 0) {
            log.warning("No putative backsplice junctions survived filtering - "
                    + "check the span bounds (-C) and that STAR ran with chimeric detection enabled");
        }
    }

    /**
     * Applies the backsplice acceptance rules to one chimeric junction record and
     * renders the two anchor intervals, or returns null when the record is rejected.
     *
     * Column layout (STAR Chimeric.out.junction): 0 donor chrom, 1 donor breakpoint,
     * 2 donor strand, 3 acceptor chrom, 4 acceptor breakpoint, 5 acceptor strand,
     * 6 junction type, 7 left repeat length, 8 right repeat length.
     */
    String chimericRecordToAnchors(String[] f) {
        String chrL = f[0];
        String chrR = f[3];
        String oL = f[2];
        String oR = f[5];

        // mono-nucleotide repeats at the junction make the breakpoint ambiguous
        if (Integer.parseInt(f[7]) > 1 || Integer.parseInt(f[8]) > 1) {
            return null;
        }
        if (!chrL.equalsIgnoreCase(chrR)) {
            return null;
        }
        if (chrL.equalsIgnoreCase("chrM") || chrL.equalsIgnoreCase("chrMT")
                || chrL.equalsIgnoreCase("M") || chrL.equalsIgnoreCase("MT")) {
            return null;
        }
        if (!oL.equalsIgnoreCase(oR)) {
            return null;
        }
        if (Integer.parseInt(f[6]) <= 0) {
            return null; // junction type -1: breakpoint falls between the mates
        }

        int left = Integer.parseInt(f[1]);
        int right = Integer.parseInt(f[4]);
        int span = Math.abs(left - right);
        if (span > maxSpan || span < minSpan) {
            return null;
        }

        if (oL.contains("-")) {
            right -= 1;
            left += 1;
            if (left > right) {
                return null;
            }
            String id = chrL + ":" + (left - 1) + "-" + right;
            // The anchor intervals are inclusive, so an arm of exactly segmentSize bases
            // ends one short of the raw offset. Without the correction each arm carried an
            // extra base and the shortest overhang a read could have was one less than
            // the configured minimum.
            int tr = right - segmentSize + 1;
            int tl = left + segmentSize - 1;
            return chrL + "\t" + tr + "\t" + right + "\t" + left + "\t" + tl + "\t" + id + "_" + oL;
        }

        left -= 1;
        right += 1;
        if (right > left) {
            return null;
        }
        String id = chrL + ":" + (right - 1) + "-" + left;
        int tr = right + segmentSize - 1;
        int tl = left - segmentSize + 1;
        return chrL + "\t" + tl + "\t" + left + "\t" + right + "\t" + tr + "\t" + id + "_" + oL;
    }

    private void read_star_canonical_junctions() throws IOException {
        File in = new File(path, CANONICAL_INPUT);
        File out = new File(path, CANONICAL_OUTPUT);

        long read = 0;
        long written = 0;
        long malformed = 0;

        log.info("Reading STAR canonical junction file " + in.getName());
        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter bw = new BufferedWriter(new FileWriter(out))) {

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') {
                    continue;
                }
                read++;

                String[] f = line.split("\t");
                if (f.length < 5) {
                    malformed++;
                    continue;
                }

                String to_print;
                try {
                    to_print = canonicalRecordToAnchors(f);
                } catch (NumberFormatException e) {
                    malformed++;
                    continue;
                }

                if (to_print != null) {
                    bw.write(to_print);
                    bw.write('\n');
                    written++;
                }
            }
        }

        if (malformed > 0) {
            log.warning("Skipped " + malformed + " malformed line(s) in " + in.getName());
        }
        log.info("Finished " + in.getName() + ": " + read + " splice junction(s) read, "
                + written + " canonical junction(s) retained");
    }

    /**
     * Renders the two anchor intervals for one SJ.out.tab record, or null when the
     * record is rejected.
     *
     * Column layout: 0 chrom, 1 intron start, 2 intron end, 3 strand
     * (0 undefined, 1 +, 2 -), 4 intron motif.
     */
    String canonicalRecordToAnchors(String[] f) {
        int motif = Integer.parseInt(f[4]);
        if (legacy) {
            // Releases up to 2.1.0 dropped motifs above 2, which admitted motif 0
            // (non-canonical) into the normalisation denominator while excluding the
            // GC/AG pair. The comment on that test said "non-canonical splice site",
            // so the intent was the rule applied below.
            if (motif > 2) {
                return null;
            }
        } else {
            if (motif < 1 || motif > MAX_MOTIF) {
                return null; // 0 is non-canonical; above 4 is AT/AC and GT/AT
            }
            // A junction with no uniquely mapping read is not a reliable denominator.
            if (f.length > 6 && Integer.parseInt(f[6]) < 1) {
                return null;
            }
        }

        int strand = Integer.parseInt(f[3]);
        if (strand != 1 && strand != 2) {
            return null; // undefined strand
        }

        String chrL = f[0];
        int left = Integer.parseInt(f[1]) - 1;
        int right = Integer.parseInt(f[2]) + 1;

        String id = chrL + ":" + left + "-" + (right - 1);
        int tr = right + segmentSize - 1;
        int tl = left - segmentSize + 1;
        String sign = strand == 1 ? "+" : "-";

        return chrL + "\t" + tl + "\t" + left + "\t" + right + "\t" + tr + "\t" + id + "_" + sign;
    }

    private static void usage() {
        System.err.println("Usage: ProcessShuffledCoordinates <working_dir> <max_span> <segment_size> "
                + "[min_span] [legacy]");
        System.err.println("  working_dir   directory containing star_Chimeric.out.junction and star_SJ.out.tab");
        System.err.println("  max_span      maximum genomic distance between backsplice anchors, in bp");
        System.err.println("  segment_size  length of each construct arm, in bp (read length - overhang)");
        System.err.println("  min_span      minimum genomic distance between anchors, in bp (default "
                + DEFAULT_MIN_SPAN + ")");
        System.err.println("  legacy        1 to reproduce the canonical motif rule used up to 2.1.0");
    }

    public static void main(String[] args) {
        if (args.length < 3 || args.length > 5) {
            usage();
            System.exit(2);
        }
        try {
            int maxSpan = Integer.parseInt(args[1]);
            int segmentSize = Integer.parseInt(args[2]);
            int minSpan = args.length >= 4 ? Integer.parseInt(args[3]) : DEFAULT_MIN_SPAN;
            boolean legacy = args.length == 5 && !"0".equals(args[4]);
            new ProcessShuffledCoordinates(args[0], maxSpan, minSpan, segmentSize, legacy);
        } catch (NumberFormatException ex) {
            System.err.println("Numeric argument expected: " + ex.getMessage());
            usage();
            System.exit(2);
        } catch (IllegalArgumentException ex) {
            System.err.println("Invalid arguments: " + ex.getMessage());
            usage();
            System.exit(2);
        } catch (IOException ex) {
            Logger.getLogger(ProcessShuffledCoordinates.class.getName()).log(Level.SEVERE, ex.getMessage(), ex);
            System.exit(1);
        }
    }
}
