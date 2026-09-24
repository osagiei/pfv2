package bio.igm.utils.filter;

import bio.igm.entities.Reads;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Discards reads that align at least as well to the unmodified genome as they do
 * to a backsplice construct.
 *
 * Reads surviving the comparison are written to unique.sam.
 *
 * @author osagie izuogu - 05/2013
 */
public class FilterGenomicHits {

    /** Construct alignments buffered before genomic.sam is scanned. */
    public static final int BATCH_SIZE = 5_000_000;

    private final File path;
    private final Logger log;
    private final Competition.Metric metric;

    Map<String, Reads> reads = new HashMap<String, Reads>();

    private long total;
    private long unparseable;
    private long discarded;
    private long retained;

    /**
     * @param _path working directory holding ptes.sam and genomic.sam
     * @throws IOException if an input is missing or an output cannot be written
     */
    public FilterGenomicHits(String _path) throws IOException {
        this(_path, Competition.Metric.ALIGNMENT_SCORE);
    }

    /**
     * @param _path  working directory holding ptes.sam and genomic.sam
     * @param metric how a construct alignment is ranked against a genomic one
     * @throws IOException if an input is missing or an output cannot be written
     */
    public FilterGenomicHits(String _path, Competition.Metric metric) throws IOException {
        this.path = new File(_path);
        this.log = PipelineFilter.log();
        this.metric = metric;

        requireReadable("ptes.sam");
        requireReadable("genomic.sam");

        log.info("Comparing reads aligned to constructs against their genomic alignments");
        read_ptes_sam();
        log.info("Genomic filter: " + total + " construct alignment(s) read, "
                + discarded + " discarded in favour of a genomic alignment, " + retained + " retained");
        if (unparseable > 0) {
            log.warning("Skipped " + unparseable + " unparseable alignment(s) in ptes.sam");
        }
        if (retained == 0) {
            log.warning("No reads survived the genomic filter");
        }
    }

    private void requireReadable(String name) throws IOException {
        File f = new File(path, name);
        if (!f.isFile() || !f.canRead()) {
            throw new IOException("Required input is missing or unreadable: " + f.getAbsolutePath());
        }
    }

    /**
     * Reads alignments to the PTES constructs in batches, scanning genomic.sam
     * once per batch so that peak memory stays bounded.
     */
    private void read_ptes_sam() throws IOException {
        File in = new File(path, "ptes.sam");

        // Opened once and truncated here rather than appended to inside the batch
        // loop, so a re-run in an existing working directory replaces the output.
        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter out = writer("genomic-filtered-out.sam");
             BufferedWriter better = writer("genomic-better.sam");
             BufferedWriter unique = writer("unique.sam")) {

            String line;
            int buffered = 0;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') {
                    continue; // SAM header
                }
                try {
                    Reads r = new Reads(line);
                    total++;
                    // Bowtie2 can report more than one alignment per read; keep the best
                    // rather than whichever happens to come last in the file.
                    Reads existing = this.reads.get(r.getId());
                    if (existing == null) {
                        this.reads.put(r.getId(), r);
                        buffered++;
                    } else if (isBetter(r, existing)) {
                        this.reads.put(r.getId(), r);
                    }
                } catch (IllegalArgumentException e) {
                    unparseable++;
                    continue;
                }

                if (buffered % BATCH_SIZE == 0) {
                    filter_genomic_reads(out, better, unique);
                    reads = new HashMap<String, Reads>();
                }
            }
            if (!reads.isEmpty()) {
                filter_genomic_reads(out, better, unique);
            }
        }
    }

    private BufferedWriter writer(String name) throws IOException {
        return new BufferedWriter(new FileWriter(new File(path, name), false));
    }

    /**
     * Scans genomic.sam and drops any buffered construct alignment that the
     * genome explains at least as well.
     */
    private void filter_genomic_reads(BufferedWriter out, BufferedWriter better, BufferedWriter unique)
            throws IOException {

        try (BufferedReader br = new BufferedReader(new FileReader(new File(path, "genomic.sam")))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') {
                    continue;
                }
                int tab = line.indexOf('\t');
                if (tab < 0) {
                    continue;
                }
                String id = line.substring(0, tab);
                Reads r = reads.get(id);
                if (r == null) {
                    continue;
                }
                Boolean constructIsBetter = Competition.constructWins(r.getLine(), line, metric);
                if (constructIsBetter == null) {
                    continue; // genomic alignment carries no comparable tags
                }
                if (!constructIsBetter) {
                    out.write(r.getLine());
                    out.write('\n');
                    better.write(line);
                    better.write('\n');
                    reads.remove(id);
                    discarded++;
                }
            }
        }

        for (Reads p : reads.values()) {
            unique.write(p.getLine());
            unique.write('\n');
            retained++;
        }
    }

    /**
     * Ranks two alignments of the same read to the same construct set.
     */
    private static boolean isBetter(Reads candidate, Reads incumbent) {
        if (candidate.getAlignmentScore() != Reads.NO_SCORE
                && incumbent.getAlignmentScore() != Reads.NO_SCORE) {
            return candidate.getAlignmentScore() > incumbent.getAlignmentScore();
        }
        return candidate.getNM() < incumbent.getNM();
    }

    /**
     * @return the reads retained by the most recent batch
     */
    public Map<String, Reads> getReads() {
        return reads;
    }

    /**
     *
     * @param reads
     */
    public void setReads(Map<String, Reads> reads) {
        this.reads = reads;
    }

    /**
     * @return the number of construct alignments retained across all batches
     */
    public long getRetained() {
        return retained;
    }

    /**
     *
     * @return
     */
    public String getPath() {
        return path.getPath();
    }
}
