package bio.igm.utils.filter;

import bio.igm.entities.Reads;
import bio.igm.utils.init.Logging;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Drives the false-positive filters over the construct alignments.
 *
 * By default a read must beat both its genomic and its transcriptomic alignment
 * before the junction span and percent identity filters are applied; either
 * reference comparison can be run on its own instead.
 *
 * @author osagie
 */
public class PipelineFilter {

    /** Which reference comparisons run before the junction filters. */
    public enum Mode {
        /** Compare against both the genome and the transcriptome (default). */
        BOTH,
        /** Compare against the genome only. */
        GENOMIC,
        /** Compare against the transcriptome only. */
        TRANSCRIPTOMIC,
        /** Apply no reference comparison; score the raw construct alignments. */
        NONE
    }

    private static Logger LOG;

    private final File path;
    private final int jspan;
    private final double pid;
    private final Mode mode;
    private final boolean legacy;

    /**
     * @param _path          working directory
     * @param _jspan         minimum junction span, in bp
     * @param _pid           minimum percent identity per flank, 0-1
     * @param _filters       run both reference comparisons
     * @param _genomic       run the genomic comparison only (when _filters is false)
     * @param _transcriptomic run the transcriptomic comparison only (when _filters is false)
     * @throws IOException if an input is missing or an output cannot be written
     */
    public PipelineFilter(String _path, int _jspan, double _pid, boolean _filters, boolean _genomic,
            boolean _transcriptomic) throws IOException {
        this(_path, _jspan, _pid, resolveMode(_filters, _genomic, _transcriptomic), false);
    }

    /**
     * @param _path          working directory
     * @param _jspan         minimum junction span, in bp
     * @param _pid           minimum percent identity per flank, 0-1
     * @param _filters       run both reference comparisons
     * @param _genomic       run the genomic comparison only (when _filters is false)
     * @param _transcriptomic run the transcriptomic comparison only (when _filters is false)
     * @param _legacy        reproduce the filter semantics of releases up to 2.1.0
     * @throws IOException if an input is missing or an output cannot be written
     */
    public PipelineFilter(String _path, int _jspan, double _pid, boolean _filters, boolean _genomic,
            boolean _transcriptomic, boolean _legacy) throws IOException {
        this(_path, _jspan, _pid, resolveMode(_filters, _genomic, _transcriptomic), _legacy);
    }

    /**
     * @param _path  working directory
     * @param _jspan minimum junction span, in bp
     * @param _pid   minimum percent identity per flank, 0-1
     * @param _mode  which reference comparisons to run
     * @throws IOException if an input is missing or an output cannot be written
     */
    public PipelineFilter(String _path, int _jspan, double _pid, Mode _mode) throws IOException {
        this(_path, _jspan, _pid, _mode, false);
    }

    /**
     * @param _path   working directory
     * @param _jspan  minimum junction span, in bp
     * @param _pid    minimum percent identity per flank, 0-1
     * @param _mode   which reference comparisons to run
     * @param _legacy reproduce the filter semantics of releases up to 2.1.0
     * @throws IOException if an input is missing or an output cannot be written
     */
    public PipelineFilter(String _path, int _jspan, double _pid, Mode _mode, boolean _legacy)
            throws IOException {
        this.path = new File(_path);
        this.jspan = _jspan;
        this.pid = _pid;
        this.mode = _mode;
        this.legacy = _legacy;

        if (!path.isDirectory()) {
            throw new IOException("Working directory does not exist: " + path.getAbsolutePath());
        }
        if (_jspan < 2) {
            throw new IllegalArgumentException("junction span must be at least 2, got " + _jspan);
        }
        if (_jspan % 2 != 0) {
            throw new IllegalArgumentException("junction span must be an even integer, got " + _jspan);
        }
        if (_pid <= 0.0 || _pid > 1.0) {
            throw new IllegalArgumentException("percent identity must be in (0, 1], got " + _pid);
        }

        LOG = Logging.forWorkingDir(_path, PipelineFilter.class);
        LOG.info("Filtering construct alignments. Mode: " + mode + ", PID: " + pid
                + ", junction span: " + jspan
                + ", ranking: " + Competition.metricFor(legacy)
                + (legacy ? ", LEGACY semantics" : ""));
        if (legacy) {
            LOG.warning("Legacy mode reproduces the filter semantics of releases up to 2.1.0: "
                    + "non-spanning perfect matches are accepted, alignments are ranked by MD and NM "
                    + "rather than by alignment score, and canonical junctions do not compete for reads");
        }

        Map<String, Reads> processed = applyReferenceFilters();

        // A read that spans a real forward splice junction is a linear read. The genomic
        // filter cannot see this, because such a read does not align contiguously to the
        // genome at all, so the canonical constructs have to compete directly.
        if (!legacy && new File(path, "canonical.sam").isFile()) {
            processed = new FilterCanonicalHits(processed, path.getPath(), jspan, pid,
                    Competition.metricFor(legacy)).getReads();
        }

        LOG.info(processed.size() + " read(s) entering the junction span and percent identity filters");

        new MDFilter(processed, jspan, pid, path.getPath(), legacy);
        processed.clear();

        LOG.info("Processing reads mapped to the flanking canonical junctions");
        new MDFilter(path.getPath(), jspan, pid);

        LOG.info("Filtering complete");
    }

    private static Mode resolveMode(boolean filters, boolean genomic, boolean transcriptomic) {
        if (filters) {
            return Mode.BOTH;
        }
        if (genomic && transcriptomic) {
            return Mode.BOTH;
        }
        if (genomic) {
            return Mode.GENOMIC;
        }
        if (transcriptomic) {
            return Mode.TRANSCRIPTOMIC;
        }
        return Mode.NONE;
    }

    /**
     * Shared logger for the filter classes.  Never returns null, so a failure
     * while parsing arguments cannot turn into a NullPointerException.
     *
     * @return the pipeline logger
     */
    public static Logger log() {
        if (LOG == null) {
            LOG = Logger.getLogger(PipelineFilter.class.getName());
        }
        return LOG;
    }

    private Map<String, Reads> applyReferenceFilters() throws IOException {
        Competition.Metric metric = Competition.metricFor(legacy);
        switch (mode) {
            case BOTH:
                new FilterGenomicHits(path.getPath(), metric);
                // unique.sam holds every batch the genomic filter retained, so it
                // is read back rather than reusing the filter's last batch.
                return new FilterTranscriptomicHits(readSam("unique.sam"), path.getPath(), metric).getReads();

            case GENOMIC:
                new FilterGenomicHits(path.getPath(), metric);
                return readSam("unique.sam");

            case TRANSCRIPTOMIC:
                return new FilterTranscriptomicHits(readSam("ptes.sam"), path.getPath(), metric).getReads();

            case NONE:
            default:
                LOG.warning("Running with no reference comparison; expect a high false positive rate");
                return readSam("ptes.sam");
        }
    }

    /**
     * Loads a SAM file into a map keyed by read name.
     *
     * @param filename file inside the working directory
     * @return the parseable alignments it contains
     * @throws IOException if the file is missing or unreadable
     */
    private Map<String, Reads> readSam(String filename) throws IOException {
        File in = new File(path, filename);
        if (!in.isFile() || !in.canRead()) {
            throw new IOException("Required input is missing or unreadable: " + in.getAbsolutePath());
        }

        Map<String, Reads> temp = new HashMap<String, Reads>();
        long unparseable = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(in))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') {
                    continue; // SAM header
                }
                try {
                    Reads r = new Reads(line);
                    temp.put(r.getId(), r);
                } catch (IllegalArgumentException e) {
                    unparseable++;
                }
            }
        }

        if (unparseable > 0) {
            LOG.warning("Skipped " + unparseable + " unparseable alignment(s) in " + filename);
        }
        return temp;
    }

    private static void usage() {
        System.err.println("Usage: PipelineFilter <working_dir> <junction_span> <pid> "
                + "<all_filters> <genomic_only> <transcriptomic_only> [legacy]");
        System.err.println("  working_dir          directory holding ptes.sam, genomic.sam, "
                + "transcriptomic.sam and canonical.sam");
        System.err.println("  junction_span        minimum junction span in bp, even integer (e.g. 8)");
        System.err.println("  pid                  minimum percent identity per flank, 0-1 (e.g. 0.85)");
        System.err.println("  all_filters          1 to run both reference comparisons, 0 otherwise");
        System.err.println("  genomic_only         1 to compare against the genome only");
        System.err.println("  transcriptomic_only  1 to compare against the transcriptome only");
        System.err.println("  legacy               1 to reproduce the filter semantics of releases up to 2.1.0");
    }

    public static void main(String[] args) {
        if (args.length < 6 || args.length > 7) {
            usage();
            System.exit(2);
        }
        try {
            String path = args[0];
            int jspan = Integer.parseInt(args[1]);
            double pid = Double.parseDouble(args[2]);
            boolean allFilters = !"0".equals(args[3]);
            boolean genomic = !"0".equals(args[4]);
            boolean transcriptomic = !"0".equals(args[5]);
            boolean legacy = args.length == 7 && !"0".equals(args[6]);

            new PipelineFilter(path, jspan, pid, allFilters, genomic, transcriptomic, legacy);
        } catch (NumberFormatException ex) {
            System.err.println("Numeric argument expected: " + ex.getMessage());
            usage();
            System.exit(2);
        } catch (IllegalArgumentException ex) {
            System.err.println("Invalid arguments: " + ex.getMessage());
            usage();
            System.exit(2);
        } catch (IOException ex) {
            log().log(Level.SEVERE, ex.getMessage(), ex);
            System.exit(1);
        }
    }
}
