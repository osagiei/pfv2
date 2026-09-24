package bio.igm.utils.filter;

import bio.igm.entities.Reads;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Discards reads that align at least as well to the transcriptome as they do to
 * a backsplice construct.
 *
 * @author osagie izuogu - 05/2013
 */
public class FilterTranscriptomicHits {

    private final File path;
    private final Logger log;
    private final Competition.Metric metric;
    private Map<String, Reads> reads;

    private long discarded;

    /**
     * @param reads reads surviving the genomic filter
     * @param _path working directory holding transcriptomic.sam
     * @throws IOException if the input is missing or an output cannot be written
     */
    FilterTranscriptomicHits(Map<String, Reads> reads, String _path) throws IOException {
        this(reads, _path, Competition.Metric.ALIGNMENT_SCORE);
    }

    /**
     * @param reads  reads surviving the genomic filter
     * @param _path  working directory holding transcriptomic.sam
     * @param metric how a construct alignment is ranked against a transcriptomic one
     * @throws IOException if the input is missing or an output cannot be written
     */
    FilterTranscriptomicHits(Map<String, Reads> reads, String _path, Competition.Metric metric)
            throws IOException {
        this.path = new File(_path);
        this.reads = reads;
        this.log = PipelineFilter.log();
        this.metric = metric;

        File in = new File(path, "transcriptomic.sam");
        if (!in.isFile() || !in.canRead()) {
            throw new IOException("Required input is missing or unreadable: " + in.getAbsolutePath());
        }

        log.info("Comparing reads aligned to constructs against their transcriptomic alignments");
        filter_refseq_reads(in);
        log.info("Transcriptomic filter: " + discarded
                + " read(s) discarded in favour of a transcriptomic alignment, "
                + reads.size() + " retained");
        if (this.reads.isEmpty()) {
            log.warning("No reads survived the transcriptomic filter");
        }
    }

    private void filter_refseq_reads(File in) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter out = writer("transcriptomic-filtered-out.sam");
             BufferedWriter better = writer("transcriptomic-better.sam");
             BufferedWriter unique = writer("transcriptomic-unique.sam")) {

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') {
                    continue; // SAM header
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
                    continue; // transcriptomic alignment carries no comparable tags
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

            for (Reads p : reads.values()) {
                unique.write(p.getLine());
                unique.write('\n');
            }
        }
    }

    private BufferedWriter writer(String name) throws IOException {
        return new BufferedWriter(new FileWriter(new File(path, name), false));
    }

    public Map<String, Reads> getReads() {
        return reads;
    }

    public void setReads(Map<String, Reads> reads) {
        this.reads = reads;
    }

    public String getPath() {
        return path.getPath();
    }
}
