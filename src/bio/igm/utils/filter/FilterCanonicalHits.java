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
 * Discards reads that a canonical (forward-spliced) junction explains at least as well
 * as a backsplice does.
 *
 * The canonical constructs are the linear alternative hypothesis: a read that spans a
 * real forward splice junction is a linear read, however well it also happens to fit a
 * backsplice construct. The genomic filter cannot catch these, because a read spanning a
 * splice junction does not align contiguously to the genome at all.
 *
 * Only canonical alignments that themselves pass the junction span and identity filters
 * count as competitors; a read sitting inside one arm of a canonical construct is just a
 * genomic read, and the genomic filter already owns that case.
 *
 * @author osagie
 */
public class FilterCanonicalHits {

    private final File path;
    private final Logger log;
    private final Map<String, Reads> reads;

    private long spanningCanonical;
    private long discarded;

    /**
     * @param reads   reads surviving the genomic and transcriptomic filters
     * @param _path   working directory holding canonical.sam
     * @param span    minimum junction span, in bp
     * @param pid     minimum percent identity per flank, 0-1
     * @param metric  how a construct alignment is ranked against a canonical one
     * @throws IOException if canonical.sam is missing or an output cannot be written
     */
    public FilterCanonicalHits(Map<String, Reads> reads, String _path, int span, double pid,
            Competition.Metric metric) throws IOException {
        this.path = new File(_path);
        this.reads = reads;
        this.log = PipelineFilter.log();

        File in = new File(path, "canonical.sam");
        if (!in.isFile() || !in.canRead()) {
            throw new IOException("Required input is missing or unreadable: " + in.getAbsolutePath());
        }

        log.info("Comparing reads aligned to constructs against spanning canonical junction alignments");
        filter(in, span, pid, metric);
        log.info("Canonical filter: " + spanningCanonical + " canonical junction-spanning alignment(s) seen, "
                + discarded + " read(s) discarded in favour of a canonical junction, "
                + reads.size() + " retained");
    }

    private void filter(File in, int span, double pid, Competition.Metric metric) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter out = new BufferedWriter(new FileWriter(new File(path, "canonical-better.sam"), false))) {

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
                Reads construct = reads.get(id);
                if (construct == null) {
                    continue; // this read is not backsplice evidence anyway
                }

                Reads canonical;
                try {
                    canonical = new Reads(line);
                } catch (IllegalArgumentException e) {
                    continue;
                }

                String[] parsedMd = MDFilter.parseMD(canonical.getMdfield(), canonical.getCigar(),
                        canonical.getRefJunction() - canonical.getStart());
                if (!MDFilter.checkJunctionSpan(canonical, parsedMd, span, pid)) {
                    continue; // does not actually span the canonical junction
                }
                spanningCanonical++;

                Boolean constructIsBetter = Competition.constructWins(construct.getLine(), line, metric);
                if (constructIsBetter != null && !constructIsBetter) {
                    out.write(line);
                    out.write('\n');
                    reads.remove(id);
                    discarded++;
                }
            }
        }
    }

    /**
     * @return the reads that survived the comparison
     */
    public Map<String, Reads> getReads() {
        return reads;
    }

    /**
     * @return the number of reads a canonical junction explained at least as well
     */
    public long getDiscarded() {
        return discarded;
    }
}
