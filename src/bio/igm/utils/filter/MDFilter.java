package bio.igm.utils.filter;

import bio.igm.entities.PTES;
import bio.igm.entities.Reads;
import bio.igm.utils.init.Logging;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;

/**
 * Applies the junction span and percent identity filters.
 *
 * A read is accepted only when it aligns across the construct junction with at
 * least {@code span / 2} positions either side, every position inside that
 * window is a match, and each flank meets the percent identity threshold.
 *
 * @author Osagie
 */
public class MDFilter {

    private static final Pattern DIGITS = Pattern.compile("[0-9]+");
    private static final Pattern LETTERS = Pattern.compile("[A-Za-z]");

    Map<String, PTES> ptes = new HashMap<String, PTES>();
    Map<String, List<Reads>> canonical = new HashMap<String, List<Reads>>();

    private long accepted;
    private long rejected;
    private long acceptedWithoutSpanning;

    /** Placeholder for output columns that are undefined for a given read. */
    private static String na(String value) {
        return value == null ? "NA" : value;
    }

    /**
     * Scores a set of reads against an existing set of putative structures.
     *
     * @param putative structures keyed by construct name
     * @param reads    reads keyed by read name
     * @param span     minimum junction span, in bp
     * @param pid      minimum percent identity per flank, 0-1
     */
    public MDFilter(Map<String, PTES> putative, Map<String, Reads> reads, int span, double pid) {
        for (Reads read : reads.values()) {
            PTES p = putative.get(read.getTarget());
            if (p == null) {
                continue;
            }
            if (read.getNM() == 0) {
                p.addRead(read);
                p.setCount(p.getReads().size());
                p.setConfirmed(true);
                putative.put(p.getId(), p);
            } else {
                String[] parsedMd = parseMD(read.getMdfield(), read.getCigar(),
                        read.getRefJunction() - read.getStart());
                if (checkJunctionSpan(read, parsedMd, span, pid)) {
                    p.addRead(read);
                    p.setCount(p.getReads().size());
                    p.setSpanned(true);
                    putative.put(p.getId(), p);
                }
            }
        }
        this.ptes = putative;
    }

    /**
     * Filters reads aligned to the backsplice constructs and writes the final
     * PTES output set.
     *
     * @param reads reads surviving the genomic and transcriptomic filters
     * @param span  minimum junction span, in bp
     * @param pid   minimum percent identity per flank, 0-1
     * @param path  working directory
     * @throws IOException if an output cannot be written
     */
    public MDFilter(Map<String, Reads> reads, int span, double pid, String path) throws IOException {
        this(reads, span, pid, path, false);
    }

    /**
     * Filters reads aligned to the backsplice constructs and writes the final
     * PTES output set.
     *
     * @param reads  reads surviving the reference comparisons
     * @param span   minimum junction span, in bp
     * @param pid    minimum percent identity per flank, 0-1
     * @param path   working directory
     * @param legacy accept a read that matches the construct perfectly even when it does
     *               not span the junction, as releases up to 2.1.0 did
     * @throws IOException if an output cannot be written
     */
    public MDFilter(Map<String, Reads> reads, int span, double pid, String path, boolean legacy)
            throws IOException {
        File dir = new File(path);
        Logger log = Logging.forWorkingDir(path, MDFilter.class);
        Map<String, Integer> counts = new LinkedHashMap<String, Integer>();

        try (BufferedWriter bpR = writer(dir, "pf-supporting-reads.tab");
             BufferedWriter bw = writer(dir, "pf-structure-counts.tsv");
             BufferedWriter bwB = writer(dir, "pf-structures.bed");
             BufferedWriter bwJ = writer(dir, "pf-junctions.fa");
             BufferedWriter bwP = writer(dir, "pf-pid.tsv");
             BufferedWriter bf = writer(dir, "pf-junctional-filtered.sam")) {

            bwP.write("Read_ID\tPTES_ID\tEdit_Distance\tLeftPID\tRightPID\n");

            for (Reads read : reads.values()) {
                String[] parsedMd = parseMD(read.getMdfield(), read.getCigar(),
                        read.getRefJunction() - read.getStart());

                boolean spans = checkJunctionSpan(read, parsedMd, span, pid);
                // Releases up to 2.1.0 trusted a perfect match to the construct outright,
                // on the basis that the reference comparisons had already discarded
                // anything a reference explained as well. But a read lying wholly inside
                // one construct arm says nothing about the junction joining the arms, so
                // spanning is now required. The old behaviour is kept for reproducing
                // earlier runs.
                boolean perfectMatch = legacy && !spans && read.getNM() == 0;

                if (spans || perfectMatch) {
                    accepted++;
                    if (perfectMatch) {
                        acceptedWithoutSpanning++;
                    }
                    bpR.write(read.getLine());
                    bpR.write('\n');
                    bwP.write(read.getId() + "\t" + read.getTarget() + "\t" + read.getEditDistance()
                            + "\t" + na(read.getLeftpid()) + "\t" + na(read.getRightpid()) + "\n");
                    bwJ.write(">" + read.getId() + "\t" + read.getTarget()
                            + "\tStart: " + read.getStart()
                            + "\tJunction:" + (read.getRefJunction() - read.getStart() + read.getJunctionShift())
                            + "\t" + na(read.getJunctionSeq())
                            + "\t" + read.getJunctionShift()
                            + "\t" + na(read.getHex())
                            + "\t" + read.getEditDistance()
                            + "\t" + read.getMdfield()
                            + "\t" + read.getCigar() + "\n"
                            + read.getSequence() + "\n"
                            + na(read.getMdTransformed()) + "\n");

                    Integer seen = counts.get(read.getTarget());
                    counts.put(read.getTarget(), seen == null ? 1 : seen + 1);
                } else {
                    rejected++;
                    bf.write(read.getLine());
                    bf.write('\n');
                }
            }

            writeCounts(counts, bw, bwB, log);
        }

        log.info("Junction filter: " + accepted + " read(s) accepted, " + rejected + " rejected, "
                + counts.size() + " backsplice structure(s) supported");
        if (acceptedWithoutSpanning > 0) {
            log.info(acceptedWithoutSpanning + " of the accepted read(s) matched the construct perfectly "
                    + "without spanning the junction");
        }
    }

    /**
     * Filters reads aligned to the flanking canonical junctions, which supply the
     * denominator for junctions-per-million normalisation.
     *
     * @param path working directory holding canonical.sam
     * @param span minimum junction span, in bp
     * @param pid  minimum percent identity per flank, 0-1
     * @throws IOException if canonical.sam cannot be read or an output written
     */
    public MDFilter(String path, int span, double pid) throws IOException {
        File dir = new File(path);
        File in = new File(dir, "canonical.sam");
        Logger log = Logging.forWorkingDir(path, MDFilter.class);

        if (!in.isFile()) {
            throw new IOException("Required input is missing: " + in.getAbsolutePath());
        }

        Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
        long malformed = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter bcR = writer(dir, "pf-flanking-canonical-reads.sam");
             BufferedWriter bw = writer(dir, "pf-flanking-canonical-junctions-counts.tsv");
             BufferedWriter bwB = writer(dir, "pf-flanking-canonical-junctions.bed")) {

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '@') {
                    continue; // SAM header
                }
                Reads read;
                try {
                    read = new Reads(line);
                } catch (IllegalArgumentException e) {
                    malformed++;
                    continue;
                }

                String[] parsedMd = parseMD(read.getMdfield(), read.getCigar(),
                        read.getRefJunction() - read.getStart());
                if (checkJunctionSpan(read, parsedMd, span, pid)) {
                    accepted++;
                    bcR.write(read.getLine());
                    bcR.write('\n');
                    Integer seen = counts.get(read.getTarget());
                    counts.put(read.getTarget(), seen == null ? 1 : seen + 1);
                } else {
                    rejected++;
                }
            }

            writeCounts(counts, bw, bwB, log);
        }

        if (malformed > 0) {
            log.warning("Skipped " + malformed + " unparseable alignment(s) in " + in.getName());
        }
        log.info("Canonical junctions: " + accepted + " read(s) accepted, " + rejected + " rejected, "
                + counts.size() + " junction(s) supported");
    }

    private static BufferedWriter writer(File dir, String name) throws IOException {
        return new BufferedWriter(new FileWriter(new File(dir, name), false));
    }

    /**
     * Writes per-structure counts as a TSV and as a BED, deriving the BED
     * interval from the construct name
     * ({@code chr4:144464659-144465123_+:41:GGTC}).
     */
    private void writeCounts(Map<String, Integer> counts, BufferedWriter tsv, BufferedWriter bed, Logger log)
            throws IOException {
        long unparseable = 0;
        for (Map.Entry<String, Integer> e : counts.entrySet()) {
            String s = e.getKey();
            tsv.write(s + "\t" + e.getValue() + "\n");
            try {
                String[] str = s.split(":");
                String[] range = str[1].split("-");
                bed.write(str[0] + "\t" + range[0] + "\t" + range[1].split("_")[0]
                        + "\t" + s + "\t" + e.getValue() + "\t" + str[1].split("_")[1] + "\t" + str[3] + "\n");
            } catch (ArrayIndexOutOfBoundsException ex) {
                unparseable++;
            }
        }
        if (unparseable > 0) {
            log.warning("Could not derive BED coordinates for " + unparseable + " construct name(s)");
        }
    }

    /**
     * Expands an MD tag into a per-position alignment string of M (match) and
     * mismatched base characters, then overlays indels from the CIGAR.
     *
     * @param md       the MD tag
     * @param cigar    the CIGAR string
     * @param junction junction offset within the read
     * @return {alignment string, junction shift, unshifted count}
     */
    static String[] parseMD(String md, String cigar, int junction) {
        StringBuilder seq = new StringBuilder();

        String[] contents = md.split(":");
        String body = contents.length > 2 ? contents[2] : md;

        Matcher m = DIGITS.matcher(body);

        List<Integer> ends = new ArrayList<Integer>();
        List<Integer> readGroups = new ArrayList<Integer>();

        while (m.find()) {
            ends.add(m.end());
            readGroups.add(Integer.parseInt(m.group()));
        }

        body = body + " ";
        for (int i = 0; i < readGroups.size(); i++) {
            for (int j = 0; j < readGroups.get(i); j++) {
                seq.append('M');
            }
            seq.append(body.charAt(ends.get(i)));
        }

        String temp = seq.toString().replaceAll("[^a-zA-Z0-9]", "");

        return parseIndels(temp, cigar, junction);
    }

    /**
     * Overlays CIGAR insertions, deletions and soft clips onto the alignment
     * string and reports the resulting shift of the junction offset.
     *
     * @param seq      alignment string from the MD tag
     * @param cigar    the CIGAR string
     * @param junction junction offset within the read
     * @return {alignment string, junction shift, unshifted count}
     */
    public static String[] parseIndels(String seq, String cigar, int junction) {
        List<Integer> result;
        int shift = 0;
        int unshift = 0;
        int counter2 = 0;
        Map<Integer, Character> index = new HashMap<Integer, Character>();

        List<Integer> ind = new ArrayList<Integer>();
        StringBuilder sb = new StringBuilder(seq);

        Matcher m = LETTERS.matcher(cigar);
        while (m.find()) {
            char c = m.group().charAt(0);
            index.put(m.end(), c);
            ind.add(m.end());
        }

        int counter = 0;
        for (Integer z : ind) {
            int pos = z;
            char c = index.get(z);
            if (c == 'S') {
                c = 's';
            }
            result = processDetail(cigar, pos);

            if (c != 'M') {
                for (int a = 1; a < result.size(); a += 2) {
                    for (int i = 0; i < result.get(a); i++) {
                        int at = result.get(a - 1) + i;
                        if (at < 0 || at > sb.length()) {
                            continue;
                        }
                        sb.insert(at, c);
                        if (c == 'D') {
                            counter--;
                        }
                        counter++;
                        counter2++;
                    }

                    if (result.get(a - 1) < junction) {
                        shift = counter;
                        unshift = counter2;
                    }
                }
            }
        }

        return new String[]{sb.toString(), Integer.toString(shift), Integer.toString(unshift)};
    }

    private static List<Integer> processDetail(String cigar, int pos) {
        List<Integer> result = new ArrayList<Integer>();

        int start = 0;
        int width = 0;

        Matcher m = DIGITS.matcher(cigar.substring(0, pos));
        while (m.find()) {
            width = Integer.parseInt(m.group());
            start += width;
        }
        result.add(start - width);
        result.add(width);

        return result;
    }

    /**
     * Tests whether a read spans its construct junction well enough to count as
     * evidence, and records the per-flank percent identities on the read.
     *
     * The junction window is the {@code span + 2} alignment positions centred on
     * the junction; every position in it must be a match, and both flanks must
     * meet the percent identity threshold.
     *
     * @param read     the alignment under test
     * @param parsedMD output of {@link #parseMD}
     * @param span     minimum junction span, in bp
     * @param pid      minimum percent identity per flank, 0-1
     * @return true when the read is accepted
     */
    static boolean checkJunctionSpan(Reads read, String[] parsedMD, int span, double pid) {
        String aln = parsedMD[0];
        int shift = Integer.parseInt(parsedMD[1]);
        read.setJunctionShift(shift);
        read.setHex(shift);

        int half = span / 2;
        // Offset of the junction within the read, corrected for indels.
        int position = read.getRefJunction() - read.getStart() + shift + 1;
        int startpos = position - half;
        int lastpos = position + half;

        // The window occupies alignment indices [startpos - 1, lastpos]. Requiring
        // it to lie inside the alignment string is exactly the requirement that the
        // read carry at least `half` aligned positions on each side of the junction,
        // and it also guarantees both flanks are non-empty so the identities below
        // are well defined.
        if (startpos < 1 || lastpos + 1 > aln.length()) {
            return false;
        }

        read.setJunctionSeq(aln.substring(startpos - 1, lastpos + 1));

        boolean junctionSpan = true;
        int counter = 0;
        for (int i = startpos; i <= lastpos; i++) {
            if (aln.charAt(i) != 'M') {
                junctionSpan = false;
            }
            counter++;
        }

        String leftMD = aln.substring(0, startpos);
        String rightMD = aln.substring(lastpos);

        double leftpid = (double) StringUtils.countMatches(leftMD, "M") / leftMD.length();
        double rightpid = (double) StringUtils.countMatches(rightMD, "M") / rightMD.length();

        StringBuilder parsed = new StringBuilder(aln);
        parsed.insert(startpos - 1, "<");
        parsed.insert(lastpos + 1, ">");

        read.setLeftpid(String.format("%.3f", leftpid));
        read.setRightpid(String.format("%.3f", rightpid));
        read.setMdTransformed(parsed + "|" + counter + "|" + position + "|" + startpos + "|" + lastpos
                + "|" + String.format("%.3f", leftpid) + "|" + String.format("%.3f", rightpid));

        if (rightpid < pid || leftpid < pid) {
            junctionSpan = false;
        }
        return junctionSpan;
    }

    /**
     * @return the number of reads accepted by this filter
     */
    public long getAccepted() {
        return accepted;
    }

    /**
     * @return the number of reads rejected by this filter
     */
    public long getRejected() {
        return rejected;
    }

    /**
     * @return reads accepted because they matched the construct perfectly, even
     *         though they do not span the junction
     */
    public long getAcceptedWithoutSpanning() {
        return acceptedWithoutSpanning;
    }

    /**
     *
     * @return
     */
    public Map<String, PTES> getPtes() {
        return this.ptes;
    }

    public Map<String, List<Reads>> getCanonical() {
        return canonical;
    }
}
