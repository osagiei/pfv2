package bio.igm.utils.filter;

import bio.igm.Assert;
import bio.igm.entities.Reads;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public final class FilterTest {

    /** Reference name whose junction offset is 85. */
    private static final String TARGET = "chrT:1000-1200_+:85:GTAG";

    public static void run() throws IOException {
        junctionSpanBounds();
        junctionSpanIdentity();
        referenceComparison();
        mdParsing();
        endToEnd();
    }

    private static Reads readAt(String target, int pos, int length) {
        StringBuilder seq = new StringBuilder();
        for (int i = 0; i < length; i++) {
            seq.append('A');
        }
        return new Reads("r\t0\t" + target + "\t" + pos + "\t42\t" + length + "M\t*\t0\t0\t"
                + seq + "\t" + seq.toString().replace('A', 'I') + "\tNM:i:1\tMD:Z:" + length);
    }

    private static String repeat(char c, int n) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(c);
        }
        return sb.toString();
    }

    private static String[] aln(String alignment) {
        return new String[]{alignment, "0", "0"};
    }

    /**
     * Regression test for the junction span guard.
     *
     * The guard used to compare the read-relative junction offset against the
     * read's absolute position in the construct, which rejected every read
     * aligning beyond roughly half the construct arm even when it spanned the
     * junction generously. For a 100 bp read against an 85 bp arm that discarded
     * everything from POS 41 onwards.
     */
    private static void junctionSpanBounds() {
        String matches = repeat('M', 100);

        Assert.isTrue("a read spanning 44/56 bases either side is accepted",
                MDFilter.checkJunctionSpan(readAt(TARGET, 41, 100), aln(matches), 8, 0.85));
        Assert.isTrue("a read starting at the construct start is accepted",
                MDFilter.checkJunctionSpan(readAt(TARGET, 1, 100), aln(matches), 8, 0.85));

        // Upstream boundary: the window needs span/2 + 1 aligned positions to its left.
        Assert.isTrue("the last read with enough upstream overhang is accepted",
                MDFilter.checkJunctionSpan(readAt(TARGET, 81, 100), aln(matches), 8, 0.85));
        Assert.isFalse("one base short of the upstream requirement is rejected",
                MDFilter.checkJunctionSpan(readAt(TARGET, 82, 100), aln(matches), 8, 0.85));

        // Downstream boundary, exercised with a read too short to reach past the junction.
        String short20 = repeat('M', 20);
        Assert.isTrue("the last read with enough downstream overhang is accepted",
                MDFilter.checkJunctionSpan(readAt(TARGET, 71, 20), aln(short20), 8, 0.85));
        Assert.isFalse("one base short of the downstream requirement is rejected",
                MDFilter.checkJunctionSpan(readAt(TARGET, 70, 20), aln(short20), 8, 0.85));

        // A mismatch anywhere inside the junction window disqualifies the read.
        String mismatchInWindow = repeat('M', 44) + "A" + repeat('M', 55);
        Assert.isFalse("a mismatch inside the junction window is rejected",
                MDFilter.checkJunctionSpan(readAt(TARGET, 41, 100), aln(mismatchInWindow), 8, 0.85));

        // A wider span demands more matched positions either side.
        Assert.isTrue("span 8 accepts this read",
                MDFilter.checkJunctionSpan(readAt(TARGET, 81, 100), aln(matches), 8, 0.85));
        Assert.isFalse("span 14 rejects the same read",
                MDFilter.checkJunctionSpan(readAt(TARGET, 81, 100), aln(matches), 14, 0.85));
    }

    private static void junctionSpanIdentity() {
        // 30 mismatches confined to the left flank, outside the junction window.
        String poorLeft = repeat('X', 30) + repeat('M', 70);

        Reads read = readAt(TARGET, 41, 100);
        Assert.isFalse("a flank below the identity threshold is rejected",
                MDFilter.checkJunctionSpan(read, aln(poorLeft), 8, 0.85));
        Assert.equals("left percent identity is recorded", "0.268", read.getLeftpid());
        Assert.equals("right percent identity is recorded", "1.000", read.getRightpid());

        Assert.isTrue("the same read passes under a permissive threshold",
                MDFilter.checkJunctionSpan(readAt(TARGET, 41, 100), aln(poorLeft), 8, 0.20));

        // Both flanks are always non-empty, so neither identity can be NaN and
        // slip past the comparison below the threshold.
        Reads edge = readAt(TARGET, 81, 100);
        MDFilter.checkJunctionSpan(edge, aln(repeat('M', 100)), 8, 0.85);
        Assert.equals("left flank identity is defined at the boundary", "1.000", edge.getLeftpid());
        Assert.equals("right flank identity is defined at the boundary", "1.000", edge.getRightpid());
    }

    private static String samLine(String id, String target, String nm, String md, String... extra) {
        StringBuilder sb = new StringBuilder();
        sb.append(id).append("\t0\t").append(target)
          .append("\t1\t42\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
        for (String e : extra) {
            sb.append('\t').append(e);
        }
        sb.append('\t').append(nm).append('\t').append(md);
        return sb.toString();
    }

    private static void referenceComparison() {
        Competition.Metric nm = Competition.Metric.EDIT_DISTANCE;
        Competition.Metric as = Competition.Metric.ALIGNMENT_SCORE;

        String construct = samLine("r1", TARGET, "NM:i:0", "MD:Z:10");
        String genomicWorse = samLine("r1", "chr1", "NM:i:3", "MD:Z:3A2A3");
        String genomicEqual = samLine("r1", "chr1", "NM:i:0", "MD:Z:10");

        Assert.equals("construct alignment with fewer mismatches wins",
                Boolean.TRUE, Competition.constructWins(construct, genomicWorse, nm));
        Assert.equals("an equally good reference alignment wins",
                Boolean.FALSE, Competition.constructWins(construct, genomicEqual, nm));

        // Tag order must not matter to the comparison.
        String reordered = "r1\t0\tchr1\t1\t42\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:3A2A3\tAS:i:-6\tNM:i:3";
        Assert.equals("tags are located by prefix, not by column",
                Boolean.TRUE, Competition.constructWins(construct, reordered, nm));

        String noTags = "r1\t0\tchr1\t1\t42\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tAS:i:0";
        Assert.equals("an untagged reference alignment is not comparable",
                null, Competition.constructWins(construct, noTags, nm));

        // Soft-clipped bases contribute to neither MD nor NM, so a short perfect genomic
        // fragment looks flawless to the edit-distance rule and beats a full-length
        // construct alignment carrying two mismatches. The alignment score sees the
        // clipping penalty and ranks them the right way round.
        String fullConstruct = "r1\t0\t" + TARGET + "\t1\t42\t100M\t*\t0\t0\t"
                + repeat('A', 100) + "\t" + repeat('I', 100) + "\tAS:i:-12\tNM:i:2\tMD:Z:40A30A28";
        String clippedGenomic = "r1\t0\tchr1\t1\t42\t40M60S\t*\t0\t0\t"
                + repeat('A', 100) + "\t" + repeat('I', 100) + "\tAS:i:-30\tNM:i:0\tMD:Z:40";
        Assert.equals("edit distance ranking loses the read to a clipped genomic fragment",
                Boolean.FALSE, Competition.constructWins(fullConstruct, clippedGenomic, nm));
        Assert.equals("alignment score ranking keeps it on the construct",
                Boolean.TRUE, Competition.constructWins(fullConstruct, clippedGenomic, as));

        // Ties go to the linear explanation.
        String tie = "r1\t0\tchr1\t1\t42\t100M\t*\t0\t0\t" + repeat('A', 100) + "\t"
                + repeat('I', 100) + "\tAS:i:-12\tNM:i:2\tMD:Z:40A30A28";
        Assert.equals("an equal alignment score is not a win", Boolean.FALSE,
                Competition.constructWins(fullConstruct, tie, as));

        // Without AS on both sides, score ranking falls back to the tags that exist.
        Assert.equals("score ranking falls back when AS is absent",
                Boolean.TRUE, Competition.constructWins(construct, genomicWorse, as));

        Assert.equals("metric follows the legacy flag",
                Competition.Metric.EDIT_DISTANCE, Competition.metricFor(true));
        Assert.equals("alignment score is the default metric",
                Competition.Metric.ALIGNMENT_SCORE, Competition.metricFor(false));
    }

    private static void mdParsing() {
        String[] clean = MDFilter.parseMD("MD:Z:100", "100M", 50);
        Assert.equals("a perfect alignment expands to all matches", repeat('M', 100), clean[0]);
        Assert.equals("a perfect alignment applies no shift", "0", clean[1]);

        String[] mismatch = MDFilter.parseMD("MD:Z:4A5", "10M", 5);
        Assert.equals("a mismatch is preserved in the alignment string", "MMMMAMMMMM", mismatch[0]);

        String[] insertion = MDFilter.parseMD("MD:Z:10", "5M2I5M", 8);
        Assert.equals("an insertion lengthens the alignment string", 12, insertion[0].length());
        Assert.equals("an insertion shifts the junction offset", "2", insertion[1]);
    }

    /**
     * Drives the PTES constructor end to end over a temporary working directory
     * and checks the BED output it produces.
     */
    private static void endToEnd() throws IOException {
        File dir = Files.createTempDirectory("pfv2-mdfilter").toFile();
        dir.deleteOnExit();

        String target = "chr4:144464659-144465123_+:41:GGTC";
        Map<String, Reads> reads = new HashMap<>();
        reads.put("r1", readAt(target, 1, 100));
        reads.put("r2", readAt(target, 5, 100));
        reads.put("r3", readAt(target, 90, 20)); // junction falls outside this read

        MDFilter filter = new MDFilter(reads, 8, 0.85, dir.getPath());
        Assert.equals("two spanning reads are accepted", 2L, filter.getAccepted());
        Assert.equals("the non-spanning read is rejected", 1L, filter.getRejected());
        Assert.equals("no read is accepted on the perfect-match shortcut here",
                0L, filter.getAcceptedWithoutSpanning());

        List<String> bed = Files.readAllLines(new File(dir, "pf-structures.bed").toPath(),
                StandardCharsets.UTF_8);
        Assert.equals("one structure is reported", 1, bed.size());
        Assert.equals("BED coordinates are derived from the construct name",
                "chr4\t144464659\t144465123\t" + target + "\t2\t+\tGGTC", bed.get(0));

        Assert.isTrue("supporting reads are written",
                new File(dir, "pf-supporting-reads.tab").length() > 0);
        Assert.isTrue("rejected reads are written",
                new File(dir, "pf-junctional-filtered.sam").length() > 0);

        // Re-running must replace the previous output rather than append to it.
        new MDFilter(reads, 8, 0.85, dir.getPath());
        bed = Files.readAllLines(new File(dir, "pf-structures.bed").toPath(), StandardCharsets.UTF_8);
        Assert.equals("a second run replaces the output instead of appending", 1, bed.size());

        perfectMatchShortcut(target);
    }

    /**
     * A read matching the construct perfectly is accepted even when the junction
     * falls outside it, on the basis that the reference comparisons have already
     * discarded anything the genome or transcriptome explains as well. The count is
     * reported separately so the contribution stays visible per run.
     */
    private static void perfectMatchShortcut(String target) throws IOException {
        File dir = Files.createTempDirectory("pfv2-shortcut").toFile();
        dir.deleteOnExit();

        Map<String, Reads> perfect = new HashMap<>();
        perfect.put("p1", new Reads("p1\t0\t" + target + "\t90\t42\t20M\t*\t0\t0\t"
                + repeat('A', 20) + "\t" + repeat('I', 20) + "\tNM:i:0\tMD:Z:20"));

        MDFilter strict = new MDFilter(perfect, 8, 0.85, dir.getPath());
        Assert.equals("by default a non-spanning read is rejected however well it matches",
                0L, strict.getAccepted());
        Assert.equals("and nothing is accepted without spanning", 0L,
                strict.getAcceptedWithoutSpanning());

        File legacyDir = Files.createTempDirectory("pfv2-legacy").toFile();
        legacyDir.deleteOnExit();
        MDFilter legacy = new MDFilter(perfect, 8, 0.85, legacyDir.getPath(), true);
        Assert.equals("legacy mode accepts it", 1L, legacy.getAccepted());
        Assert.equals("and counts it separately", 1L, legacy.getAcceptedWithoutSpanning());

        List<String> pid = Files.readAllLines(new File(legacyDir, "pf-pid.tsv").toPath(),
                StandardCharsets.UTF_8);
        Assert.equals("undefined identities are written as NA, not null",
                true, pid.get(1).endsWith("\tNA\tNA"));
    }
}
