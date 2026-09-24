package bio.igm.utils.discovery;

import bio.igm.Assert;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Map;

public final class DiscoveryTest {

    public static void run() throws IOException {
        genomeParser();
        strandFlip();
        assemblyGaps();
        reverseComplement();
        chimericFiltering();
        canonicalFiltering();
    }

    /**
     * Regression test: the previous parser discarded the first record in the file
     * and everything after the last even-numbered header, so chr1 and the trailing
     * chromosomes silently vanished from the genome map.
     */
    private static void genomeParser() throws IOException {
        File fasta = File.createTempFile("pfv2-genome", ".fa");
        fasta.deleteOnExit();
        Files.write(fasta.toPath(), String.join("\n",
                ">chr1",
                "acgtacgtac",
                "GGGG",
                ">chr2 dna:chromosome chromosome:GRCh38:2:1:100:1",
                "TTTTTTTTTT",
                ">chr3",
                "CCCCCCCCCC",
                ">chr4",
                "AAAAAAAAAA",
                ">chr5",
                "GGGGGGGGGG",
                "").getBytes(StandardCharsets.UTF_8));

        Map<String, String> genome = GenerateSequenceConstructsGenome.loadGenome(fasta);

        Assert.equals("every record is loaded", 5, genome.size());
        Assert.equals("the first record is not dropped", "ACGTACGTACGGGG", genome.get("chr1"));
        Assert.equals("the last record is not dropped", "GGGGGGGGGG", genome.get("chr5"));
        Assert.equals("descriptive headers resolve to their first token",
                "TTTTTTTTTT", genome.get("chr2"));
        Assert.isTrue("middle records survive", genome.containsKey("chr3") && genome.containsKey("chr4"));
    }

    private static void assemblyGaps() {
        Assert.isFalse("a clean sequence has no gap",
                GenerateSequenceConstructsGenome.hasAssemblyGap("ACGTNNNACGT"));
        Assert.isTrue("a run of ten Ns is an assembly gap",
                GenerateSequenceConstructsGenome.hasAssemblyGap("ACGT" + "N".repeat(10) + "ACGT"));
        Assert.isTrue("lowercase Ns count too",
                GenerateSequenceConstructsGenome.hasAssemblyGap("n".repeat(12)));
    }

    private static void strandFlip() {
        Assert.equals("plus becomes minus", "chr1:100-200_-:85",
                GenerateSequenceConstructsGenome.flipStrand("chr1:100-200_+:85"));
        Assert.equals("minus becomes plus", "chr1:100-200_+:85",
                GenerateSequenceConstructsGenome.flipStrand("chr1:100-200_-:85"));
        Assert.equals("anything else is left alone", "chr1:100-200_?:85",
                GenerateSequenceConstructsGenome.flipStrand("chr1:100-200_?:85"));
    }

    private static void reverseComplement() {
        Assert.equals("plain reverse complement", "ACGT",
                GenerateSequenceConstructsGenome.reverse_complement_sequence("ACGT"));
        Assert.equals("splice signal complement", "CTAC",
                GenerateSequenceConstructsGenome.reverse_complement_sequence("GTAG"));
        Assert.equals("lowercase input is handled", "ACGT",
                GenerateSequenceConstructsGenome.reverse_complement_sequence("acgt"));
        // IUPAC ambiguity codes used to produce the literal text "null" in the output.
        Assert.equals("IUPAC codes complement properly", "YRN",
                GenerateSequenceConstructsGenome.reverse_complement_sequence("NYR"));
        Assert.equals("unknown symbols become N", "N",
                GenerateSequenceConstructsGenome.reverse_complement_sequence("?"));
    }

    private static String[] chimeric(String chrA, int brkA, String strandA,
            String chrB, int brkB, String strandB, int type, int repL, int repR) {
        return new String[]{chrA, Integer.toString(brkA), strandA, chrB, Integer.toString(brkB),
            strandB, Integer.toString(type), Integer.toString(repL), Integer.toString(repR)};
    }

    private static void chimericFiltering() throws IOException {
        File dir = Files.createTempDirectory("pfv2-psc").toFile();
        dir.deleteOnExit();
        Files.write(new File(dir, "star_Chimeric.out.junction").toPath(), new byte[0]);
        Files.write(new File(dir, "star_SJ.out.tab").toPath(), new byte[0]);

        ProcessShuffledCoordinates p = new ProcessShuffledCoordinates(dir.getPath(), 100000, 50, 85);

        // A backsplice runs donor -> acceptor against the direction of transcription:
        // on the plus strand the acceptor sits at the lower coordinate, on the minus
        // strand the donor does.
        Assert.equals("plus strand anchors",
                "chr1\t4915\t4999\t1001\t1085\tchr1:1000-4999_+",
                p.chimericRecordToAnchors(chimeric("chr1", 5000, "+", "chr1", 1000, "+", 1, 0, 0)));
        Assert.equals("minus strand anchors",
                "chr1\t4915\t4999\t1001\t1085\tchr1:1000-4999_-",
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 5000, "-", 1, 0, 0)));

        // The same coordinates in linear order are a normal splice, not a backsplice.
        Assert.equals("a forward-ordered plus strand junction is rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "+", "chr1", 5000, "+", 1, 0, 0)));
        Assert.equals("a forward-ordered minus strand junction is rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 5000, "-", "chr1", 1000, "-", 1, 0, 0)));

        Assert.equals("different chromosomes are rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr2", 5000, "-", 1, 0, 0)));
        Assert.equals("different strands are rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 5000, "+", 1, 0, 0)));
        Assert.equals("the mitochondrial genome is excluded", null,
                p.chimericRecordToAnchors(chimeric("chrM", 1000, "-", "chrM", 5000, "-", 1, 0, 0)));
        Assert.equals("repeats at the breakpoint are rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 5000, "-", 1, 2, 0)));
        Assert.equals("junction type -1 is rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 5000, "-", -1, 0, 0)));
        Assert.equals("spans below the minimum are rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 1040, "-", 1, 0, 0)));
        Assert.equals("spans above the maximum are rejected", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 200000, "-", 1, 0, 0)));

        // The maximum span must follow the configured value rather than a hardcoded
        // 1 Mb ceiling that used to override it.
        ProcessShuffledCoordinates wide = new ProcessShuffledCoordinates(dir.getPath(), 5000000, 50, 85);
        Assert.isTrue("a span above 1 Mb is retained when the maximum allows it",
                wide.chimericRecordToAnchors(
                        chimeric("chr1", 1000, "-", "chr1", 2000000, "-", 1, 0, 0)) != null);
        Assert.equals("the same span is rejected under the default maximum", null,
                p.chimericRecordToAnchors(chimeric("chr1", 1000, "-", "chr1", 2000000, "-", 1, 0, 0)));
    }

    private static void canonicalFiltering() throws IOException {
        File dir = Files.createTempDirectory("pfv2-psc2").toFile();
        dir.deleteOnExit();
        Files.write(new File(dir, "star_Chimeric.out.junction").toPath(), new byte[0]);
        Files.write(new File(dir, "star_SJ.out.tab").toPath(), new byte[0]);

        ProcessShuffledCoordinates p = new ProcessShuffledCoordinates(dir.getPath(), 100000, 50, 85);

        Assert.equals("forward strand canonical junction",
                "chr1\t915\t999\t2001\t2085\tchr1:999-2000_+",
                p.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "1", "1", "0", "9", "0", "40"}));
        Assert.equals("reverse strand canonical junction",
                "chr1\t915\t999\t2001\t2085\tchr1:999-2000_-",
                p.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "2", "2", "0", "9", "0", "40"}));
        Assert.equals("undefined strand is rejected", null,
                p.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "0", "1", "0", "9", "0", "40"}));
        // Motif codes 1-4 are GT/AG, CT/AC, GC/AG and CT/GC.
        Assert.isTrue("the GC-AG motif is accepted",
                p.canonicalRecordToAnchors(
                        new String[]{"chr1", "1000", "2000", "1", "3", "0", "9", "0", "40"}) != null);
        Assert.equals("the non-canonical motif is rejected", null,
                p.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "1", "0", "0", "9", "0", "40"}));
        Assert.equals("motifs above CT-GC are rejected", null,
                p.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "1", "5", "0", "9", "0", "40"}));
        Assert.equals("a junction with no uniquely mapping read is rejected", null,
                p.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "1", "1", "0", "0", "9", "40"}));

        // Legacy mode reproduces the old rule: motif 0 admitted, GC-AG excluded.
        ProcessShuffledCoordinates old = new ProcessShuffledCoordinates(dir.getPath(), 100000, 50, 85, true);
        Assert.isTrue("legacy mode admits the non-canonical motif",
                old.canonicalRecordToAnchors(
                        new String[]{"chr1", "1000", "2000", "1", "0", "0", "9", "0", "40"}) != null);
        Assert.equals("legacy mode excludes the GC-AG motif", null,
                old.canonicalRecordToAnchors(new String[]{"chr1", "1000", "2000", "1", "3", "0", "9", "0", "40"}));
    }
}
