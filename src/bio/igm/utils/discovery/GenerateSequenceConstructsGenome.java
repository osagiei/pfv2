package bio.igm.utils.discovery;

import bio.igm.entities.PutativeStructure;
import bio.igm.utils.init.Logging;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Builds the junction-spanning sequence constructs that reads are re-mapped to.
 *
 * Emits Constructs.fa for the putative backsplices and Can.fa for the flanking
 * canonical junctions; only canonical constructs are required to carry a GT-AG
 * (or CT-AC) splice signal.
 *
 * @author osagie izuogu
 */
public class GenerateSequenceConstructsGenome {

    /** Structures buffered before a batch is written out. */
    public static final int BATCH_SIZE = 5_000_000;

    private static final String PUTATIVE_INPUT = "putative_structures.txt";
    private static final String CANONICAL_INPUT = "canonical_structures.txt";
    private static final String PUTATIVE_OUTPUT = "Constructs.fa";
    private static final String CANONICAL_OUTPUT = "Can.fa";

    /** IUPAC complement table, indexed by ASCII code; unknown symbols map to N. */
    private static final char[] COMPLEMENT = buildComplementTable();

    /**
     * Splice signals accepted for a canonical junction, read in transcript orientation:
     * GT-AG, GC-AG and AT-AC.
     */
    private static final Set<String> CANONICAL_SIGNALS =
            new HashSet<>(Arrays.asList("GTAG", "GCAG", "ATAC"));

    /**
     * Signals accepted up to release 2.1.0. CTAC is GT-AG read on the opposite strand, so
     * accepting it alongside GTAG admitted junctions whose motif contradicts the strand
     * STAR assigned them.
     */
    private static final Set<String> LEGACY_SIGNALS =
            new HashSet<>(Arrays.asList("GTAG", "CTAC"));

    /**
     * Runs of N this long mark an assembly gap. An alignment to a construct built over
     * one carries no information, so the construct is not emitted.
     */
    private static final int MAX_N_RUN = 10;

    private final File path;
    private final Logger log;
    private final boolean legacy;
    private final boolean normaliseStrand;
    private final Map<String, String> genome;

    private Map<String, Map<String, PutativeStructure>> structures = new HashMap<>();

    private long written;
    private long outOfBounds;
    private long missingChromosome;
    private long signalRejected;
    private long assemblyGap;
    private long strandNormalised;
    private final Map<String, Long> unknownChromosomes = new LinkedHashMap<>();

    /**
     * @param _path      working directory holding the structure files
     * @param _gen_fasta genome reference in FASTA format
     * @throws IOException if an input is missing or an output cannot be written
     */
    public GenerateSequenceConstructsGenome(String _path, String _gen_fasta) throws IOException {
        this(_path, _gen_fasta, false);
    }

    /**
     * @param _path      working directory holding the structure files
     * @param _gen_fasta genome reference in FASTA format
     * @param _legacy    accept the canonical splice signal set used up to release 2.1.0
     * @throws IOException if an input is missing or an output cannot be written
     */
    public GenerateSequenceConstructsGenome(String _path, String _gen_fasta, boolean _legacy)
            throws IOException {
        this(_path, _gen_fasta, _legacy, false);
    }

    /**
     * @param _path            working directory holding the structure files
     * @param _gen_fasta       genome reference in FASTA format
     * @param _legacy          accept the canonical splice signal set used up to release 2.1.0
     * @param _normaliseStrand report each backsplice on the strand its splice motif implies
     * @throws IOException if an input is missing or an output cannot be written
     */
    public GenerateSequenceConstructsGenome(String _path, String _gen_fasta, boolean _legacy,
            boolean _normaliseStrand) throws IOException {
        this.normaliseStrand = _normaliseStrand;
        this.legacy = _legacy;
        this.path = new File(_path);
        this.log = Logging.forWorkingDir(_path, GenerateSequenceConstructsGenome.class);

        File fasta = new File(_gen_fasta);
        if (!fasta.isFile() || !fasta.canRead()) {
            throw new IOException("Genome FASTA is missing or unreadable: " + fasta.getAbsolutePath());
        }
        requireReadable(PUTATIVE_INPUT);
        requireReadable(CANONICAL_INPUT);

        log.info("Loading genome from " + fasta.getAbsolutePath());
        this.genome = loadGenome(fasta);
        log.info("Loaded " + genome.size() + " sequence(s) from the genome FASTA");

        buildConstructs("ptes", PUTATIVE_INPUT, PUTATIVE_OUTPUT);
        buildConstructs("canonical", CANONICAL_INPUT, CANONICAL_OUTPUT);
    }

    private void requireReadable(String name) throws IOException {
        File f = new File(path, name);
        if (!f.isFile() || !f.canRead()) {
            throw new IOException("Required input is missing or unreadable: " + f.getAbsolutePath());
        }
    }

    /**
     * Streams a FASTA file into a map of sequence name to uppercase sequence.
     *
     * The sequence name is the first whitespace-delimited token of the header,
     * matching how STAR and Bowtie2 name their references, so that descriptive
     * headers such as {@code >1 dna:chromosome chromosome:GRCh38:1:...} still
     * resolve.
     *
     * @param fasta the genome FASTA
     * @return every record in the file
     * @throws IOException if the file cannot be read or contains no records
     */
    static Map<String, String> loadGenome(File fasta) throws IOException {
        Map<String, String> genome = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(fasta))) {
            String name = null;
            StringBuilder seq = new StringBuilder();
            String line;

            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }
                if (line.charAt(0) == '>') {
                    if (name != null) {
                        genome.put(name, seq.toString());
                    }
                    name = headerToName(line);
                    seq = new StringBuilder();
                } else if (name != null) {
                    seq.append(line.trim().toUpperCase());
                }
            }
            if (name != null) {
                genome.put(name, seq.toString());
            }
        }

        if (genome.isEmpty()) {
            throw new IOException("No FASTA records found in " + fasta.getAbsolutePath());
        }
        return genome;
    }

    private static String headerToName(String header) {
        String body = header.substring(1).trim();
        int space = body.indexOf(' ');
        int tab = body.indexOf('\t');
        int cut = space < 0 ? tab : (tab < 0 ? space : Math.min(space, tab));
        return cut < 0 ? body : body.substring(0, cut);
    }

    private void buildConstructs(String type, String inputName, String outputName) throws IOException {
        File in = new File(path, inputName);
        File out = new File(path, outputName);

        written = 0;
        outOfBounds = 0;
        missingChromosome = 0;
        signalRejected = 0;
        assemblyGap = 0;
        strandNormalised = 0;
        unknownChromosomes.clear();
        structures = new HashMap<>();

        long read = 0;
        long malformed = 0;

        log.info("Generating " + type + " constructs from " + in.getName());
        // Opened once and truncated, so that re-running in an existing working
        // directory replaces the constructs rather than appending to them.
        try (BufferedReader br = new BufferedReader(new FileReader(in));
             BufferedWriter bw = new BufferedWriter(new FileWriter(out, false))) {

            String line;
            int buffered = 0;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }
                read++;
                PutativeStructure structure;
                try {
                    structure = new PutativeStructure(line);
                } catch (IllegalArgumentException e) {
                    malformed++;
                    continue;
                }

                Map<String, PutativeStructure> byId = structures.get(structure.getChromosome());
                if (byId == null) {
                    byId = new HashMap<>();
                    structures.put(structure.getChromosome(), byId);
                }
                byId.put(structure.getId(), structure);
                buffered++;

                if (buffered % BATCH_SIZE == 0) {
                    flush(bw, type);
                }
            }
            flush(bw, type);
        }

        if (malformed > 0) {
            log.warning("Skipped " + malformed + " malformed line(s) in " + in.getName());
        }
        if (!unknownChromosomes.isEmpty()) {
            log.severe("The genome FASTA has no sequence named " + unknownChromosomes.keySet()
                    + "; " + missingChromosome + " " + type + " structure(s) were dropped. "
                    + "The FASTA sequence names must match the STAR index (e.g. 'chr1' vs '1').");
        }
        if (outOfBounds > 0) {
            log.info("Dropped " + outOfBounds + " " + type
                    + " structure(s) whose anchors fall outside the sequence bounds");
        }
        if (signalRejected > 0) {
            log.info("Dropped " + signalRejected + " canonical structure(s) without an accepted splice signal "
                    + (legacy ? LEGACY_SIGNALS : CANONICAL_SIGNALS));
        }
        if (assemblyGap > 0) {
            log.info("Dropped " + assemblyGap + " " + type + " structure(s) spanning an assembly gap");
        }
        if (strandNormalised > 0) {
            log.info("Reported " + strandNormalised + " " + type
                    + " structure(s) on the strand implied by their splice motif");
        }
        log.info("Finished " + type + " constructs: " + read + " structure(s) read, "
                + written + " construct(s) written to " + out.getName());

        if (read > 0 && written == 0) {
            throw new IOException("No " + type + " constructs could be generated from " + read
                    + " structure(s). This usually means the genome FASTA does not match the "
                    + "reference used to build the STAR index.");
        }
    }

    private void flush(BufferedWriter bw, String type) throws IOException {
        for (Map.Entry<String, Map<String, PutativeStructure>> entry : structures.entrySet()) {
            String chr = entry.getKey();
            String sequence = genome.get(chr);
            if (sequence == null) {
                missingChromosome += entry.getValue().size();
                unknownChromosomes.merge(chr, (long) entry.getValue().size(), Long::sum);
                continue;
            }
            generate_sequences(sequence, entry.getValue(), type, bw);
        }
        structures = new HashMap<>();
    }

    private void generate_sequences(String sequence, Map<String, PutativeStructure> temp_structures,
            String _type, BufferedWriter bw) throws IOException {

        boolean requireSpliceSignal = !_type.equalsIgnoreCase("ptes");

        for (PutativeStructure structure : temp_structures.values()) {
            if (!structure.fitsWithin(sequence.length())) {
                outOfBounds++;
                continue;
            }

            // -1 because Java indexes from 0, unlike the 1-based coordinates
            String seq1 = sequence.substring(structure.getStart1() - 1, structure.getStop1());
            String seq2 = sequence.substring(structure.getStart2() - 1, structure.getStop2());
            String signal = sequence.substring(structure.getStop1(), structure.getStop1() + 2)
                    + sequence.substring(structure.getStart2() - 3, structure.getStart2() - 1);

            String seq = seq1 + seq2;
            String id = structure.getId();
            if (!id.contains("+")) {
                seq = reverse_complement_sequence(seq2) + reverse_complement_sequence(seq1);
                signal = reverse_complement_sequence(signal);
            }

            // STAR reports the strand a chimeric segment aligned to, which for an
            // unstranded library carries no information about the strand the host gene is
            // on. The splice motif does: a junction whose signal reads as the reverse
            // complement of a canonical motif belongs to the other strand. Reporting it
            // that way makes the backsplice constructs consistent with the canonical ones,
            // which take their strand from the motif via SJ.out.tab.
            if (normaliseStrand && !requireSpliceSignal
                    && !CANONICAL_SIGNALS.contains(signal)
                    && CANONICAL_SIGNALS.contains(reverse_complement_sequence(signal))) {
                seq = reverse_complement_sequence(seq);
                signal = reverse_complement_sequence(signal);
                id = flipStrand(id);
                strandNormalised++;
            }

            if (requireSpliceSignal && !(legacy ? LEGACY_SIGNALS : CANONICAL_SIGNALS).contains(signal)) {
                signalRejected++;
                continue;
            }

            if (hasAssemblyGap(seq)) {
                assemblyGap++;
                continue;
            }

            bw.write('>');
            bw.write(id);
            bw.write(':');
            bw.write(signal);
            bw.write('\n');
            bw.write(seq);
            bw.write('\n');
            written++;
        }
    }

    /**
     * @return true when the sequence contains a run of at least {@link #MAX_N_RUN} Ns
     */
    static boolean hasAssemblyGap(String sequence) {
        int run = 0;
        for (int i = 0; i < sequence.length(); i++) {
            char c = sequence.charAt(i);
            if (c == 'N' || c == 'n') {
                if (++run >= MAX_N_RUN) {
                    return true;
                }
            } else {
                run = 0;
            }
        }
        return false;
    }

    /**
     * Flips the {@code _+} / {@code _-} suffix of a structure id.
     */
    static String flipStrand(String id) {
        int underscore = id.lastIndexOf('_');
        if (underscore < 0 || underscore + 1 >= id.length()) {
            return id;
        }
        char strand = id.charAt(underscore + 1);
        char flipped = strand == '+' ? '-' : (strand == '-' ? '+' : strand);
        return id.substring(0, underscore + 1) + flipped + id.substring(underscore + 2);
    }

    private static char[] buildComplementTable() {
        char[] table = new char[128];
        for (int i = 0; i < table.length; i++) {
            table[i] = 'N';
        }
        String from = "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn";
        String to = "TGCAAKYWSRMBDHVNTGCAAKYWSRMBDHVN";
        for (int i = 0; i < from.length(); i++) {
            table[from.charAt(i)] = to.charAt(i);
        }
        return table;
    }

    /**
     * Reverse complements a nucleotide sequence.  IUPAC ambiguity codes are
     * complemented properly and anything else becomes N, so that unexpected
     * characters cannot leak into the construct.
     *
     * @param sequence the sequence to reverse complement
     * @return the reverse complement, in uppercase
     */
    static String reverse_complement_sequence(String sequence) {
        char[] out = new char[sequence.length()];
        for (int i = 0, j = sequence.length() - 1; j >= 0; i++, j--) {
            char c = sequence.charAt(j);
            out[i] = c < COMPLEMENT.length ? COMPLEMENT[c] : 'N';
        }
        return new String(out);
    }

    private static void usage() {
        System.err.println("Usage: GenerateSequenceConstructsGenome <working_dir> <genome_fasta> "
                + "[legacy] [normalise]");
        System.err.println("  working_dir   directory containing putative_structures.txt and canonical_structures.txt");
        System.err.println("  genome_fasta  genome reference in FASTA format, named as in the STAR index");
        System.err.println("  legacy        1 to accept the splice signal set used up to 2.1.0");
        System.err.println("  normalise     1 to report each backsplice on the strand its motif implies");
    }

    public static void main(String[] args) {
        if (args.length < 2 || args.length > 4) {
            usage();
            System.exit(2);
        }
        try {
            boolean legacy = args.length >= 3 && !"0".equals(args[2]);
            boolean normalise = args.length == 4 && !"0".equals(args[3]);
            new GenerateSequenceConstructsGenome(args[0], args[1], legacy, normalise);
        } catch (IOException ex) {
            Logger.getLogger(GenerateSequenceConstructsGenome.class.getName())
                    .log(Level.SEVERE, ex.getMessage(), ex);
            System.exit(1);
        }
    }
}
