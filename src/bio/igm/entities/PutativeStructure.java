package bio.igm.entities;

/**
 * One candidate backsplice (or canonical) junction, as emitted by
 * {@code ProcessShuffledCoordinates} into putative_structures.txt /
 * canonical_structures.txt.
 *
 * Expected line layout (tab separated):
 * {@code <chrom> <start1> <stop1> <start2> <stop2> <chrom:start-stop_strand>}
 *
 * @author osagie - 05/2013
 */
public class PutativeStructure {

    /** Number of tab separated columns in a structure line. */
    public static final int FIELDS = 6;

    String chromosome, id, strand, sequence, splice_signal;
    int start1, stop1, start2, stop2, junction;

    public PutativeStructure(String line) {
        String[] f = line.split("\t");
        if (f.length < FIELDS) {
            throw new IllegalArgumentException(
                    "structure line has " + f.length + " fields, expected " + FIELDS + ": '" + line + "'");
        }
        chromosome = f[0];
        start1 = parseField(f[1], "start1", line);
        stop1 = parseField(f[2], "stop1", line);
        start2 = parseField(f[3], "start2", line);
        stop2 = parseField(f[4], "stop2", line);
        // The anchor interval is inclusive at both ends, so the first arm spans
        // (stop1 - start1 + 1) bases and the seam falls immediately after them. Recording
        // the difference alone put the junction one base upstream of the real seam, which
        // skewed the purity window applied downstream.
        junction = Math.abs(stop1 - start1) + 1;
        id = f[5] + ":" + junction;

        int underscore = id.indexOf('_');
        if (underscore < 0 || underscore + 1 >= id.length()) {
            throw new IllegalArgumentException("structure id '" + id + "' carries no _<strand> suffix");
        }
        strand = id.substring(underscore + 1).split(":")[0];
    }

    private static int parseField(String value, String name, String line) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                    "structure field " + name + " is not an integer ('" + value + "') in line: '" + line + "'");
        }
    }

    /**
     * @return true when both anchors fall inside a sequence of the given length,
     *         leaving room for the two splice signal dinucleotides
     */
    public boolean fitsWithin(int sequenceLength) {
        return start1 >= 1 && start2 >= 3
                && stop1 + 2 <= sequenceLength
                && stop2 <= sequenceLength
                && start1 <= stop1 && start2 <= stop2;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public int getStart1() {
        return start1;
    }

    public void setStart1(int start1) {
        this.start1 = start1;
    }

    public int getStop1() {
        return stop1;
    }

    public void setStop1(int stop1) {
        this.stop1 = stop1;
    }

    public int getStart2() {
        return start2;
    }

    public void setStart2(int start2) {
        this.start2 = start2;
    }

    public int getStop2() {
        return stop2;
    }

    public void setStop2(int stop2) {
        this.stop2 = stop2;
    }

    public int getJunction() {
        return junction;
    }

    public void setJunction(int junction) {
        this.junction = junction;
    }

    public String getSplice_signal() {
        return splice_signal;
    }

    public void setSplice_signal(String splice_signal) {
        this.splice_signal = splice_signal;
    }
}
