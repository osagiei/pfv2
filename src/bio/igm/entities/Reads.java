package bio.igm.entities;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A single alignment record parsed from a SAM file.
 *
 * Targets vary between SAM files, e.g. NM_006699.3.2, chr12 or
 * gi|329999235|ref|NM_....  For the PTES and canonical construct SAM files the
 * reference name carries the junction offset, e.g.
 * {@code chr4:144464659-144465123_+:41:GGTC}.
 *
 * @author Osagie - 05/2013
 */
public class Reads {

    /** Number of mandatory fields in a SAM alignment line. */
    public static final int MANDATORY_FIELDS = 11;

    /** Sentinel for an alignment that carries no AS tag. */
    public static final int NO_SCORE = Integer.MIN_VALUE;

    private static final Pattern DIGITS = Pattern.compile("[0-9]+");

    String id;
    int orientation;
    String target;
    String targetRaw;
    String locus;
    int start;
    String cigar;
    String mdfield;
    String editDistance;
    String sequence;
    String quality;
    int refJunction;
    boolean spansJunction;
    String junctionSeq;
    String hex;
    String line;
    int junctionShift;
    int aligned = 0;
    String mdTransformed;
    boolean genomicMatch = false;
    String genAlignment;
    String leftpid;
    String rightpid;
    int alignmentScore = NO_SCORE;

    /**
     * Expects lines from SAM files.
     *
     * @param line a SAM alignment line
     * @throws IllegalArgumentException if the line is not a usable alignment
     */
    public Reads(String line) {
        this.line = line;
        setAttributes(line);
    }

    private void setAttributes(String line) {
        String[] attributes = line.split("\t");
        if (attributes.length < MANDATORY_FIELDS) {
            throw new IllegalArgumentException(
                    "SAM line has " + attributes.length + " fields, expected at least " + MANDATORY_FIELDS);
        }

        String md = findTag(attributes, "MD:Z:");
        String nm = findTag(attributes, "NM:i:");
        if (md == null || nm == null) {
            throw new IllegalArgumentException("SAM line is missing the MD:Z: and/or NM:i: optional tag");
        }

        setId(attributes[0]);
        setOrientation(parseIntField(attributes[1], "FLAG"));
        setTargetRaw(attributes[2]);
        setTarget(attributes[2]);
        setStart(parseIntField(attributes[3], "POS"));
        setCigar(attributes[5]);
        setMdfield(md);
        setAligned(md);
        setEditDistance(nm);
        setQuality(attributes[10]);
        setRefJunction(parseRefJunction(attributes[2]));
        setSequence(attributes[9]);
        this.alignmentScore = alignmentScore(attributes);
    }

    /**
     * Locates an optional SAM tag by its {@code TAG:TYPE:} prefix rather than by
     * column position.  Bowtie2 happens to emit MD and NM in a fixed order, but
     * the SAM specification does not guarantee it and other aligners do not.
     *
     * @param fields the tab-split SAM line
     * @param prefix the tag prefix to look for, e.g. {@code "MD:Z:"}
     * @return the whole tag field, or null when absent
     */
    public static String findTag(String[] fields, String prefix) {
        for (int i = MANDATORY_FIELDS; i < fields.length; i++) {
            if (fields[i].startsWith(prefix)) {
                return fields[i];
            }
        }
        return null;
    }

    /**
     * Sums the matched-base run lengths encoded in an MD tag.
     *
     * @param mdTag a full MD tag, e.g. {@code MD:Z:31A68}
     * @return the number of aligned (matched) nucleotides
     */
    public static int alignedFromMd(String mdTag) {
        int total = 0;
        String body = mdTag.startsWith("MD:Z:") ? mdTag.substring(5) : mdTag;
        Matcher m = DIGITS.matcher(body);
        while (m.find()) {
            total += Integer.parseInt(m.group());
        }
        return total;
    }

    /**
     * Extracts the numeric value of an {@code NM:i:} tag.
     *
     * @param nmTag a full NM tag, e.g. {@code NM:i:2}
     * @return the edit distance
     */
    public static int editDistanceFromNm(String nmTag) {
        return Integer.parseInt(nmTag.substring(5));
    }

    private static int parseIntField(String value, String name) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("SAM field " + name + " is not an integer: '" + value + "'");
        }
    }

    /**
     * Reads the junction offset out of a construct reference name of the form
     * {@code <chrom>:<start>-<stop>_<strand>:<junction>:<signal>}.
     */
    private static int parseRefJunction(String referenceName) {
        String[] temp = referenceName.split(":");
        if (temp.length < 2) {
            throw new IllegalArgumentException(
                    "Reference name '" + referenceName + "' does not carry a junction offset; "
                    + "expected <chrom>:<start>-<stop>_<strand>:<junction>:<signal>");
        }
        return parseIntField(temp[temp.length - 2], "junction offset");
    }

    /**
     * @return the edit distance of this alignment as an integer
     */
    public int getNM() {
        return editDistanceFromNm(this.editDistance);
    }

    /**
     * Bowtie2's alignment score, which unlike the MD and NM tags accounts for soft
     * clipping and gap penalties.
     *
     * @return the AS tag value, or {@link #NO_SCORE} when the alignment carries none
     */
    public int getAlignmentScore() {
        return this.alignmentScore;
    }

    /**
     * Extracts the numeric value of an {@code AS:i:} tag.
     *
     * @param fields the tab-split SAM line
     * @return the alignment score, or {@link #NO_SCORE} when absent or unparseable
     */
    public static int alignmentScore(String[] fields) {
        String tag = findTag(fields, "AS:i:");
        if (tag == null) {
            return NO_SCORE;
        }
        try {
            return Integer.parseInt(tag.substring(5));
        } catch (NumberFormatException e) {
            return NO_SCORE;
        }
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Reads other = (Reads) obj;
        if ((this.target == null) ? (other.target != null) : !this.target.equals(other.target)) {
            return false;
        }
        if ((this.sequence == null) ? (other.sequence != null) : !this.sequence.equals(other.sequence)) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 59 * hash + (this.target != null ? this.target.hashCode() : 0);
        return hash;
    }

    /**
     *
     * @return cigar
     */
    public String getCigar() {
        return this.cigar;
    }

    /**
     *
     * @param cigar
     */
    public void setCigar(String cigar) {
        this.cigar = cigar;
    }

    /**
     *
     * @return editDistance
     */
    public String getEditDistance() {
        return this.editDistance;
    }

    /**
     *
     * @param editDistance
     */
    public void setEditDistance(String editDistance) {
        this.editDistance = editDistance;
    }

    /**
     *
     * @return
     */
    public String getId() {
        return this.id;
    }

    private void setId(String id) {
        this.id = id;
    }

    /**
     *
     * @return
     */
    public String getMdfield() {
        return this.mdfield;
    }

    /**
     *
     * @param mdfield
     */
    public void setMdfield(String mdfield) {
        this.mdfield = mdfield;
    }

    /**
     *
     * @return
     */
    public int getOrientation() {
        return this.orientation;
    }

    /**
     *
     * @param orientation
     */
    public void setOrientation(int orientation) {
        this.orientation = orientation;
    }

    /**
     *
     * @return quality
     */
    public String getQuality() {
        return this.quality;
    }

    /**
     *
     * @param quality
     */
    public void setQuality(String quality) {
        this.quality = quality;
    }

    /**
     *
     * @return
     */
    public String getSequence() {
        return this.sequence;
    }

    /**
     *
     * @param sequence
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     *
     * @return
     */
    public int getStart() {
        return this.start;
    }

    /**
     *
     * @param start
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     *
     * @return
     */
    public String getLocus() {
        return this.locus;
    }

    /**
     *
     * @param locus
     */
    public void setLocus(String locus) {
        this.locus = locus;
    }

    /**
     *
     * @return
     */
    public String getTarget() {
        return this.target;
    }

    /**
     *
     * @param target
     */
    public void setTarget(String target) {
        this.target = target;
    }

    /**
     *
     * @return
     */
    public String getJunctionSeq() {
        return this.junctionSeq;
    }

    /**
     *
     * @param junctionSeq
     */
    public void setJunctionSeq(String junctionSeq) {
        this.junctionSeq = junctionSeq;
    }

    /**
     *
     * @return
     */
    public int getRefJunction() {
        return this.refJunction;
    }

    /**
     *
     * @param refJunction
     */
    public void setRefJunction(int refJunction) {
        this.refJunction = refJunction;
    }

    /**
     *
     * @return
     */
    public boolean isSpansJunction() {
        return this.spansJunction;
    }

    /**
     *
     * @param spansJunction
     */
    public void setSpansJunction(boolean spansJunction) {
        this.spansJunction = spansJunction;
    }

    /**
     *
     * @return
     */
    public String getLine() {
        return this.line;
    }

    /**
     *
     * @return
     */
    public int getJunctionShift() {
        return this.junctionShift;
    }

    /**
     *
     * @param junctionShift
     */
    public void setJunctionShift(int junctionShift) {
        this.junctionShift = junctionShift;
    }

    /**
     *
     * @return
     */
    public String getTargetRaw() {
        return this.targetRaw;
    }

    /**
     *
     * @param targetRaw
     */
    public void setTargetRaw(String targetRaw) {
        this.targetRaw = targetRaw;
    }

    /**
     *
     * @return the six bases centred on the junction, or null when the junction
     *         sits too close to either end of the read to extract them
     */
    public String getHex() {
        return this.hex;
    }

    /**
     *
     * @param hex
     */
    public void setHex(String hex) {
        this.hex = substringOrNull(hex, junctionPosition());
    }

    /**
     * Records the six bases centred on the junction.  Leaves the value null when
     * the junction is too close to either end of the read for the window to fit.
     *
     * @param shift indel shift applied to the junction offset
     */
    public void setHex(int shift) {
        this.hex = substringOrNull(getSequence(), junctionPosition());
    }

    private int junctionPosition() {
        return getRefJunction() > getStart() + 3
                ? getRefJunction() - getStart() + 1 + getJunctionShift()
                : getRefJunction() + getJunctionShift();
    }

    private static String substringOrNull(String source, int pos) {
        if (source == null || pos < 3 || pos + 3 > source.length()) {
            return null;
        }
        return source.substring(pos - 3, pos + 3);
    }

    /**
     *
     * @return
     */
    public int getAligned() {
        return this.aligned;
    }

    /**
     * Calculates the number of aligned nucleotides from the MD field.
     *
     * @param aligned a full MD tag
     */
    public void setAligned(String aligned) {
        this.aligned = alignedFromMd(aligned);
    }

    /**
     *
     * @return mdTransformed
     */
    public String getMdTransformed() {
        return this.mdTransformed;
    }

    /**
     *
     * @param mdTransformed
     */
    public void setMdTransformed(String mdTransformed) {
        this.mdTransformed = mdTransformed;
    }

    /**
     *
     * @return
     */
    public boolean isGenomicMatch() {
        return this.genomicMatch;
    }

    /**
     *
     * @param genomicMatch
     */
    public void setGenomicMatch(boolean genomicMatch) {
        this.genomicMatch = genomicMatch;
    }

    /**
     *
     * @return
     */
    public String getGenAlignment() {
        return this.genAlignment;
    }

    /**
     *
     * @param genAlignment
     */
    public void setGenAlignment(String genAlignment) {
        this.genAlignment = genAlignment;
    }

    /**
     *
     * @return
     */
    public String getLeftpid() {
        return this.leftpid;
    }

    /**
     *
     * @param leftpid
     */
    public void setLeftpid(String leftpid) {
        this.leftpid = leftpid;
    }

    /**
     *
     * @return
     */
    public String getRightpid() {
        return this.rightpid;
    }

    /**
     *
     * @param rightpid
     */
    public void setRightpid(String rightpid) {
        this.rightpid = rightpid;
    }
}
