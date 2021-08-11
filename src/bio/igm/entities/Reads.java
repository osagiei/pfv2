package bio.igm.entities;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Osagie - 05/2013
 */
public class Reads {

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

    /**
     * Expects lines from SAM files
     * @param line
     */
    public Reads(String line) {
        this.line = line;
        setAttributes(line);
    }
    
    /*
     * NOTE: Targets vary between SAM files. Eg. NM_006699.3.2, chr12 and
     * gi|329999235|ref|NM_....
     */
    private void setAttributes(String line) {
        String[] attributes = line.split("\t");
        setId(attributes[0]);
        setOrientation(Integer.parseInt(attributes[1]));
        setTargetRaw(attributes[2]);
        setTarget(attributes[2]);

        setStart(Integer.parseInt(attributes[3]));
        setCigar(attributes[5]);
        setMdfield(attributes[(attributes.length - 2)]);
        setAligned(attributes[(attributes.length - 2)]);
        setEditDistance(attributes[(attributes.length - 3)]);
        setQuality(attributes[10]);
        String[] temp = attributes[2].split(":");
        
        setRefJunction(Integer.parseInt(temp[temp.length - 2]));
        setSequence(attributes[9]);
        
    }

    private static String processID(String string) {
        String id = "";
        String[] target = string.split("ref");
        id = target[1].replaceAll("\\..*|[^a-zA-Z0-9_]", "");
        return id;
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
     * @return
     */
    public String getHex() {
        return this.hex;
    }

    /**
     *
     * @param hex
     */
    public void setHex(String hex) {
        int pos = getRefJunction() > getStart() + 3 ? getRefJunction() - getStart() + 1 + getJunctionShift() : getRefJunction() + getJunctionShift();
        this.hex = hex.substring(pos - 3, pos + 3);
    }

    /**
     *
     * @param shift
     */
    public void setHex(int shift) {
        try {
            int pos = getRefJunction() > getStart() + 3 ? getRefJunction() - getStart() + 1 + getJunctionShift() : getRefJunction() + getJunctionShift();
            this.hex = getSequence().substring(pos - 3, pos + 3);
        } catch (Exception e) {
        }
    }

    /**
     *
     * @return
     */
    public int getAligned() {
        return this.aligned;
    }

    /**
     * Method calculates number of aligned nucleotides from
     * MD field
     * @param aligned
     */
    public void setAligned(String aligned) {
        String pattern = "[0-9]+";

        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(aligned.split(":")[2]);

        while (m.find()) {
            this.aligned += Integer.parseInt(m.group());
        }
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
