package bio.igm.entities;

import bio.igm.Assert;

public final class PutativeStructureTest {

    public static void run() {
        PutativeStructure s = new PutativeStructure("chr4\t100\t185\t200\t285\tchr4:99-185_+");
        Assert.equals("chromosome", "chr4", s.getChromosome());
        // The anchor interval 100..185 is inclusive, so the arm is 86 bases and the seam
        // sits after the 86th; recording 85 put the junction one base upstream of it.
        // ProcessShuffledCoordinates now emits intervals sized so this equals the segment.
        Assert.equals("junction offset is the length of the first arm", 86, s.getJunction());
        Assert.equals("id carries the junction offset", "chr4:99-185_+:86", s.getId());
        Assert.equals("strand", "+", s.getStrand());

        Assert.isTrue("fits inside a long enough sequence", s.fitsWithin(1000));
        Assert.isFalse("does not fit when the sequence is too short", s.fitsWithin(200));

        // Anchors that run off the start of a chromosome must be rejected, not clamped.
        PutativeStructure offStart = new PutativeStructure("chr1\t-40\t45\t60\t145\tchr1:-41-45_+");
        Assert.isFalse("negative anchor start does not fit", offStart.fitsWithin(1000));

        Assert.throwsException("short line is rejected", IllegalArgumentException.class,
                () -> new PutativeStructure("chr4\t100\t185"));
        Assert.throwsException("non-numeric coordinate is rejected", IllegalArgumentException.class,
                () -> new PutativeStructure("chr4\tx\t185\t200\t285\tchr4:99-185_+"));
    }
}
