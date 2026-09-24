package bio.igm.entities;

import bio.igm.Assert;

public final class ReadsTest {

    private static final String TARGET = "chr4:144464659-144465123_+:41:GGTC";

    private static String sam(String cigar, String... tags) {
        StringBuilder sb = new StringBuilder();
        sb.append("read1\t0\t").append(TARGET).append("\t1\t42\t").append(cigar)
          .append("\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
        for (String tag : tags) {
            sb.append('\t').append(tag);
        }
        return sb.toString();
    }

    public static void run() {
        Reads r = new Reads(sam("10M", "AS:i:0", "XN:i:0", "NM:i:2", "MD:Z:4A5"));
        Assert.equals("read id", "read1", r.getId());
        Assert.equals("target", TARGET, r.getTarget());
        Assert.equals("junction offset from reference name", 41, r.getRefJunction());
        Assert.equals("edit distance", 2, r.getNM());
        Assert.equals("aligned bases from MD", 9, r.getAligned());

        // Tags are located by prefix, so a different optional-tag order still parses.
        Reads reordered = new Reads(sam("10M", "MD:Z:4A5", "NM:i:2", "AS:i:0", "YT:Z:UU"));
        Assert.equals("MD found regardless of column position", "MD:Z:4A5", reordered.getMdfield());
        Assert.equals("NM found regardless of column position", 2, reordered.getNM());

        Assert.throwsException("too few fields is rejected", IllegalArgumentException.class,
                () -> new Reads("read1\t0\tchr1\t1"));
        Assert.throwsException("missing MD/NM is rejected", IllegalArgumentException.class,
                () -> new Reads(sam("10M", "AS:i:0")));
        Assert.throwsException("reference name without a junction offset is rejected",
                IllegalArgumentException.class,
                () -> new Reads("read1\t0\tchr1\t1\t42\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tMD:Z:10"));

        // The hex window is reported as null rather than throwing when it does not fit.
        Reads edge = new Reads(sam("10M", "NM:i:0", "MD:Z:10"));
        edge.setHex(0);
        Assert.equals("hex is null when the window overruns the read", null, edge.getHex());
    }
}
