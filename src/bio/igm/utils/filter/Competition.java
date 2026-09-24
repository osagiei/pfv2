package bio.igm.utils.filter;

import bio.igm.entities.Reads;

/**
 * Decides whether a construct alignment beats a competing reference alignment.
 *
 * A read is only evidence for a backsplice if the construct explains it better than the
 * genome, the transcriptome, or a canonical junction does. A tie goes to the competitor:
 * if a linear explanation accounts for the read equally well, the linear explanation is
 * the parsimonious one.
 *
 * @author osagie
 */
public final class Competition {

    /** How two alignments of the same read are ranked. */
    public enum Metric {
        /**
         * Bowtie2's AS tag. Preferred, because it accounts for soft clipping and gap
         * penalties.
         */
        ALIGNMENT_SCORE,

        /**
         * Aligned bases from MD plus edit distance from NM, as PFv2 has always compared
         * them. Soft-clipped bases contribute to neither tag, so a clipped partial
         * alignment can outrank a full-length one; retained for reproducing earlier runs.
         */
        EDIT_DISTANCE
    }

    private Competition() {
    }

    /**
     * @param constructLine  the SAM line of the alignment to a PFv2 construct
     * @param competitorLine the SAM line of the same read against a reference
     * @param metric         how to rank the two
     * @return true when the construct wins, false when the competitor wins or ties, and
     *         null when the pair cannot be compared at all
     */
    public static Boolean constructWins(String constructLine, String competitorLine, Metric metric) {
        String[] cf = constructLine.split("\t");
        String[] xf = competitorLine.split("\t");

        if (metric == Metric.ALIGNMENT_SCORE) {
            int constructAs = Reads.alignmentScore(cf);
            int competitorAs = Reads.alignmentScore(xf);
            if (constructAs != Reads.NO_SCORE && competitorAs != Reads.NO_SCORE) {
                return constructAs > competitorAs;
            }
            // Bowtie2 always emits AS, but another aligner may not; rather than treat a
            // missing score as a loss, fall back to the tags that are present.
        }

        return byEditDistance(cf, xf);
    }

    private static Boolean byEditDistance(String[] cf, String[] xf) {
        String constructMd = Reads.findTag(cf, "MD:Z:");
        String competitorMd = Reads.findTag(xf, "MD:Z:");
        String constructNm = Reads.findTag(cf, "NM:i:");
        String competitorNm = Reads.findTag(xf, "NM:i:");

        if (constructMd == null || competitorMd == null || constructNm == null || competitorNm == null) {
            return null;
        }

        try {
            int constructAligned = Reads.alignedFromMd(constructMd);
            int competitorAligned = Reads.alignedFromMd(competitorMd);
            int constructEdits = Reads.editDistanceFromNm(constructNm);
            int competitorEdits = Reads.editDistanceFromNm(competitorNm);
            return constructAligned >= competitorAligned && constructEdits < competitorEdits;
        } catch (NumberFormatException e) {
            return null;
        }
    }

    /**
     * @param legacy true to reproduce the historical comparison
     * @return the metric to apply
     */
    public static Metric metricFor(boolean legacy) {
        return legacy ? Metric.EDIT_DISTANCE : Metric.ALIGNMENT_SCORE;
    }
}
