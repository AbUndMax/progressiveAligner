package toolClasses;

import java.util.ArrayList;

/**
 * Implements methods for sequence alignment.
 */
public class SequenceAlignment {

    /**
     * Adapted version of the Needleman-Wunsch algorithm to compute optimal global sequence alignments of sequences that may already contain gaps.
     *
     * @param sequence1     The first sequence to align.
     * @param sequence2     The second sequence to align.
     * @return A {@link AlignedSequences} record that stores the aligned sequences, their score and novel gaps inserted.
     */
    public static AlignedSequences adaptedNeedlemanWunsch(String sequence1, String sequence2) {

        int[][] dpMatrix = calculateDPmatrix(sequence1, sequence2);

        // Run traceback.
        StringBuilder alignedSequenceBuilder1 = new StringBuilder();
        StringBuilder alignedSequenceBuilder2 = new StringBuilder();
        ArrayList<Integer> gapsAlignedSequence1 = new ArrayList<>();
        ArrayList<Integer> gapsAlignedSequence2 = new ArrayList<>();
        int i = sequence1.length();
        int j = sequence2.length();
        int score = dpMatrix[i][j];
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && dpMatrix[i][j] == (dpMatrix[i - 1][j - 1] + (sequence1.charAt(i - 1) == sequence2.charAt(j - 1) ? ScoreValues.MATCH_SCORE.value() : (sequence1.charAt(i - 1) == '-' || sequence2.charAt(j - 1) == '-' ? -ScoreValues.GAP_PENALTY.value() : ScoreValues.MIS_MATCH_SCORE.value())))) {
                alignedSequenceBuilder1.insert(0, sequence1.charAt(i - 1));
                alignedSequenceBuilder2.insert(0, sequence2.charAt(j - 1));
                i--;
                j--;
            } else if (j > 0 && dpMatrix[i][j] == (dpMatrix[i][j - 1] - ScoreValues.GAP_PENALTY.value())) {
                alignedSequenceBuilder1.insert(0, "-");
                alignedSequenceBuilder2.insert(0, sequence2.charAt(j - 1));
                gapsAlignedSequence1.add(0, j - 1);
                j--;
            } else if (i > 0 && dpMatrix[i][j] == (dpMatrix[i - 1][j] - ScoreValues.GAP_PENALTY.value())) {
                alignedSequenceBuilder1.insert(0, sequence1.charAt(i - 1));
                alignedSequenceBuilder2.insert(0, "-");
                gapsAlignedSequence2.add(0, i - 1);
                i--;
            }
        }
        // Return aligned sequences.
        return new AlignedSequences(alignedSequenceBuilder1.toString(), alignedSequenceBuilder2.toString(), score, gapsAlignedSequence1, gapsAlignedSequence2);
    }

    /**
     * calculates the dpMatrix of to sequences
     * @param sequence1 the first sequence to align with the second
     * @param sequence2 the second sequence to align with the first
     * @return the dpMatrix int[][]
     */
    public static int[][] calculateDPmatrix(String sequence1, String sequence2) {
        // Initialize the DP matrix.
        int[][] dpMatrix = new int[sequence1.length() + 1][sequence2.length() + 1];
        for (int i = 0; i <= sequence1.length(); i++) {
            dpMatrix[i][0] = i * -ScoreValues.GAP_PENALTY.value();
        }
        for (int j = 0; j <= sequence2.length(); j++) {
            dpMatrix[0][j] = dpMatrix[0][j] = j * -ScoreValues.GAP_PENALTY.value();
        }
        // Fill the DP matrix.
        for (int i = 1; i <= sequence1.length(); i++) {
            for (int j = 1; j <= sequence2.length(); j++) {
                dpMatrix[i][j] = Math.max(dpMatrix[i - 1][j] - ScoreValues.GAP_PENALTY.value(),
                                          Math.max(dpMatrix[i][j - 1] - ScoreValues.GAP_PENALTY.value(), dpMatrix[i - 1][j - 1] +
                                                  (sequence1.charAt(i - 1) == sequence2.charAt(j - 1) ? ScoreValues.MATCH_SCORE.value() : (sequence1.charAt(i - 1) == '-' || sequence2.charAt(j - 1) == '-' ? -ScoreValues.GAP_PENALTY.value() : ScoreValues.MIS_MATCH_SCORE.value()))));
            }
        }

        return dpMatrix;
    }

    /**
     * calculates the alignmentScore by using a dpMatrix and calling the lowest right entry
     * @param sequence1 first sequence to compute the alignmentScore for alignment with second sequence
     * @param sequence2 second sequence to compute the alignment score for alignment with first sequence
     * @return the alignmentScore of two sequences
     */
    public static int computeAlignmentScore(String sequence1, String sequence2) {
        return calculateDPmatrix(sequence1, sequence2)[sequence1.length()][sequence2.length()];
    }

    /**
     * aligns two sequences in position i (of profileI) and j (of profileJ) and returns the combined Profile of booth
     * (including propagation of the gap inside the Profile)
     * @param profileI Profile in which sequence in position i should be aligned.
     * @param profileJ Profile in which sequence in position j should be aligned.
     * @return combined profile of ProfileI and ProfileJ after alignment.
     */
    public static Profile pairGuidedAlignment(Profile profileI, Profile profileJ) {

        String sequence_i;
        String sequence_j;

        sequence_i = profileI.getConsensusSequence();
        sequence_j = profileJ.getConsensusSequence();

        AlignedSequences alignmentOutput = adaptedNeedlemanWunsch(sequence_i, sequence_j);

        return  Profile.combineProfiles(profileI, profileJ,
                                        alignmentOutput.gapsAlignedSequence1(),
                                        alignmentOutput.gapsAlignedSequence2());
    }

}
