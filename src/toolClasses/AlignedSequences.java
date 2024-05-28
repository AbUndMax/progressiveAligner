package toolClasses;

import java.util.ArrayList;

/**
 * Stores the results of pairwise sequence alignment.
 */
public record AlignedSequences(String alignedSequence1, String alignedSequence2, int alignmentScore,
                               ArrayList<Integer> gapsAlignedSequence1,
                               ArrayList<Integer> gapsAlignedSequence2) {
}
