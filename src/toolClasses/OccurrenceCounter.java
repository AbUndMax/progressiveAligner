package toolClasses;

/**
 * This class is an occurrence counter that counts the occurrence of aminoAcids
 */
public class OccurrenceCounter {
    private final char[] aminoAcidCodes = {
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
            'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'B', 'Z', 'X', 'J',
            'U', 'O', '-'
    };
    private final int[] occurrences = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0
    };

    /**
     * returns the position of the aminoAcid in the Code array which can than be used to access the counterValue of the
     * same aminoAcid
     * @param aminoAcid the oneLetter Code of the AA which
     * @return the index of the AA in the Code array
     * @throws IllegalArgumentException if given oneLetter AA is not supported
     */
    private int getPositionOf(char aminoAcid) throws IllegalArgumentException{
        for (int i = 0; i < aminoAcidCodes.length; i++) {
            if (aminoAcidCodes[i] == aminoAcid) return i;
        }

        throw new IllegalArgumentException("<<<<<<! the given aminoAcid is not supported by FASTA !>>>>>>");
    }

    /**
     * searches the index with the highest value
     * @return index with the highest value in occurrences array
     */
    private int getPositionOfMaximum() {
        int max = 0;
        int maxPosition = 0;
        for (int i = 0; i < occurrences.length; i++) {
            if (occurrences[i] > max) maxPosition = i;
        }

        return maxPosition;
    }

    /**
     * increases the counter of an AA by one
     * @param aminoAcid the AA as oneLetterCode which should be incremented by one
     */
    public void increaseByOne(char aminoAcid) {
        int index = getPositionOf(aminoAcid);
        occurrences[index] += 1;
    }

    /**
     * computes the most frequent AminoAcid
     * @return most frequent AA as oneLetter Code
     */
    public char getMostFrequentAminoAcid() {
        int indexOfMaximum = getPositionOfMaximum();
        return aminoAcidCodes[indexOfMaximum];
    }

    /**
     * calculates the frequency of the most Frequent AA
     * @return frequency of most frequent AA
     */
    public double getFrequencyOfMostFrequentAA() {
        double totalAAs = 0;
        int numberOfOccurrencesOfMostFrequentAA = occurrences[getPositionOf(getMostFrequentAminoAcid())];

        for (int occurrence : occurrences) {
            totalAAs += occurrence;
        }

        return numberOfOccurrencesOfMostFrequentAA / totalAAs;
    }

    /**
     * resets the counter all to 0
     */
    public void resetCounter() {
        for (int i = 0; i < occurrences.length; i++) {
            occurrences[i] = 0;
        }
    }
}
