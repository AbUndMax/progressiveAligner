import java.util.*;

/**
 * Stores one unaligned or multiple aligned sequences.
 */
public class Profile {

    /**
     * The sequences currently stored in this profile.
     */
    private final ArrayList<String> sequences = new ArrayList<>();

    /**
     * Returns an empty sequence profile.
     */
    public Profile() {

    }

    /**
     * Constructor which initializes a Profile with one given sequence
     * @param sequence the sequence which should be added.
     */
    public Profile(String sequence) {
        sequences.add(sequence);
    }

    /**
     * Adds a new sequence to the end of the sequence list.
     *
     * @param sequence The sequence to add.
     */
    public void addSequenceToProfile(String sequence) {
        this.sequences.add(sequence);
    }

    /**
     *
     * @return a list of all sequences in this profile
     */
    public ArrayList<String> getSequenceList() {
        return this.sequences;
    }

    /**
     * this method calculates the consensus sequence of the Profile
     * @return consensus sequence of the profile
     */
    public String getConsensusSequence() {

        OccurrenceCounter oC = new OccurrenceCounter();
        StringBuilder consensusSequence = new StringBuilder();

        for (int i = 0; i < sequences.getFirst().length(); i++) {

            for (String sequence : sequences) {
                char aminoAcid = sequence.charAt(i);
                oC.increaseByOne(aminoAcid);
            }

            consensusSequence.append(oC.getMostFrequentAminoAcid());
            oC.resetCounter();
        }

        return consensusSequence.toString();
    }

    /**
     * this method prints the Profile to the console
     */
    public void printProfile(){
        for (String sequence : sortedSequences()) {
            System.out.println(sequence);
        }
        printMatches();
    }

    /**
     * this method prints "*" on 100% identity positions and "." on Positions in which the most abundant AA occurs more
     * than 80% of the time!
     */
    public void printMatches() {
        StringBuilder matches = new StringBuilder();
        OccurrenceCounter oC = new OccurrenceCounter();

        for (int i = 0; i < sequences.getFirst().length(); i++) {
            for (String sequence : sequences) {
                oC.increaseByOne(sequence.charAt(i));
            }

            if (oC.getFrequencyOfMostFrequentAA() == 1.0) matches.append("*");
            else if (oC.getFrequencyOfMostFrequentAA() >= 0.8) matches.append(".");
            else matches.append(" ");

            oC.resetCounter();
        }

        System.out.println(matches);
    }

    /**
     * combines two profiles while propagating gaps found in an alignment of two sequences each from one Profile.
     * @param profile1 the profile to merge with profile 2
     * @param profile2 the profile to merge with profile 1
     * @param gapsProfile1 gaps inside sequence i in profile 1 found while aligning with sequence j in profile 2
     * @param gapsProfile2 gaps inside sequence j in profile 2 found while aligning with sequence i in profile 1
     * @return a combined Profile with all gaps propagated
     */
    @SuppressWarnings("DuplicatedCode")
    public static Profile combineProfiles(Profile profile1, Profile profile2, ArrayList<Integer> gapsProfile1, ArrayList<Integer> gapsProfile2) {
        Profile combinedProfile = new Profile();

        for (String sequence : profile1.getSequenceList()) {
            StringBuilder sequenceBuilder = new StringBuilder(sequence);
            for (int gapIndex : gapsProfile1) {
                sequenceBuilder.insert(gapIndex, '-');
            }

            combinedProfile.addSequenceToProfile(sequenceBuilder.toString());
        }

        for (String sequence : profile2.getSequenceList()) {
            StringBuilder sequenceBuilder = new StringBuilder(sequence);
            for (int gapIndex : gapsProfile2) {
                sequenceBuilder.insert(gapIndex, '-');
            }

            combinedProfile.addSequenceToProfile(sequenceBuilder.toString());
        }

        return combinedProfile;
    }

    /**
     * sorts the sequences by number of gaps inside them. This is just used for a better looking printout
     * @return sorted sequences list
     */
    private List<String> sortedSequences() {
        List<String> list = new ArrayList<>(sequences);

        list.sort((s1, s2) -> {
            long count1 = s1.chars().filter(ch -> ch == '-').count();
            long count2 = s2.chars().filter(ch -> ch == '-').count();
            return Long.compare(count1, count2);
        });

        return list;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null || getClass() != obj.getClass()) {
            return false;
        }
        Profile profile = (Profile) obj;
        return sequences.equals(profile.sequences);
    }

    @Override
    public int hashCode() {
        return sequences.hashCode();
    }

    @Override
    public String toString() {
        return "Profile{" +
                "sequences=" + sequences +
                '}';
    }
}
