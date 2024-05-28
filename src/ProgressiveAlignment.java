import toolClasses.*;

import java.util.*;

/**
 * Computes pair-guided progressive multiple sequence alignment.
 */
public class ProgressiveAlignment {

    private static boolean verbose = false;

    /**
     * Parses the input FASTA file to our desired profile objects. Each Sequence gets its own Profile.
     * @param filePath path to the FASTA file (need to have at least 2 sequences!)
     * @return a list of profiles which all holds exactly one sequence
     * @throws IllegalArgumentException if the FASTA holds only 1 sequence!
     */
    private static ArrayList<Profile> parseProfileListFromFasta(String filePath) throws IllegalArgumentException {
        ArrayList<Profile> parsedSequences = new ArrayList<>();

        List<Fasta> loadedFasta = FastaIO.readInFasta(filePath);

        if (loadedFasta.size() == 1) throw new IllegalArgumentException("This FASTA holds only 1 sequence!");

        for (Fasta fasta : loadedFasta) {
            parsedSequences.add(new Profile(fasta.sequence()));
        }

        return parsedSequences;
    }

    /**
     * computes the MSA as learned in the lectures.
     * profile-profile technique:
     * -> by default, it is computed by using random sequences
     * -> if consensusMode is enabled, consensus sequences are used!
     * @param profiles all Profils that should be aligned (initially all profils hold one sequence)
     * @return a profile with the result of the MSA
     */
    public static Profile computeMSA(ArrayList<Profile> profiles) {

        while (profiles.size() != 1) {

            int highScore = 0;
            int indexProfileI = 0;
            int indexProfileJ = 0;

            if(verbose) {
                System.out.println("## start of iteration:");
                System.out.println("number of profiles: " + profiles.size());
                System.out.println("profiles:");
                System.out.println(profiles + "\n");
            }

            // find the Profile alignment with the highest score
            for (int i = 0; i < profiles.size(); i++) {
                for (int j = i + 1; j < profiles.size(); j++) {
                    // "random sequence" picking or alternatively use consensus sequence
                    String sequence1;
                    String sequence2;

                    sequence1 = profiles.get(i).getConsensusSequence();
                    sequence2 = profiles.get(j).getConsensusSequence();

                    int profileAlignScore = SequenceAlignment.computeAlignmentScore(sequence1, sequence2);

                    if(verbose) {
                        System.out.println("current i: " + i);
                        System.out.println("current j: " + j);
                        System.out.println("high-score: " + highScore);
                        System.out.println("current score: " + profileAlignScore + "\n");
                    }


                    if (profileAlignScore > highScore) {
                        highScore = profileAlignScore;
                        indexProfileI = i;
                        indexProfileJ = j;
                    }
                }
            }

            if(verbose) {
                System.out.println("index of highest profil I: " + indexProfileI);
                System.out.println("index of highest profil J: " + indexProfileJ + "\n");
            }

            Profile profile1;
            Profile profile2;

            // the Profile with the higher index needs to be removed first, since otherwise it would move the desired
            // Profile at the lower index away from its postion
            if (indexProfileJ > indexProfileI) {
                profile2 = profiles.remove(indexProfileJ);
                profile1 = profiles.remove(indexProfileI);
            } else {
                profile1 = profiles.remove(indexProfileI);
                profile2 = profiles.remove(indexProfileJ);
            }

            // Since we are allowed to choose "random" sequences as representative for a Profile, we decided to just use
            // always the first sequence since this allows us to predict the outcome better than just picking one by random!
            profiles.add(SequenceAlignment.pairGuidedAlignment(profile1, profile2));

            if(verbose) System.out.println("## end of this iteration\n");
        }

        return profiles.getFirst();
    }

    public static void main(String[] args) {

        String pathToFasta = "";

        try {
            ScoreValues.MATCH_SCORE.setValue(Integer.parseInt(args[0]));
            ScoreValues.MIS_MATCH_SCORE.setValue(Integer.parseInt(args[1]));
            int gapPenalty = Integer.parseInt(args[2]);
            if (gapPenalty < 0) gapPenalty = gapPenalty * -1;
            ScoreValues.GAP_PENALTY.setValue(gapPenalty);
            pathToFasta = args[3];

            // optional arguments:
            for (int i = 4; i < args.length; i++) {
                if (args[i].equals("-v")) {
                    verbose = true;
                } else if (args[i].equals("-cs")) {
                    continue;
                }
            }

            System.out.println("#### The MSA is computed by comparing !CONSENSUS! sequences!\n");

        } catch (Exception e) {
            System.out.println("<<<<<<! ERROR: given arguments aren't valid !>>>>>>");
            System.out.println("<<<<<<! use: <matchScore> <misMatchScore> <gapPenalty> <pathToFasta> optionals: <-cs> <-v>");
        }

        ArrayList<Profile> initialProfiles = parseProfileListFromFasta(pathToFasta);

        Profile msaOutput = computeMSA(initialProfiles);

        msaOutput.printProfile();
    }
}