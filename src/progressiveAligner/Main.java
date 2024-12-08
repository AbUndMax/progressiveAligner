package progressiveAligner;

import progressiveAligner.MainComponents.ProgressiveAlignment;
import progressiveAligner.MainComponents.Profile;
import progressiveAligner.ToolClasses.ScoreValues;

import java.util.LinkedList;
public class Main {

    private static boolean verbose = false;

    public static boolean verbose() {
        return verbose;
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
                    break;
                }
            }

            System.out.println("#### The MSA is computed by comparing !CONSENSUS! sequences!\n");

        } catch (Exception e) {
            System.out.println("<<<<<<! ERROR: given arguments aren't valid !>>>>>>");
            System.out.println("<<<<<<! use: <matchScore> <misMatchScore> <gapPenalty> <pathToFasta> optionals: <-cs> <-v>");
        }

        LinkedList<Profile> initialProfiles = ProgressiveAlignment.parseProfileListFromFasta(pathToFasta);

        // Profile msaOutput = consensusMSA(initialProfiles);
        Profile njOutput = ProgressiveAlignment.neighbourJoiningGuidedMSA(initialProfiles);

        // msaOutput.printProfile();
        njOutput.printProfile();
    }
}
