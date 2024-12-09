package progressiveAligner;

import ArgsParser.*;
import progressiveAligner.MainComponents.ProgressiveAlignment;
import progressiveAligner.MainComponents.Profile;

import java.util.LinkedList;

import static progressiveAligner.MainComponents.ProgressiveAlignment.consensusMSA;

public class Main {

    public static boolean verbose = false;
    public static int matchScore;
    public static int mismatchScore;
    public static int gapPenalty;

    public static boolean verbose() {
        return verbose;
    }

    public static void main(String[] args) {

        ArgsParser parser = new ArgsParser();
        Parameter<String> pathToFasta = parser.addMandatoryStringParameter("fastaPath", "fp", "specify the path to the fasta that holds at least 2 sequences");
        Parameter<Integer> matchScore = parser.addDefaultIntegerParameter("matchScore", "ms", "positive value of the matchScore", 4);
        Parameter<Integer> misMatchScore = parser.addDefaultIntegerParameter("misMatchScore", "mms", "positive value of the misMatchScore", 2);
        Parameter<Integer> gapPenalty = parser.addDefaultIntegerParameter("gapPenalty", "gp", "positive value of the gapPenalty", 1);

        Command useConensus = parser.addCommand("Consensus", "c", "specify to use everytime newly computed distances between Profiles using consensus sequences to decide which Profiles to align next");
        Command useNJ = parser.addCommand("NeighbourJoining", "nj", "specify to use Neighbour Joining to build a guiding Tree for Profile-Profile alignment order");
        parser.toggle(useNJ, useConensus);

        parser.parse(args);

        LinkedList<Profile> initialProfiles = ProgressiveAlignment.parseProfileListFromFasta(pathToFasta.getArgument());
        Main.matchScore = matchScore.getArgument();
        Main.mismatchScore = misMatchScore.getArgument();
        Main.gapPenalty = (gapPenalty.getArgument() > 0) ? gapPenalty.getArgument() : -1 * gapPenalty.getArgument();

        Profile result = null;
        if (useConensus.isProvided()){
            result = consensusMSA(initialProfiles);
        } else {
            result = ProgressiveAlignment.neighbourJoiningGuidedMSA(initialProfiles);
        }

        // msaOutput.printProfile();
        result.printProfile();
    }
}
