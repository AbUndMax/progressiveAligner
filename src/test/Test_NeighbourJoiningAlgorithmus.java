package test;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

import toolClasses.*;

import java.util.LinkedList;

public class Test_NeighbourJoiningAlgorithmus {

    @Test
    public void testNeighbourJoining() {

        LinkedList<Profile> profiles = new LinkedList<>() {{
            add(new Profile("MKTPL-VGAIQV--LQNDRAKIV"));
            add(new Profile("MKTPLVG--IQVELQND-RAKIV"));
            add(new Profile("MKTP-VGAIQVELQNDR--AKIV"));
            add(new Profile("M-TPLVGAIQV-LQNDRAKIV--"));
            add(new Profile("MKTPLVG--IQV-LQNDR-KIV-"));
            add(new Profile("MKTP--VGAIQVELQNDRA--IV"));
            add(new Profile("MK-TPLVGAIQVE-LQNDRAKIV"));
            add(new Profile("MKTPLVG-AIQVELQN--RAKIV"));
            add(new Profile("MKTP--VG-IQVELQNDRAKIV-"));
        }};

        ScoreValues.MATCH_SCORE.setValue(2);
        ScoreValues.MIS_MATCH_SCORE.setValue(1);
        ScoreValues.GAP_PENALTY.setValue(1);

        NeighbourJoiningAlgorithmus nj = new NeighbourJoiningAlgorithmus(profiles);
        nj.runAlgorithm();

        nj.printInitialDistanceMatrix();

    }

}
