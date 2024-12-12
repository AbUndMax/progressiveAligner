import org.junit.jupiter.api.Test;
import progressiveAligner.ToolClasses.ScoreMatrix;

import static org.junit.jupiter.api.Assertions.*;

public class TestClass {

    @Test
    public void test() {

        double percent = ScoreMatrix.computePercentIdentity("ATCGGACGTA", "ATCGGGCGTA");
        assertEquals(0.9, percent, 0.000001);
    }
}
