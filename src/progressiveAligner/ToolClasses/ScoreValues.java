package progressiveAligner.ToolClasses;

/**
 * This enum class is simply used as "global value" holder for the given scores at the start of the program, so they
 * can be accessed from everywhere and don't have to be propagated inside composed functions
 */
public enum ScoreValues {
    MATCH_SCORE,
    MIS_MATCH_SCORE,
    GAP_PENALTY;

    private int value;
    private boolean isSet = false;

    /**
     * sets the values of the scores
     * @param value the value of which the score should be set
     */
    public void setValue(int value) {
        if (! isSet) {
            this.value = value;
            this.isSet = true;
        } else {
            throw new IllegalStateException("value already set for: " + this.name());
        }
    }

    /**
     * gets the values of the scores
     * @return scoreValue
     */
    public int value() {
        if (isSet) {
            return this.value;
        } else {
            throw new IllegalStateException("Value was not set for: " + this.name());
        }
    }
}
