package main.hmm.profil;

import main.fastaparser.Sequence;

public class ViterbiPath {
    private final Sequence sequence;
    private final double score;
    private final char[] statePath;

    public ViterbiPath(Sequence sequence, double score, char[] statePath) {
        this.sequence = sequence;
        this.score = score;
        this.statePath = statePath;
    }

    public Sequence getSequence() {
        return sequence;
    }

    public double getScore() {
        return score;
    }

    public char[] getStatePath() {
        return statePath;
    }

    @Override
    public String toString() {
        return String.valueOf(statePath);
    }
}
