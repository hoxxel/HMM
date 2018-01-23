package main.hmm.profil;

import main.fastaparser.Sequence;

import java.util.Arrays;

public class ViterbiPath implements Comparable<ViterbiPath> {
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

    /**
     * liefert true zurueck falls die Werte des uebergebenen Objekts mit dem aufrufenden identisch sind. Ansonsten false
     *
     * @param o zu vergleichendes Objekt
     * @return true zurueck falls die Werte des uebergebenen Objekts mit dem aufrufenden identisch sind. Ansonsten false
     */
    @Override
    public boolean equals(Object o) {
        if (o.getClass() != this.getClass())
            return false;
        ViterbiPath viterbiPath = (ViterbiPath) o;
        if (viterbiPath.score != this.score || !Arrays.equals(viterbiPath.statePath, this.statePath) || !viterbiPath.sequence.equals(this.sequence))
            return false;
        return true;
    }

    /**
     * Liefert
     * 1, falls der score des uebergebenen Objekts größer,
     * -1, falls der score des uebergebenen Objekts kleiner,
     * oder 0 zurueck, falls der score des uebergebenen Objekts gleich
     * dem des aufrufenden Objekts ist.
     * Wird für die sortierung anhand des scores benoetigt
     *
     * @param viterbiPath zu vergleichendes {@link ViterbiPath}-Objekt
     * @return 1, falls der score des uebergebenen Objekts größer,
     * -1, falls der score des uebergebenen Objekts kleiner,
     * oder 0 zurueck, falls der score des uebergebenen Objekts gleich
     * dem des aufrufenden Objekts ist.
     */
    @Override
    public int compareTo(ViterbiPath viterbiPath) {
        double c = score - viterbiPath.getScore();
        if (c > 0)
            return 1;
        if (c < 0)
            return -1;
        return 0;
    }

    @Override
    public String toString() {
        return String.valueOf(statePath);
    }
}
