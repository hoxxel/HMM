package main.hmm.profil;

import main.fastaparser.Sequence;

/**
 * Wrapper fuer Sequenz {@link Sequence}, die zusaetzlich das Ergebnis des Viterbi-Algo aus {@link ProfilHMM} haelt.
 *
 * @author Soeren Metje
 */
public class ViterbiPath {
    /**
     * Sequenz {@link Sequence}
     */
    private final Sequence sequence;
    /**
     * Bewertung
     */
    private final double score;
    /**
     * Zustands-Pfad
     */
    private final char[] statePath;

    /**
     * Konstruktor
     *
     * @param sequence  Sequenz
     * @param score     Bewertung
     * @param statePath Zustands-Pfad
     */
    public ViterbiPath(Sequence sequence, double score, char[] statePath) {
        this.sequence = sequence;
        this.score = score;
        this.statePath = statePath;
    }

    /**
     * liefert Sequenz {@link Sequence} zurueck
     *
     * @return Sequenz {@link Sequence}
     */
    public Sequence getSequence() {
        return sequence;
    }

    /**
     * liefert Bewertung zurueck
     *
     * @return Bewertung
     */
    public double getScore() {
        return score;
    }

    /**
     * liefert Zustands-Pfad zurueck
     *
     * @return Zustands-Pfad
     */
    public char[] getStatePath() {
        return statePath;
    }

    /**
     * liefert String mit Infos ueber den Zustands-Pfad zurueck
     *
     * @return String mit Infos ueber den Zustands-Pfad
     */
    @Override
    public String toString() {
        return String.valueOf(statePath);
    }
}
