package main.hmm.profil.viterbi;

import main.fastaparser.Sequence;
import main.hmm.profil.RNAProfilHMM;

/**
 * Wrapper fuer Sequenz {@link Sequence}, die zusaetzlich das Ergebnis des Viterbi-Algo aus {@link RNAProfilHMM} haelt.
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
     * Liefert Sequenz {@link Sequence} zurueck
     *
     * @return Sequenz {@link Sequence}
     */
    public Sequence getSequence() {
        return sequence;
    }

    /**
     * Liefert Bewertung zurueck
     *
     * @return Bewertung
     */
    public double getScore() {
        return score;
    }

    /**
     * Liefert Zustands-Pfad zurueck
     *
     * @return Zustands-Pfad
     */
    public char[] getStatePath() {
        return statePath;
    }

    /**
     * Liefert String mit Infos ueber den Zustands-Pfad zurueck
     *
     * @return String mit Infos ueber den Zustands-Pfad
     */
    @Override
    public String toString() {
        return "path " + sequence.getDescription();
    }
}
