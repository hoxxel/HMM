package main.hmm.casino;

/**
 * <p>
 * Hidden Markov Model fuer ein unehrliches Casino. F = Fair, L = Loaded.
 * Ruft Konstruktor {@link HMM} mit Konstanten auf.
 * </p>
 *
 * @author Soeren Metje
 */
public class CasinoHMM extends HMM {

    /**
     * beobachtbare Ereignisse
     */
    private static final char[] OBSERVATION_SPACE = {'1', '2', '3', '4', '5', '6'};

    /**
     * Zustaende in denen sich das Modell befindet
     */
    private static final char[] STATE_CHAR = {'F', 'L'};

    /**
     * Uebergangswahrscheinlichen aus dem Startzustand in die jeweiligen Zustaende
     */
    private static final double[] INIT_PROBABILITIES = new double[]{.5d, .5d};

    /**
     * Uebergangswahrscheinlichen zwischen den Zustaenden
     */
    private static final double[][] TRANSITION_MATRIX = new double[][]{{.95d, .05d}, {.1d, .9d}};

    /**
     * Beobachtungswahrscheinlichketen der Ereignisse in den jeweiligen Zustaenden
     */
    private static final double[][] EMISSION_MATRIX = new double[][]{{1d / 6, 1d / 6, 1d / 6, 1d / 6, 1d / 6, 1d / 6}, {.1d, .1d, .1d, .1d, .1d, .5d}};

    /**
     * Konstruktor.
     * Konvertiert die Matrizen fuer die Uebergangswahrscheinlichen und Beobachtungswahrscheinlichketen in den logarithmischen Raum
     */
    public CasinoHMM() {
        super(OBSERVATION_SPACE, STATE_CHAR, INIT_PROBABILITIES, TRANSITION_MATRIX, EMISSION_MATRIX);
    }
}
