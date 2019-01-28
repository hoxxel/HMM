package main.hmm.casino;

import main.hmm.HMMFunc;

/**
 * <p>
 * Hidden Markov Model.
 * Anhand  https://en.wikipedia.org/wiki/Viterbi_algorithm  implementiert.
 * </p>
 * <p>
 * Enthaelt die Implementation des Viterbi-Algorithmus.
 * Dieser generiert aus einer uebergebenen Sequenz einen Zustands-Pfad.
 * </p>
 *
 * @author Soeren Metje
 */
public class HMM {

    /**
     * beobachtbare Ereignisse
     */
    protected final char[] observationSpace;

    /**
     * Zustaende in denen sich das Modell befindet
     */
    protected final char[] stateChar;

    /**
     * Anzahl der Zustaende
     */
    protected final int stateCount;

    /**
     * Uebergangswahrscheinlichen aus dem Startzustand in die jeweiligen Zustaende
     */
    protected final double[] initProbabilities;

    /**
     * Uebergangswahrscheinlichen zwischen den Zustaenden
     */
    protected final double[][] transitionMatrix;

    /**
     * Beobachtungswahrscheinlichketen der Ereignisse in den jeweiligen Zustaenden
     */
    protected final double[][] emissionMatrix;

    /**
     * Konstruktor.
     * Konvertiert die Matrizen fuer die Uebergangswahrscheinlichen und Beobachtungswahrscheinlichketen in den logarithmischen Raum
     *
     * @param observationSpace
     * @param stateChar
     * @param initProbabilities
     * @param transitionMatrix
     * @param emissionMatrix
     */
    public HMM(char[] observationSpace, char[] stateChar, double[] initProbabilities, double[][] transitionMatrix, double[][] emissionMatrix) {
        this.observationSpace = observationSpace;
        this.stateChar = stateChar;
        this.stateCount = stateChar.length;
        this.initProbabilities = initProbabilities;
        this.transitionMatrix = transitionMatrix;
        this.emissionMatrix = emissionMatrix;

        // calc log for each element in all matrices (can be done before Viterbi-Algo is running)
        HMMFunc.logspace(initProbabilities);
        HMMFunc.logspace(transitionMatrix);
        HMMFunc.logspace(emissionMatrix);
    }

    /**
     * Implementation des Viterbi-Algorithmus fuer bereits logarithmierte Werte.
     * Liefert den wahrscheinlichsten Zustands-Pfad bei uebergebenen Beobachtungen zurueck.
     *
     * @param observations Beobachtungsfolge
     * @return Zustands-Pfad
     * @throws IllegalArgumentException falls uebergebenes Feld == null oder
     *                                  Beobachtung nicht im Feld gefunden wird
     */
    public char[] viterbi(final char[] observations) throws IllegalArgumentException {
        if (observations == null)
            throw new IllegalArgumentException("observations is null");

        // init fields
        int[] observationIndices = observationsToIndices(observations);
        int length = observationIndices.length;

        double[][] viterbiVar = new double[stateCount][length];
        int[][] viterbiArg = new int[stateCount][length];

        // init
        for (int stateIndex = 0; stateIndex < stateCount; stateIndex++) {

            viterbiVar[stateIndex][0] = initProbabilities[stateIndex] + emissionMatrix[stateIndex][observationIndices[0]]; // log-space
            viterbiArg[stateIndex][0] = -1;
        }

        // iterate observations indices
        for (int i = 1; i < length; i++) {
            // iterate states indices
            for (int j = 0; j < stateCount; j++) {


                //find max
                double maxProb = Double.NEGATIVE_INFINITY;
                int maxArg = -1; // maximizing argument

                for (int stateIndex = 0; stateIndex < stateCount; stateIndex++) {

                    double prob = viterbiVar[stateIndex][i - 1] + transitionMatrix[stateIndex][j] + emissionMatrix[j][observationIndices[i]]; // log-space
                    if (prob > maxProb) {
                        maxProb = prob;
                        maxArg = stateIndex;
                    }
                }

                viterbiVar[j][i] = maxProb;
                viterbiArg[j][i] = maxArg;
            }
        }

        // backtrace init
        int zLast = -1;
        {
            double probLast = Double.NEGATIVE_INFINITY;
            for (int stateIndex = 0; stateIndex < stateCount; stateIndex++) {

                double prob = viterbiVar[stateIndex][length - 1];
                if (prob > probLast) {
                    zLast = stateIndex;
                    probLast = prob;
                }
            }
        }

        int[] z = new int[length]; // stateIndexPath
        char[] x = new char[length]; // statePath

        z[length - 1] = zLast;
        x[length - 1] = stateChar[zLast];

        // backtrace iterate
        for (int i = length - 1; i > 0; i--) {
            int m = viterbiArg[z[i]][i];
            z[i - 1] = m;
            x[i - 1] = stateChar[m];
        }

        return x;
    }

    /**
     * mappt Beaobachtung-Folge auf entsprechende Index-Folge
     *
     * @param observations Beobachtungs-Folge
     * @return entsprechende Index-Folge
     */
    private int[] observationsToIndices(final char[] observations) {
        return HMMFunc.charsToIndices(observationSpace, observations);
    }

    /**
     * mappt Beaobachtung auf den entsprechenden Index
     *
     * @param observation Beobachtung
     * @return entsprechender Index
     */
    private int obesrvationToIndex(final char observation) {
        return HMMFunc.charToIndex(observationSpace, observation);
    }
}
