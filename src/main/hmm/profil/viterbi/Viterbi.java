package main.hmm.profil.viterbi;

import main.fastaparser.Sequence;
import main.hmm.profil.ProfilHMM;
import main.logger.Log;

import java.util.LinkedList;
import java.util.List;

/**
 * Enthaelt die Implementation des Viterbi-Algorithmus fuer logarithmische Werte.
 * Dieser generiert aus einer uebergebenen Sequenz einen Zustands-Pfad. M = Match, D = Delete, I = Insert.
 *
 * @author Soeren Metje
 */
public class Viterbi {

    /**
     * Implementation des Viterbi-Algorithmus fuer bereits logarithmierte Werte.
     * Liefert den wahrscheinlichsten Zustands-Pfad mit score bei uebergebenen Beobachtungen zurueck.
     *
     * @param model    Profil Hidden Markov Model
     * @param sequence Beobachtungsfolge
     * @return Zustands-Pfad
     * @throws IllegalArgumentException falls uebergebene Sequenz {@link Sequence} == null
     *                                  oder falls Beobachtung nicht im Feld entsprechenden gefunden wird
     */
    public static ViterbiPath viterbi(final ProfilHMM model, final Sequence sequence) throws IllegalArgumentException {
        if (sequence == null)
            throw new IllegalArgumentException("sequence is null");

        // init
        char[] observations = sequence.getNucleotideSequence().toCharArray();
        int[] observationIndices = model.observationsToIndices(observations);
        int length = observationIndices.length + 1;

        // FILL MATRIX ----------------------------------------------------------------------------------
        int lengthModel = model.getLengthModel();
        double[][][] viterbiVar = new double[ProfilHMM.STATE_COUNT][length][lengthModel];
        int[][][] viterbiArg = new int[ProfilHMM.STATE_COUNT][length][lengthModel];

        // init
        {
            double initValue = Double.NEGATIVE_INFINITY;

            int stateIndex = ProfilHMM.STATE_MATCH_INDEX;
            viterbiVar[stateIndex][0][0] = 0d;
            viterbiArg[stateIndex][0][0] = -1;
            for (int j = 1; j < lengthModel; j++) {
                viterbiVar[stateIndex][0][j] = initValue;
            }
            for (int i = 1; i < length; i++) {
                viterbiVar[stateIndex][i][0] = initValue;
            }

            stateIndex = ProfilHMM.STATE_INSERT_INDEX;
            for (int j = 0; j < lengthModel; j++) {
                viterbiVar[stateIndex][0][j] = initValue;
            }

            stateIndex = ProfilHMM.STATE_DELETE_INDEX;
            for (int i = 0; i < length; i++) {
                viterbiVar[stateIndex][i][0] = initValue;
            }
        }

        // iterate observations indices
        for (int i = 0; i < length; i++) {
            // iterate model indices
            for (int j = 0; j < lengthModel; j++) {
                // iterate states indices
                for (int s = 0; s < ProfilHMM.STATE_COUNT; s++) { // order of iteration-loops is relevant!
                    char state = ProfilHMM.STATES[s];

                    int iShift = i, jShift = j;
                    double[][] emissionProbMatrix = null;
                    if (state == ProfilHMM.STATE_MATCH) {
                        iShift -= 1;
                        jShift -= 1;
                        emissionProbMatrix = model.getEmissionProbMatch();
                    } else if (state == ProfilHMM.STATE_INSERT) {
                        iShift -= 1;
                        emissionProbMatrix = model.getEmissionProbInsert();
                    } else if (state == ProfilHMM.STATE_DELETE) {
                        jShift -= 1;
                    } else
                        throw new RuntimeException("no valid state");

                    // calc must be possible for Delete-State in first column ans Insert-State in first row
                    if (iShift >= 0 && jShift >= 0) {
                        //find max
                        double maxProb = Double.NEGATIVE_INFINITY;
                        int maxArg = -1; // maximizing argument

                        for (int stateIndex = 0; stateIndex < ProfilHMM.STATE_COUNT; stateIndex++) {

                            double[][][] transitionProb = model.getTransitionProb();
                            double prob = viterbiVar[stateIndex][iShift][jShift] + transitionProb[stateIndex][s][jShift]; // log-space
                            if (prob > maxProb) {
                                maxProb = prob;
                                maxArg = stateIndex;
                            }
                        }

                        double emissionProb = 0d; // 0 is neutral element of addition (log-space)

                        if (emissionProbMatrix != null) {
                            emissionProb = emissionProbMatrix[j][observationIndices[i - 1]];
                        }
                        viterbiVar[s][i][j] = emissionProb + maxProb;
                        viterbiArg[s][i][j] = maxArg;
                    }
                }
            }
        }

        if (Log.isPrintDebug()) {
            // Debug output viterbi 3d-matrix
            StringBuilder outViterbiVar = new StringBuilder("ViterbiVar: \n");
            // [STATE_COUNT][length][lengthModel]
            for (int j = 0; j < lengthModel; j++) {
                for (int k = 0; k < ProfilHMM.STATES.length; k++) {
                    outViterbiVar.append("\u001B[37m").append(k == 0 ? String.format("j%3d%s ", j, ProfilHMM.STATES[k]) : "    " + ProfilHMM.STATES[k] + " ").append("\u001B[0m");
                    for (int i = 0; i < length; i++) {
                        outViterbiVar.append(String.format("%.5s ", String.format("%f", viterbiVar[k][i][j])));
                    }
                    outViterbiVar.append('\n');
                }
                outViterbiVar.append('\n');
            }
            Log.dLine(outViterbiVar.toString());

            StringBuilder outViterbiArg = new StringBuilder("ViterbiArg: \n");
            // [STATE_COUNT][length][lengthModel]
            for (int j = 0; j < lengthModel; j++) {
                for (int k = 0; k < ProfilHMM.STATES.length; k++) {
                    outViterbiArg.append("\u001B[37m").append(k == 0 ? String.format("j%3d%s ", j, ProfilHMM.STATES[k]) : "    " + ProfilHMM.STATES[k] + " ").append("\u001B[0m");
                    for (int i = 0; i < length; i++) {
                        int a = viterbiArg[k][i][j];
                        outViterbiArg.append(String.format("%5s ", (a >= 0 ? String.valueOf(ProfilHMM.STATES[a]) : a)));
                    }
                    outViterbiArg.append('\n');
                }
                outViterbiArg.append('\n');
            }
            Log.dLine(outViterbiArg.toString());
        }

        // BACKTRACE -------------------------------------------------------------------------------------
        List<Character> listStatePath = new LinkedList<>();
        double score = Double.NEGATIVE_INFINITY;
        {
            // backtrace init / Find path with max prob
            int i = length - 1, j = lengthModel - 1;
            int stateIndexEnd = -1;
            {
                for (int stateIndex = 0; stateIndex < ProfilHMM.STATE_COUNT; stateIndex++) {

                    double[][][] transitionProb = model.getTransitionProb();
                    double prob = viterbiVar[stateIndex][i][j] + transitionProb[stateIndex][ProfilHMM.STATE_MATCH_INDEX][j]; // log-space
                    if (prob > score) {
                        stateIndexEnd = stateIndex;
                        score = prob;
                    }
                }
            }
            viterbiVar = null; // no reference left -> allow GC to trash

            listStatePath.add(ProfilHMM.STATES[stateIndexEnd]);

            // backtrace iterate
            try {
                while (i >= 0 && j >= 0 && (i > 1 || j > 1)) { // FIXME correct?!
                    int stateIndex = viterbiArg[stateIndexEnd][i][j];
                    char state = ProfilHMM.STATES[stateIndex];

                    listStatePath.add(0, state);

                    if (state == ProfilHMM.STATE_MATCH) {
                        i--;
                        j--;
                    } else if (state == ProfilHMM.STATE_INSERT) {
                        i--;
                    } else if (state == ProfilHMM.STATE_DELETE) {
                        j--;
                    }
                }

            } catch (ArrayIndexOutOfBoundsException e) {
                throw new ArrayIndexOutOfBoundsException(e.getMessage() + " i=" + i + " j=" + j + " seq=" + sequence.getDescription());
            }
        }

        // copy into char array
        char[] statePath = new char[listStatePath.size()];
        {
            int i = 0;
            for (Character c : listStatePath) {
                statePath[i] = c;
                i++;
            }
        }

        return new ViterbiPath(sequence, score, statePath);
    }
}
