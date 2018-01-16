package main.hmm.profil;

import main.fastaparser.Sequence;
import main.hmm.casino.CasinoHMM;
import main.logger.Log;

import java.util.List;

public class ProfilHMM {

    /*
    PROFIL HIDDEN MARKOV MODEL
               __              __                __                __                __                                  __                __
              (  )            (  )              (  )              (  )              (  )                                (  )              (  )
              /---\           /---\             /---\             /---\             /---\                               /---\             /---\
             /  I  \         /  I  \           /  I  \           /  I  \           /  I  \                             /  I  \           /  I  \
             \     /         \     /           \     /           \     /           \     /                             \     /           \     /
            />\---/\        />\---/\          />\---/\          />\---/\          />\---/\                            />\---/\          />\---/\
           /        \      /         \       /         \       /         \       /         \       /                           \       /         \
    |-----|         |------|          |------|          |------|          |------|          |------|                            |------|          |-----|
    |  B  |-------->|   M  |--------->|   M  |--------->|   M  |--------->|   M  |--------->|   M  |    ##   ##  ##   --------->|   M  |--------->|  E  |
    |-----|\        |------|\       />|------|\       />|------|\       />|------|\       />|------|    ##   ##  ##   \       />|------|        />|-----|
             \                \    /            \    /            \    /            \    /                              \    /                 /
               \                \/                \/                \/                \/                                  \/                 /
                 \     __     /    \     __     /    \     __     /    \     __     /    \     __                       /    \     __     /
                   \>/    \ /        \>/    \ /        \>/    \ /        \>/    \ /        \>/    \                   /        \>/    \ /
                    (  D   )          (  D   )          (  D   )          (  D   )          (  D   )                            (  D   )
                     \ __ /            \ __ /            \ __ /            \ __ /            \ __ /                              \ __ /

     */

    public static final int PSEUDO_COUNT_EMISSION = 1;
    public static final int PSEUDO_COUNT_TRANSITION = 1;
    public static final double THRESHOLD_MATCHSTATE = .5d; // min amount of nucleotides (no gap) to count column as match-state

    private static final char GAP = '-';
    private static final char[] BASES = {'A', 'C', 'G', 'U'};
    private static final char STATE_MATCH = 'M', STATE_INSERT = 'I', STATE_DELETE = 'D', STATE_END = 'E', STATE_IGNORE = ' ';
    private static final char[] STATES = {STATE_MATCH, STATE_INSERT, STATE_DELETE};
    private static final int STATE_COUNT = STATES.length;
    private static final Trans[] TRANSITIONS = {new Trans(STATE_MATCH, STATE_MATCH), new Trans(STATE_MATCH, STATE_INSERT), new Trans(STATE_MATCH, STATE_DELETE),
            new Trans(STATE_INSERT, STATE_MATCH), new Trans(STATE_INSERT, STATE_INSERT), new Trans(STATE_INSERT, STATE_DELETE),
            new Trans(STATE_DELETE, STATE_MATCH), new Trans(STATE_DELETE, STATE_INSERT), new Trans(STATE_DELETE, STATE_DELETE)};
    private static final Trans[] END_TRANSITIONS = {new Trans(STATE_MATCH, STATE_END), new Trans(STATE_INSERT, STATE_END), new Trans(STATE_DELETE, STATE_END)};


    private double[][] emissionProbMatch;
    private double[][] emissionProbInsert;
    private double[][][] transitionProb;
    private int lengthModel; // interpreting beginning-state als match-state

    public ProfilHMM() {
    }

    public ProfilHMM(List<Sequence> sequencesTrain) {
        buildModel(sequencesTrain);

        // convert into logspace (can be done before Viterbi-Algo is running)
        CasinoHMM.convertToLogspace(transitionProb);
        CasinoHMM.convertToLogspace(emissionProbMatch);
        CasinoHMM.convertToLogspace(emissionProbInsert);
    }

    public synchronized void buildModel(List<Sequence> sequencesTrain) {
        Log.iLine("Building ProfilHMM -----------------------------");

        // find Match or Insertion-States and Model length --------------------------------------------------------
        int length = sequencesTrain.get(0).getSequence().length();
        int[] gapCounts = new int[length];
        int[][] baseCounts = new int[length][BASES.length];

        for (Sequence seq : sequencesTrain) {
            for (int i = 0; i < length; i++) {
                char base = seq.getSequence().charAt(i);
                if (base == GAP) {
                    gapCounts[i]++;
                } else {
                    baseCounts[i][baseToIndex(base)]++;
                }
            }
        }

        boolean[] matchState = new boolean[length];
        lengthModel = 1;

        {
            int countThres = (int) (sequencesTrain.size() * (1d - THRESHOLD_MATCHSTATE));
            for (int i = 0; i < gapCounts.length; i++) {
                if (gapCounts[i] <= countThres) {
                    lengthModel++;
                    matchState[i] = true;
                }
            }
        }

        Log.iLine("Sequence length = " + length);
        Log.iLine("Model length = " + lengthModel + " (interpreting beginning-state als match-state)");


        // calc Emission Prob ---------------------------------------------------------------------------

        emissionProbMatch = new double[lengthModel][BASES.length]; // emission-probabilities for match-states
        emissionProbInsert = new double[lengthModel][BASES.length]; // emission-probabilities for insert-states


        int[] baseCountsInsert = new int[BASES.length];
        for (int i = 0, iModel = 1; i < length + 1; i++) {

            boolean end = i >= length;

            if (end || matchState[i]) {
                // set Insert-State Prob
                {
                    double[] emissionProbVector = emissionProbInsert[iModel - 1];

                    int sumEmissionCount = 0;

                    for (int baseCount : baseCountsInsert) {
                        sumEmissionCount += baseCount + PSEUDO_COUNT_EMISSION;
                    }

                    for (int j = 0; j < emissionProbVector.length; j++) {
                        emissionProbVector[j] = ((double) (PSEUDO_COUNT_EMISSION + baseCountsInsert[j])) / sumEmissionCount;
                    }
                }

            }

            if (!end) {
                if (matchState[i]) {
                    // set Match-State Prob
                    {
                        double[] emissionProbVector = emissionProbMatch[iModel];
                        int[] baseCountVector = baseCounts[i];

                        int sumEmissionCount = 0;

                        for (int baseCount : baseCountVector) {
                            sumEmissionCount += baseCount + PSEUDO_COUNT_EMISSION;
                        }

                        for (int j = 0; j < emissionProbVector.length; j++) {
                            emissionProbVector[j] = ((double) (PSEUDO_COUNT_EMISSION + baseCountVector[j])) / sumEmissionCount;
                        }
                    }

                    // reset Inset-State-Counts
                    for (int j = 0; j < baseCountsInsert.length; j++) {
                        baseCountsInsert[j] = 0;
                    }

                    iModel++;
                }
                // Insert-State
                else {
                    for (int j = 0; j < baseCountsInsert.length; j++) {
                        baseCountsInsert[j] += baseCounts[i][j];
                    }
                }
            }
        }


        // output
        Log.dLine("Gap-Counts (m=Match-State \u001B[32mI\u001B[0m=Insert-State):");
        for (int i = 0; i < length; i++) {
            Log.d(matchState[i] ? "   m " : "\u001B[32m   I \u001B[0m");
        }
        Log.dLine();

        for (int i = 0; i < gapCounts.length; i++) {
            Log.d(String.format("%4d ", gapCounts[i]));
        }
        Log.dLine();
        Log.dLine("Nucleotide-Counts: ");
        for (int i = 0; i < baseCounts[0].length; i++) {
            for (int j = 0; j < baseCounts.length; j++) {

                Log.d(String.format("%4d ", baseCounts[j][i]));
            }
            Log.dLine();
        }

        Log.dLine();
        Log.dLine("Emission Prob Match: (Pseudo-Count = " + PSEUDO_COUNT_EMISSION + ")");
        for (int i = 0; i < emissionProbMatch[0].length; i++) {
            for (int j = 0; j < emissionProbMatch.length; j++) {

                Log.d(String.format("%.2f ", emissionProbMatch[j][i]));
            }
            Log.dLine();
        }

        Log.dLine();
        Log.dLine("Emission Prob Insert: (Pseudo-Count = " + PSEUDO_COUNT_EMISSION + ")");
        for (int i = 0; i < emissionProbInsert[0].length; i++) {
            for (int j = 0; j < emissionProbInsert.length; j++) {

                Log.d(String.format("%.2f ", emissionProbInsert[j][i]));
            }
            Log.dLine();
        }

        Log.dLine();

        // calc Transition Count ----------------------------------------------------------------------
        Log.dLine("States and Transition Counts:");

        int[][][] transitionCount = new int[STATES.length][STATES.length][lengthModel];
        int[] endTransitionCount = new int[END_TRANSITIONS.length];


        Log.d("    ");
        for (int i = 0; i < length; i++) {
            Log.d(matchState[i] ? "  m" : "\u001B[32m  I\u001B[0m");
        }
        Log.dLine();
        for (int i1 = 0; i1 < sequencesTrain.size(); i1++) {
            Sequence sequence = sequencesTrain.get(i1);
            String sequenceString = sequence.getSequence();

            // output sequence
            Log.dLine("\u001B[37m");
            Log.d(String.format("%4s", sequence.getDescription()));
            for (int i = 0; i < length; i++) {
                Log.d("  " + sequenceString.charAt(i));
            }
            Log.dLine("\u001B[0m");


            // get States and Transitions
            char lastState = STATE_MATCH;
            int insertCount = 0;

            Log.d("    ");
            for (int i = 0, iModel = 0; i < length + 1; i++) {
                char state = getState(sequenceString, matchState, i);
                if (state != STATE_IGNORE) {
                    if (state == STATE_END) {
                        endTransitionCount[endTransitionToIndex(lastState, state)]++;
                    } else {

                        if (state == STATE_INSERT && lastState == STATE_INSERT) {
                            insertCount++;
                        } else if (state == STATE_INSERT) {
                            transitionCount[stateToIndex(lastState)][stateToIndex(state)][iModel]++;
                        }
                        // End of InsertStates
                        else {
                            if (lastState == STATE_INSERT) { // set Insert-Insert
                                transitionCount[stateToIndex(STATE_INSERT)][stateToIndex(lastState)][iModel] += insertCount;
                                insertCount = 0;
                            }

                            transitionCount[stateToIndex(lastState)][stateToIndex(state)][iModel]++;

                            iModel++;
                        }
                    }

                    lastState = state;
                }
                Log.d("  " + state);
            }
            Log.dLine();
        }


        // output Transition Counts
        Log.dLine();
        /*
        for (int i = 0; i < initTransitionCount.length; i++) {
            Log.dLine(INIT_TRANSITIONS[i] + ":  " + String.format("%4d ", initTransitionCount[i]));
        }
        */
        for (int i = 0; i < STATES.length; i++) {
            for (int j = 0; j < STATES.length; j++) {
                Log.d(STATES[i] + "" + STATES[j] + ":  ");
                for (int k = 0; k < transitionCount[i][j].length; k++) {
                    Log.d(String.format("%4d ", transitionCount[i][j][k]));
                }
                Log.dLine();
            }
        }
        for (int i = 0; i < endTransitionCount.length; i++) {
            Log.dLine(END_TRANSITIONS[i] + ":  " + String.format("%4d ", endTransitionCount[i]));
        }


        // calc Transition Prob ----------------------------------------------------------------------
        /*
        initTransitionProb = new double[INIT_TRANSITIONS.length];
        {
            int sumInitTrans = 0;
            for (int anInitTransitionCount : initTransitionCount) {
                sumInitTrans += anInitTransitionCount + PSEUDO_COUNT_TRANSITION;
            }

            for (int i = 0; i < initTransitionCount.length; i++) {
                initTransitionProb[i] = ((double) initTransitionCount[i] + PSEUDO_COUNT_TRANSITION) / sumInitTrans;
            }
        }
        */

        transitionProb = new double[STATES.length][STATES.length][lengthModel];
        for (int i = 0; i < lengthModel; i++) {

            int[] sumTrans = new int[STATES.length]; // vergleichsgroessen
            for (int j = 0; j < STATES.length; j++) {
                for (int k = 0; k < STATES.length; k++) {

                    sumTrans[j] += transitionCount[j][k][i] + PSEUDO_COUNT_TRANSITION;
                }
            }

            for (int j = 0; j < STATES.length; j++) {
                for (int s = 0; s < STATES.length; s++) {
                    transitionProb[j][s][i] = ((double) transitionCount[j][s][i] + PSEUDO_COUNT_TRANSITION) / sumTrans[j];
                }
            }
        }


        // output Transition Prob
        Log.iLine();
        Log.iLine("Transition Prob: (Pseudo-Count = " + PSEUDO_COUNT_TRANSITION + ")");

        /*
        for (int i = 0; i < initTransitionProb.length; i++) {
            Log.iLine(String.format("%s:  %.2f", INIT_TRANSITIONS[i].toString(), initTransitionProb[i]));
        }
        */
        for (int i = 0; i < transitionProb.length; i++) {
            double[][] transProbMatrix = transitionProb[i];
            for (int j = 0; j < transProbMatrix.length; j++) {
                double[] transProbVector = transProbMatrix[j];
                Log.i(String.format("%s:  ", STATES[i] + "" + STATES[j]));
                for (double prob : transProbVector) {
                    Log.i(String.format("%.2f ", prob));
                }
                Log.iLine();
            }
        }
    }


    /**
     * Implementation des Viterbi-Algorithmus für den logarithmischen Raum
     *
     * @param observations Beobachtungsfolge
     * @return Zustands-Pfad
     */
    public synchronized char[] viterbi(final char[] observations) { // FIXME
        Log.iLine("Viterbi ProfilHMM -----------------------------");

        // init
        int[] observationIndices = basesToIndices(observations);
        int length = observationIndices.length + 1;

        double[][][] viterbiVar = new double[STATE_COUNT][length][lengthModel];
        int[][][] viterbiArg = new int[STATE_COUNT][length][lengthModel];

        // init
        for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {
            char state = STATES[stateIndex];
            if (state == STATE_MATCH) {
                viterbiVar[stateIndex][0][0] = 0;
                for (int j = 1; j < lengthModel; j++) {
                    viterbiVar[stateIndex][0][j] = Double.NEGATIVE_INFINITY;
                }
                for (int i = 1; i < length; i++) {
                    viterbiVar[stateIndex][i][0] = Double.NEGATIVE_INFINITY;
                }
            } else if (state == STATE_INSERT) {
                for (int j = 0; j < lengthModel; j++) {
                    viterbiVar[stateIndex][0][j] = Double.NEGATIVE_INFINITY;
                }
            } else if (state == STATE_DELETE) {
                for (int i = 0; i < length; i++) {
                    viterbiVar[stateIndex][i][0] = Double.NEGATIVE_INFINITY;
                }
            }
        }


        Log.iLine("ViterbiVar Init: "); //TODO remove
        // [STATE_COUNT][length][lengthModel]
        for (int j = 0; j < lengthModel; j++) {
            for (int k = 0; k < STATES.length; k++) {
                Log.i((k == 0 ? String.format("j%3d%s ", j, STATES[k]) : "    " + STATES[k] + " "));
                for (int i = 0; i < length; i++) {
                    Log.i(String.format("%.4s ", String.format("%f", viterbiVar[k][i][j])));
                }
                Log.iLine();
            }
            Log.iLine();
        }
        Log.iLine();

        // iterate observations indices
        for (int i = 1; i < length; i++) {
            // iterate model indices
            for (int j = 1; j < lengthModel; j++) {
                // iterate states indices
                for (int s = 0; s < STATE_COUNT; s++) {

                    char state = STATES[s];

                    double emissionProb = 0d; // 0 is neutral element of addition (log-space)

                    int nextI = i;
                    int nextJ = j;

                    if (state == STATE_MATCH) {
                        nextI = i - 1;
                        nextJ = j - 1;
                        emissionProb = emissionProbMatch[j][observationIndices[i - 1]];
                    } else if (state == STATE_INSERT) {
                        nextI = i - 1;
                        emissionProb = emissionProbInsert[j][observationIndices[i - 1]];
                    } else if (state == STATE_DELETE) {
                        nextJ = j - 1;
                    }

                    //find max
                    double maxProb = Double.NEGATIVE_INFINITY;
                    int maxArg = -1; // maximizing argument

                    for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                        double prob = viterbiVar[stateIndex][nextI][nextJ] + transitionProb[s][stateIndex][nextJ]; // log-space
                        if (prob > maxProb) {
                            maxProb = prob;
                            maxArg = stateIndex;
                        }
                    }

                    viterbiVar[s][i][j] = emissionProb + maxProb;
                    viterbiArg[s][i][j] = maxArg;
                }
            }
        }

        Log.iLine("ViterbiVar: "); // TODO change to debug
        // [STATE_COUNT][length][lengthModel]
        for (int j = 0; j < lengthModel; j++) {
            for (int k = 0; k < STATES.length; k++) {
                Log.i((k == 0 ? String.format("j%3d%s ", j, STATES[k]) : "    " + STATES[k] + " "));
                for (int i = 0; i < length; i++) {
                    Log.i(String.format("%.4s ", String.format("%f", viterbiVar[k][i][j])));
                }
                Log.iLine();
            }
            Log.iLine();
        }
        Log.iLine();

        // backtrace init
        int zLast = -1;
        {
            double probLast = Double.NEGATIVE_INFINITY;
            for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                double prob = viterbiVar[stateIndex][length - 1][lengthModel - 1];
                if (prob > probLast) {
                    zLast = stateIndex;
                    probLast = prob;
                }
            }
        }

        int[] stateIndexPath = new int[length]; // z
        char[] statePath = new char[length]; // x

        stateIndexPath[length - 1] = zLast;
        statePath[length - 1] = STATES[zLast];

        // backtrace iterate // FIXME
        {
            int i = length - 1, j = lengthModel - 1;
            while (i > 0 && j > 0) {
                int stateIndex = viterbiArg[stateIndexPath[i]][i][j];
                char state = STATES[stateIndex];

                stateIndexPath[i - 1] = stateIndex;
                statePath[i - 1] = state;
                if (state == STATE_MATCH) {
                    i--;
                    j--;
                } else if (state == STATE_INSERT) {
                    i--;
                } else if (state == STATE_DELETE) {
                    j--;
                }
            }
        }

        /*
        for (int i = length - 1; i > 0; i--) {
            for (int j = lengthModel - 1; j > 0; j--) {
                int m = viterbiArg[stateIndexPath[i]][i][j];
                stateIndexPath[i - 1] = m;
                statePath[i - 1] = STATES[m];
            }
        }
        */

        for (int i = 0; i < stateIndexPath.length; i++) {
            Log.d(stateIndexPath[i]);
        }
        Log.dLine();

        return statePath;
    }


    /**
     * mappt Beaobachtung-Folge auf entsprechende Index-Folge
     *
     * @param observations Beobachtungs-Folge
     * @return entsprechende Index-Folge
     */
    private static int[] basesToIndices(final char[] observations) {
        int length = observations.length;
        int[] ret = new int[length];

        for (int i = 0; i < length; i++) {
            ret[i] = baseToIndex(observations[i]);
        }

        return ret;
    }


    /**
     * gibt den Zustand zurueck, in dem sich das HMM an uebergebenem index in uebergebener Sequenz befindet
     * oder ein Leerzeichen, falls sich das Modell an der Stelle im Insert-Zustand befindet, aber kein Zeichen in der Sequenz vorhanden ist.
     * Das uebergebene Feld matchState muss auskunft darueber geben, ob sich das Modell an uebergebenem index im Insert- oder Match-Zustand befindet (true falls Match-Zustand).
     *
     * @param seq        Sequenz
     * @param matchState Insert- oder Match-Zustaende
     * @param index      position Sequenz
     * @return Zustand
     */
    private static char getState(String seq, boolean[] matchState, int index) {

        if (index >= seq.length())
            return STATE_END;

        boolean match = matchState[index];

        if (!match) {
            if (seq.charAt(index) != GAP) {
                return STATE_INSERT;
            } else {
                return STATE_IGNORE;
            }
        }

        // at i is a Match-state
        if (seq.charAt(index) == GAP)
            return STATE_DELETE;
        return STATE_MATCH;
    }

    private static int transitionToIndex(char start, char end) {
        return transitionToIndex(TRANSITIONS, start, end);

    }

    private static int endTransitionToIndex(char start, char end) {
        return transitionToIndex(END_TRANSITIONS, start, end);
    }

    private static int transitionToIndex(Trans[] array, char start, char end) {
        for (int i = 0, arrayLength = array.length; i < arrayLength; i++) {
            Trans trans = array[i];
            if (trans.equals(start, end))
                return i;
        }
        throw new IllegalArgumentException("Transition " + start + "" + end + " not found");
    }


    private static int baseToIndex(char base) {
        return charToIndex(BASES, base);
    }

    private static int stateToIndex(char state) {
        return charToIndex(STATES, state);
    }

    private static int charToIndex(char[] array, char c) {
        for (int i = 0, basesLength = array.length; i < basesLength; i++) {
            if (array[i] == c)
                return i;
        }
        throw new IllegalArgumentException("Character " + c + " not found");
    }
}
