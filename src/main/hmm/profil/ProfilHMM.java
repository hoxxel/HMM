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
    private static final char STATE_MATCH = 'M', STATE_INSERT = 'I', STATE_DELETE = 'D', STATE_BEGIN = 'B', STATE_END = 'E', STATE_IGNORE = ' ';
    private static final char[] STATES = {STATE_MATCH, STATE_INSERT, STATE_DELETE};
    private static final int STATE_COUNT = STATES.length;
    private static final Trans[] TRANSITIONS = {new Trans(STATE_MATCH, STATE_MATCH), new Trans(STATE_MATCH, STATE_INSERT), new Trans(STATE_MATCH, STATE_DELETE),
            new Trans(STATE_INSERT, STATE_MATCH), new Trans(STATE_INSERT, STATE_INSERT), new Trans(STATE_INSERT, STATE_DELETE),
            new Trans(STATE_DELETE, STATE_MATCH), new Trans(STATE_DELETE, STATE_INSERT), new Trans(STATE_DELETE, STATE_DELETE)};
    private static final Trans[] INIT_TRANSITIONS = {new Trans(STATE_BEGIN, STATE_MATCH), new Trans(STATE_BEGIN, STATE_INSERT), new Trans(STATE_BEGIN, STATE_DELETE)};
    private static final Trans[] END_TRANSITIONS = {new Trans(STATE_MATCH, STATE_END), new Trans(STATE_INSERT, STATE_END), new Trans(STATE_DELETE, STATE_END)};


    private double[][] emissionProbMatch;
    private double[][] emissionProbInsert;
    private double[] initTransitionProb;
    private double[][] transitionProb;

    public ProfilHMM() {
    }

    public ProfilHMM(List<Sequence> sequencesTrain) {
        buildModel(sequencesTrain);

        // convert into logspace (can be done before Viterbi-Algo is running)
        CasinoHMM.convertToLogspace(initTransitionProb);
        CasinoHMM.convertToLogspace(transitionProb);
        CasinoHMM.convertToLogspace(emissionProbMatch);
        CasinoHMM.convertToLogspace(emissionProbInsert);
    }

    public void buildModel(List<Sequence> sequencesTrain) {

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
        int lengthModel = 0;

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
        Log.iLine("Model length = " + lengthModel);


        // calc Emission Prob ---------------------------------------------------------------------------

        emissionProbMatch = new double[lengthModel + 1][BASES.length]; // emission-probabilities for match-states
        emissionProbInsert = new double[lengthModel + 1][BASES.length]; // emission-probabilities for insert-states


        int[] baseCountsInsert = new int[BASES.length];
        for (int i = 0, iModel = 0; i < length + 1; i++) { // FIXME distinguish between match and insert

            boolean end = i >= length;

            if (end || matchState[i]) {
                // set Insert-State Prob
                {
                    double[] emissionProbVector = emissionProbInsert[iModel];

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

        int[][] transitionCount = new int[TRANSITIONS.length][lengthModel + 1];
        int[] initTransitionCount = new int[INIT_TRANSITIONS.length];
        int[] endTransitionCount = new int[END_TRANSITIONS.length];


        Log.d("    ");
        for (int i = 0; i < length; i++) {
            Log.d(matchState[i] ? "  m" : "\u001B[32m  I\u001B[0m");
        }
        Log.dLine();
        for (int i1 = 0; i1 < 10; i1++) {
            Sequence sequence = sequencesTrain.get(i1);
            String sequenceString = sequence.getSequence();

            // output sequence
            Log.dLine("\u001B[37m");
            Log.d(String.format("%5s", sequence.getDescription()));
            for (int i = 0; i < length; i++) {
                Log.d("  " + sequenceString.charAt(i));
            }
            Log.dLine("\u001B[0m");


            // get States and Transitions
            char lastState = STATE_BEGIN;
            int insertCount = 0;

            Log.d("    ");
            for (int i = 0, iModel = 0; i < length + 1; i++) {
                char state = getState(sequenceString, matchState, i);
                if (state != STATE_IGNORE) {
                    if (state == STATE_END) {
                        endTransitionCount[endTransitionToIndex(lastState, state)]++;
                    } else if (lastState == STATE_BEGIN) {
                        initTransitionCount[initTransitionToIndex(lastState, state)]++;
                        if (state == STATE_MATCH || state == STATE_DELETE) // no Insert-State in first column
                            iModel++;
                    } else {

                        if (state == STATE_INSERT && lastState == STATE_INSERT) {
                            insertCount++;
                        } else if (state == STATE_INSERT) {
                            transitionCount[transitionToIndex(lastState, state)][iModel]++;
                        }
                        // End of InsertStates
                        else {
                            if (lastState == STATE_INSERT) { // set Insert-Insert
                                transitionCount[transitionToIndex(STATE_INSERT, lastState)][iModel] += insertCount;
                                insertCount = 0;
                            }

                            transitionCount[transitionToIndex(lastState, state)][iModel]++;

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
        for (int i = 0; i < initTransitionCount.length; i++) {
            Log.dLine(INIT_TRANSITIONS[i] + ":  " + String.format("%4d ", initTransitionCount[i]));
        }
        for (int i = 0; i < transitionCount.length; i++) {
            Log.d(TRANSITIONS[i] + ":  ");
            for (int j = 0; j < transitionCount[i].length; j++) {
                Log.d(String.format("%4d ", transitionCount[i][j]));
            }
            Log.dLine();

        }
        for (int i = 0; i < endTransitionCount.length; i++) {
            Log.dLine(END_TRANSITIONS[i] + ":  " + String.format("%4d ", endTransitionCount[i]));
        }


        // calc Transition Prob ----------------------------------------------------------------------
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

        transitionProb = new double[TRANSITIONS.length][lengthModel + 1];
        for (int i = 0; i < lengthModel + 1; i++) {

            int[] sumTrans = new int[STATES.length]; // vergleichsgroessen
            for (int j = 0; j < transitionCount.length; j++) {
                Trans trans = TRANSITIONS[j];
                sumTrans[stateToIndex(trans.getStartState())] += transitionCount[j][i] + PSEUDO_COUNT_TRANSITION;
            }

            for (int j = 0; j < transitionCount.length; j++) {
                Trans trans = TRANSITIONS[j];
                transitionProb[j][i] = ((double) transitionCount[j][i] + PSEUDO_COUNT_TRANSITION) / sumTrans[stateToIndex(trans.getStartState())];
            }
        }


        // output Transition Prob
        Log.iLine();
        Log.iLine("Transition Prob: (Pseudo-Count = " + PSEUDO_COUNT_TRANSITION + ")");

        for (int i = 0; i < initTransitionProb.length; i++) {
            Log.iLine(String.format("%s:  %.2f", INIT_TRANSITIONS[i].toString(), initTransitionProb[i]));
        }
        for (int i = 0; i < transitionProb.length; i++) {
            double[] transProbVector = transitionProb[i];
            Log.i(String.format("%s:  ", TRANSITIONS[i].toString()));
            for (double prob : transProbVector) {
                Log.i(String.format("%.2f ", prob));
            }
            Log.iLine();
        }
    }


    /**
     * Implementation des Viterbi-Algorithmus für den logarithmischen Raum
     *
     * @param observations Beobachtungsfolge
     * @return Zustands-Pfad
     */
    public char[] viterbi(final char[] observations) {

        // init
        int[] observationIndices = basesToIndices(observations);
        int length = observationIndices.length;

        double[][] viterbiVar = new double[STATE_COUNT][length];
        int[][] viterbiArg = new int[STATE_COUNT][length];

        // init
        for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

            viterbiVar[stateIndex][0] = initTransitionProb[stateIndex] + emissionProbMatch[stateIndex][observationIndices[0]]; // log-space  // FIXME distinguish between match and insert
            viterbiArg[stateIndex][0] = -1;
        }

        // iterate observations indices
        for (int i = 1; i < length; i++) {
            // iterate states indices
            for (int j = 0; j < STATE_COUNT; j++) {


                //find max
                double maxProb = Double.NEGATIVE_INFINITY;
                int maxArg = -1; // maximizing argument

                for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                    double prob = viterbiVar[stateIndex][i - 1] + transitionProb[stateIndex][j] + emissionProbMatch[j][observationIndices[i]]; // log-space  // FIXME distinguish between match and insert
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
            for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

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
        x[length - 1] = STATES[zLast];

        // backtrace iterate
        for (int i = length - 1; i > 0; i--) {
            int m = viterbiArg[z[i]][i];
            z[i - 1] = m;
            x[i - 1] = STATES[m];
        }

        return x;
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

    private static int initTransitionToIndex(char start, char end) {
        return transitionToIndex(INIT_TRANSITIONS, start, end);
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
