package main.hmm.profil;

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
    private static final Trans[] TRANSITIONS = {new Trans(STATE_MATCH, STATE_MATCH), new Trans(STATE_MATCH, STATE_INSERT), new Trans(STATE_MATCH, STATE_DELETE),
            new Trans(STATE_INSERT, STATE_MATCH), new Trans(STATE_INSERT, STATE_INSERT), new Trans(STATE_INSERT, STATE_DELETE),
            new Trans(STATE_DELETE, STATE_MATCH), new Trans(STATE_DELETE, STATE_INSERT), new Trans(STATE_DELETE, STATE_DELETE)};
    private static final Trans[] INIT_TRANSITIONS = {new Trans(STATE_BEGIN, STATE_MATCH), new Trans(STATE_BEGIN, STATE_INSERT), new Trans(STATE_BEGIN, STATE_DELETE)};
    private static final Trans[] END_TRANSITIONS = {new Trans(STATE_MATCH, STATE_END), new Trans(STATE_INSERT, STATE_END), new Trans(STATE_DELETE, STATE_END)};


    private double[][] emissionProb;
    private double[] initTransitionProb;
    private double[][] transitionProb;

    public ProfilHMM() {
    }

    public ProfilHMM(List<String> sequencesTrain) {
        buildModel(sequencesTrain);
    }

    public void buildModel(List<String> sequencesTrain) {

        // find Match or Insertion-States and Model length --------------------------------------------------------
        int length = sequencesTrain.get(0).length();
        int[] gapCounts = new int[length];
        int[][] baseCounts = new int[length][BASES.length];

        for (String seq : sequencesTrain) {
            for (int i = 0; i < length; i++) {
                char base = seq.charAt(i);
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

        emissionProb = new double[length][BASES.length]; // emission-probability for match- and insert states

        for (int i = 0; i < length; i++) {
            double[] emissionProbVector = emissionProb[i];
            int[] baseCountVector = baseCounts[i];

            int sumEmissionCount = 0;

            for (int baseCount : baseCountVector) {
                sumEmissionCount += baseCount + PSEUDO_COUNT_EMISSION;
            }

            for (int j = 0; j < emissionProbVector.length; j++) {
                emissionProbVector[j] = ((double) (PSEUDO_COUNT_EMISSION + baseCountVector[j])) / sumEmissionCount;
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
        Log.dLine("Nucleotide-Counts and Emission Prob:");
        for (int i = 0; i < baseCounts[0].length; i++) {
            for (int j = 0; j < baseCounts.length; j++) {

                Log.d(String.format("%4d ", baseCounts[j][i]));
            }
            Log.dLine();
            for (int j = 0; j < emissionProb.length; j++) {

                Log.d(String.format("%.2f ", emissionProb[j][i]));
            }
            Log.dLine();
        }

        Log.dLine();

        // calc Transition Count ----------------------------------------------------------------------

        int[][] transitionCount = new int[TRANSITIONS.length][lengthModel + 1];
        int[] initTransitionCount = new int[INIT_TRANSITIONS.length];
        int[] endTransitionCount = new int[END_TRANSITIONS.length];


        Log.d("    ");
        for (int i = 0; i < length; i++) {
            Log.d(matchState[i] ? "  m" : "\u001B[32m  I\u001B[0m");
        }
        Log.dLine();
        for (int i1 = 0; i1 < 10; i1++) {
            String seq = sequencesTrain.get(i1);

            // output sequence
            Log.dLine("\u001B[37m");
            Log.d(String.format("%3d:", i1));
            for (int i = 0; i < length; i++) {
                Log.d("  " + seq.charAt(i));
            }
            Log.dLine("\u001B[0m");


            // get States and Transitions
            char lastState = STATE_BEGIN;
            int insertCount = 0;

            Log.d("    ");
            for (int i = 0, iModel = 0; i < length + 1; i++) {
                char state = getState(seq, matchState, i);
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
                sumInitTrans += anInitTransitionCount;
            }

            for (int i = 0; i < initTransitionCount.length; i++) {
                initTransitionProb[i] = ((double) initTransitionCount[i]) / sumInitTrans;
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
        Log.dLine();
        Log.dLine("Transition Prob: ");
        for (int i = 0; i < transitionProb.length; i++) {
            double[] transProbVector = transitionProb[i];
            Log.d(String.format("%s:  ", TRANSITIONS[i].toString()));
            for (double prob : transProbVector) {
                Log.d(String.format("%.2f ", prob));
            }
            Log.dLine();
        }
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
