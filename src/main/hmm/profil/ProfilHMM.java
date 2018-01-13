package main.hmm.profil;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;
import main.logger.Log;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
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
    private static final char[] STATES = {'M', 'I', 'D'};
    private static final String[] TRANSITIONS = {"MM", "MI", "MD", "IM", "II", "ID", "DM", "DI", "DD"};
    private static final String[] INIT_TRANSITIONS = {"BM", "BI", "BD"};
    private static final String[] END_TRANSITIONS = {"ME", "IE", "DE"};

    public static void main(String[] args) {
        // set up Parameter
        ParameterSet parameterSet = new ParameterSet();
        Setting filePathTrain = new Setting("filetrain", true);
        Setting filePathTestF = new Setting("filetestshort", false); //TODO req
        Setting filePathTestS = new Setting("filetrainfull", false); // TODO req
        parameterSet.addSetting(filePathTrain);
        parameterSet.addSetting(filePathTestF);
        parameterSet.addSetting(filePathTestS);

        try {
            ArgumentParser parser = new ArgumentParser(parameterSet);
            parser.parseArgs(args);
        } catch (ArgumentParserException e) { // if parameter is missing or not intended
            System.err.println(e.getMessage());
            System.exit(1);
        }

        List<String> sequencesTrain = readFile(filePathTrain.getValue());
        List<String> sequencesTestF;
        List<String> sequencesTestS;

        if (filePathTestF.isSet()) {
            sequencesTestF = readFile(filePathTestF.getValue());
        }
        if (filePathTestS.isSet()) {
            sequencesTestS = readFile(filePathTestS.getValue());
        }

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

        Log.dLine(length);
        Log.dLine("Model length = " + lengthModel);


        // calc Emission Prob ---------------------------------------------------------------------------

        double[][] emissionProb = new double[length][BASES.length]; // emission-probability for match- and insert states

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

        Log.dLine("Gap-Counts:");
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

/*
        for (int i = 0; i < baseCounts.length; i++) {
            Log.d(gapCounts[i]);
            Log.d("  " + Arrays.toString(baseCounts[i]));
            Log.d("  " + Arrays.toString(emissionProb[i]));

            Log.dLine();
        }
        */

        // calc Transition Prob ----------------------------------------------------------------------

        int[][] transitionCount = new int[TRANSITIONS.length][lengthModel + 1];
        int[] initTransitionCount = new int[INIT_TRANSITIONS.length];
        int[] endTransitionCount = new int[END_TRANSITIONS.length];

        double[][] transitionProb = new double[TRANSITIONS.length][lengthModel + 1];


        for (int i = 0; i < length; i++) {
            Log.d(matchState[i] ? "m  " : "\u001B[32mI  \u001B[0m");
        }
        Log.dLine();
        for (int i1 = 0; i1 < 10; i1++) {
            String seq = sequencesTrain.get(i1);

            // output sequence
            Log.dLine("\u001B[37m");
            for (int i = 0; i < length; i++) {
                Log.d(seq.charAt(i) + "  ");
            }
            Log.dLine("\u001B[0m");


            // get States and Transitions
            char lastState = 'B';
            int insertCount = 0;

            for (int i = 0, iModel = 0; i < length + 1; i++) {
                char state = getState(seq, matchState, i);
                if (state != ' ') {
                    if (state == 'E') {
                        endTransitionCount[endTransitionToIndex(lastState + "" + state)]++;
                    } else if (lastState == 'B') {
                        initTransitionCount[initTransitionToIndex(lastState + "" + state)]++;
                        if (state == 'M' || state == 'D') // no Insert-State in first column
                            iModel++;
                    } else {

                        if (state == 'I') {
                            insertCount++;
                        }
                        // End of InsertStates
                        else {
                            if (lastState == 'I') {
                                transitionCount[transitionToIndex("II")][iModel] += insertCount;
                                insertCount = 0;
                            }
                            String transition = lastState + "" + state;
                            transitionCount[transitionToIndex(transition)][iModel]++;

                            iModel++;
                        }
                    }

                    lastState = state;
                }
                Log.d(state + "  ");
            }
            Log.dLine();
        }


        Log.dLine();
        for (int i = 0; i < initTransitionCount.length; i++) {
            Log.dLine(INIT_TRANSITIONS[i] + " " + initTransitionCount[i]);
        }
        for (int i = 0; i < transitionCount.length; i++) {
            Log.dLine(TRANSITIONS[i] + " " + Arrays.toString(transitionCount[i]));
        }
        for (int i = 0; i < endTransitionCount.length; i++) {
            Log.dLine(END_TRANSITIONS[i] + " " + endTransitionCount[i]);
        }

    }

    private static char getState(String seq, boolean[] matchState, int index) {

        if (index >= seq.length())
            return 'E';

        boolean match = matchState[index];

        if (!match) {
            if (seq.charAt(index) != GAP) {
                return 'I';
            } else {
                return ' ';
            }
        }

        // at i is a Match-state
        if (seq.charAt(index) == GAP)
            return 'D';
        return 'M';
    }

    private static int transitionToIndex(String transition) {
        return transitionToIndex(TRANSITIONS, transition);

    }

    private static int initTransitionToIndex(String transition) {
        return transitionToIndex(INIT_TRANSITIONS, transition);
    }

    private static int endTransitionToIndex(String transition) {
        return transitionToIndex(END_TRANSITIONS, transition);
    }

    private static int transitionToIndex(String[] array, String transition) {
        for (int i = 0, arrayLength = array.length; i < arrayLength; i++) {
            if (array[i].equals(transition))
                return i;
        }
        throw new IllegalArgumentException("Transition " + transition + " not found");
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

    private static List<String> readFile(String filePath) {

        List<String> ret = null;

        Log.dLine("reading " + filePath);
        try {
            ret = FastaParser.parseFile(filePath);
        } catch (FileNotFoundException e) {
            System.err.println("ERROR: file " + filePath + " not found");
            System.exit(1);
        } catch (IOException e) {
            System.err.println("ERROR: while reading file " + filePath);
            System.exit(1);
        }

        Log.dLine("successfully finished reading file");

        return ret;
    }

}
