package main.hmm;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;

import java.io.*;

public class HMM {

    private static final char[] OBSERVATION_SPACE = {'1', '2', '3', '4', '5', '6'};

    private enum STATE {
        FAIR, UNFAIR
    }

    private static final char[] STATE_CHAR = {'F', 'U'};

    private static final double[] INIT_PROBABILITIES = new double[]{.5d, .5d};
    private static final int STATE_COUNT = INIT_PROBABILITIES.length;


    private static final double[][] TRANSITION_MATRIX = new double[][]{{.95d, .05d}, {.1d, .9d}};
    private static final double[][] EMISSION_MATRIX = new double[][]{{1d / 6, 1d / 6, 1d / 6, 1d / 6, 1d / 6, 1d / 6}, {.1d, .1d, .1d, .1d, .1d, .5d}};

    public static void main(String[] args) {

        ParameterSet parameterSet = new ParameterSet();
        Setting filePath = new Setting("file", true);
        parameterSet.addSetting(filePath);

        try {
            ArgumentParser parser = new ArgumentParser(parameterSet);
            parser.parseArgs(args);
        } catch (ArgumentParserException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }

        // reading file
        File file = new File(filePath.getValue());
        BufferedReader bufferedReader = null;
        try {
            bufferedReader = new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException e) {
            System.err.println("ERROR: file " + file + " not found");
            System.exit(1);
        }


        System.out.println("reading " + file);

        String inSeqRolls = null, inSeqDice = null, inSeqViterbi = null;
        try {
            inSeqRolls = bufferedReader.readLine(); // first line is dice sequence
            inSeqDice = bufferedReader.readLine();
            inSeqViterbi = bufferedReader.readLine();
        } catch (IOException e) {
            System.err.println("ERROR: while reading file " + file);
            System.exit(1);
        } finally {
            try {
                bufferedReader.close();
            } catch (IOException e) {
                System.err.println("ERROR: while closing reader");
                System.exit(1);
            }
        }
        System.out.println("successfully finished reading file");

        System.out.println(inSeqRolls);
        System.out.println(inSeqDice);
        System.out.println(inSeqViterbi);

        char[] observations = inSeqRolls.toCharArray();

        char[] statePath = viterbi(observations);
        System.out.println(String.valueOf(statePath));
    }


    private static char[] viterbi(final char[] observations) {

        int[] observationIndices = observationsToIndices(observations);
        int length = observationIndices.length;

        double[][] tOne = new double[STATE_COUNT][length];
        int[][] tTwo = new int[STATE_COUNT][length];

        for (STATE state : STATE.values()) {
            int stateIndex = state.ordinal();

            tOne[stateIndex][1] = INIT_PROBABILITIES[stateIndex] * EMISSION_MATRIX[stateIndex][observationIndices[1]];
            tTwo[stateIndex][1] = -1;
        }

        for (int i = 1; i < length; i++) {

            for (STATE state : STATE.values()) {
                int stateIndex = state.ordinal();

                Result result = maxState(tOne, observationIndices, i, stateIndex);

                tOne[stateIndex][i] = result.getProbability();
                tTwo[stateIndex][i] = result.getArgument();
            }
        }

        // ENDE
        int zT = -1;
        for (STATE state : STATE.values()) {
            int stateIndex = state.ordinal();

            double prob = tOne[stateIndex][length - 1];
            if (prob > zT) {
                zT = stateIndex;
            }
        }

        int[] z = new int[length];
        char[] x = new char[length];

        z[length - 1] = zT;
        x[length - 1] = STATE_CHAR[zT];

        for (int i = length - 1; i > 0; i--) {
            int m = tTwo[z[i]][i];
            z[i - 1] = m;
            x[i - 1] = STATE_CHAR[m];
        }


        return x;
    }

    private static Result maxState(final double[][] tOne, int[] observationIndices, int i, int j) {
        double maxProbability = -1d;
        int maxArg = -1;

        for (STATE state : STATE.values()) {
            int stateIndex = state.ordinal();

            double prob = tOne[stateIndex][i - 1] * TRANSITION_MATRIX[stateIndex][j] * EMISSION_MATRIX[j][observationIndices[i]];
            if (prob > maxProbability) {
                maxProbability = prob;
                maxArg = stateIndex;
            }
        }
        return new Result(maxProbability, maxArg);
    }

    private static int[] observationsToIndices(char[] observcations) {
        int length = observcations.length;
        int[] ret = new int[length];

        for (int i = 0; i < length; i++) {
            ret[i] = obesrvationToIndex(observcations[i]);
        }

        return ret;
    }

    private static int obesrvationToIndex(char observation) {
        for (int i = 0, observation_spaceLength = OBSERVATION_SPACE.length; i < observation_spaceLength; i++) {
            char c = OBSERVATION_SPACE[i];
            if (c == observation)
                return i;
        }
        return -1;
    }

    /*
    private static double prob() {
        int fairState = charToFairState(inSeqRolls.charAt(lenght - 1));
        int unfairState = charToUnfairState(inSeqRolls.charAt(lenght - 1));
        return Math.max(viterbi(fairState, lenght) * a(fairState, 0), viterbi(unfairState, lenght) * a(unfairState, 0));
    }

    private static double viterbi(int state, int pos) {
        System.out.println("viterbi: state=" + state + " pos=" + pos);
        if (pos == 0) {
            if (state == 0)
                return 1d;
            else
                return 0d;
        }
        int fairState = charToFairState(inSeqRolls.charAt(pos - 1));
        int unfairState = charToUnfairState(inSeqRolls.charAt(pos - 1));

        return e(state, inSeqRolls.charAt(pos - 1)) * Math.max(viterbi(fairState, pos - 2) * a(fairState, state), viterbi(unfairState, pos - 2) * a(unfairState, state));
    }

    private static int charToFairState(char c) {
        return c - '1' + 1;
    }

    private static int charToUnfairState(char c) {
        return charToFairState(c) + 6;
    }

    /**
     * zustandswechsel
     *
     * @return
     */
    /*
    private static double a(int stateStart, int stateEnd) {
        if (stateStart == 0 || stateEnd == 0)
            return 1d;
        if (stateStart < 7) {
            if (stateEnd < 7)
                return 0.95d; // fair to fair
            else
                return .05d; // fair to unfair
        } else {
            if (stateEnd < 7)
                return 0.1d; // unfair to fair
            else
                return 0.9d; // unfair to unfair
        }
    }

    private static double e(int state, char emission) {
        return INIT_PROBABILITIES[state];
    }
    */
}
