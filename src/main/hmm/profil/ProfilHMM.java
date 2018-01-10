package main.hmm.profil;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class ProfilHMM {

    public static final int PSEUDO_COUNT_EMISSION = 1;
    public static final int PSEUDO_COUNT_TRANSITION = 1;
    public static final double THRESHOLD_MATCHSTATE = .5d; // min amount of nucleotides (no gap) to count column as match-state

    private static final char GAP = '-';
    private static final char[] BASES = {'A', 'C', 'G', 'U'};
    private static final char[] STATES = {'M', 'I', 'D'};
    private static final String[] TRANSITIONS = {"MM", "MI", "MD", "IM", "II", "DM", "DD"};

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

        System.out.println(length);
        System.out.println("Model length = " + lengthModel);

        /*
        System.out.println(GAP + " " + Arrays.toString(gapCounts));


        for (int i = 0; i < baseCounts.length; i++) {
            int[] a = baseCounts[i];
            System.out.println(BASES[i] + " " + Arrays.toString(a));
        }
        */

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

        int[][] transitionCount = new int[length - 1][TRANSITIONS.length];
        double[][] transitionProb = new double[length - 1][TRANSITIONS.length];


        for (int i = 0; i < length; i++) {
            System.out.print(matchState[i] ? "m  " : "i  ");
        }
        for (int i1 = 0; i1 < 10; i1++) {
            System.out.println();
            String seq = sequencesTrain.get(i1);
            char lastState = getState(seq, matchState, 0);

            for (int i = 0; i < length; i++) {
                System.out.print(seq.charAt(i) + "  ");
            }
            System.out.println();
            System.out.print(" ");
            for (int i = 1, indexModel = 0; i < length; i++) {

                String transition = "  ";
                char state = getState(seq, matchState, i);
                if (matchState[i]) {
                    transition = lastState + "" + state;
                    lastState = state;
                    indexModel++;
                } else {
                    if (state == 'I' && lastState == 'I' && seq.charAt(i - 1) == GAP) {
                        transition = "  ";
                        lastState = state;
                    } else if (state == 'I') {
                        transition = lastState + "" + state;
                        lastState = state;
                    }
                }

                int pos = i - 1;
                int indexTrans = transitionToIndex(transition);
                if (indexTrans > 0) {
                    transitionCount[pos][indexTrans]++;
                }


                System.out.print(transition + " ");

            }
        }

        System.out.println();
        System.out.println(Arrays.toString(TRANSITIONS));
        for (int i = 0; i < 20; i++) {
            System.out.println(Arrays.toString(transitionCount[i]));
        }


        /*

        // output
        for (int i = 0; i < baseCounts.length; i++) {
            System.out.print(gapCounts[i]);
            System.out.print("  " + Arrays.toString(baseCounts[i]));
            System.out.print("  " + Arrays.toString(emissionProb[i]));

            System.out.println();
        }
        */

        // test
        String[] strings = new String[]{"UACAAU", "UA-AAU", "U--AAU", "U-CAAU", "U--AAU", "U--AAU", "UACAAU", "U--AAU", "U--AAU"};
        for (String s : strings) {
            System.out.println(getState(s, new boolean[]{true, false, false, true, true, true}, 0));
        }

    }

    private static char getNextState(String seq, boolean[] matchState, int index) {
        index++;

        if (index >= seq.length())
            return 'E';

        return getState(seq, matchState, index);
    }

    private static char getState(String seq, boolean[] matchState, int index) {


        int length = seq.length();

        // iterate Insert-states
        while (!matchState[index]) {
            if (seq.charAt(index) != GAP) {
                return 'I';
            }
            index++;
            if (index >= length)
                return 'E';
        }

        // at i is a Match-state
        if (seq.charAt(index) == GAP)
            return 'D';
        return 'M';
    }

    private static int baseToIndex(char base) {
        return charToIndex(BASES, base);
    }

    private static int stateToIndex(char base) {
        return charToIndex(STATES, base);
    }

    private static int transitionToIndex(String s) {
        for (int i = 0, basesLength = TRANSITIONS.length; i < basesLength; i++) {
            if (TRANSITIONS[i].equals(s))
                return i;
        }
        return -1;
    }

    private static int charToIndex(char[] array, char c) {
        for (int i = 0, basesLength = array.length; i < basesLength; i++) {
            if (array[i] == c)
                return i;
        }
        return -1;
    }

    private static List<String> readFile(String filePath) {

        List<String> ret = null;

        System.out.println("reading " + filePath);
        try {
            ret = FastaParser.parseFile(filePath);
        } catch (FileNotFoundException e) {
            System.err.println("ERROR: file " + filePath + " not found");
            System.exit(1);
        } catch (IOException e) {
            System.err.println("ERROR: while reading file " + filePath);
            System.exit(1);
        }

        System.out.println("successfully finished reading file");

        return ret;
    }

}
