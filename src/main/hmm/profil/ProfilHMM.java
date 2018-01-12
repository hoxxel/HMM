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
    private static final String[] TRANSITIONS = {"MM", "MI", "MD", "IM", "II", "ID", "DM", "DI", "DD"};
    private static final String[] INIT_TRANSITIONS = {"AM", "AI", "AD"};

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

        int[][] transitionCount = new int[TRANSITIONS.length][lengthModel + 1];
        int[] initTransitionCount = new int[INIT_TRANSITIONS.length];
        double[][] transitionProb = new double[TRANSITIONS.length][lengthModel + 1];


        for (int i = 0; i < length; i++) {
            System.out.print(matchState[i] ? "m  " : "\u001B[32mI  \u001B[0m");
        }
        System.out.println();
        for (int i1 = 0; i1 < 7; i1++) {
            String seq = sequencesTrain.get(i1);

            // output sequence
            System.out.println("\u001B[37m");
            for (int i = 0; i < length; i++) {
                System.out.print(seq.charAt(i) + "  ");
            }
            System.out.println("\u001B[0m");


            // get States and Transitions
            char lastState = 'A';
            int insertCount = 0;

            for (int i = 0, iModel = 1; i < length; i++) {
                char state = getState(seq, matchState, i);
                if (state != ' ') {


                    if (lastState == 'A') {
                        initTransitionCount[initTransitionToIndex(lastState + "" + state)]++;
                    } else {

                        if (state == 'I') {
                            insertCount++;
                        }
                        // End of InsertStates
                        else {
                            if (lastState == 'I') {
                                transitionCount[transitionToIndex("II")][iModel - 1] += insertCount;
                                insertCount = 0;
                            }
                            String transition = lastState + "" + state;
                            int matrixIndex = iModel;
                            if (transition.equals("IM") || transition.equals("ID"))
                                matrixIndex--;
                            transitionCount[transitionToIndex(transition)][matrixIndex]++;


                            iModel++;
                        }
                    }

                    lastState = state;
                }
                System.out.print(state + "  ");
            }
            System.out.println();
        }


        System.out.println();
        for (int i = 0; i < initTransitionCount.length; i++) {
            System.out.println(INIT_TRANSITIONS[i] + " " + initTransitionCount[i]);
        }
        for (int i = 0; i < transitionCount.length; i++) {
            System.out.println(TRANSITIONS[i] + " " + Arrays.toString(transitionCount[i]));
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

    }

    private static char getNextState(String seq, boolean[] matchState, int index) {
        index++;

        if (index >= seq.length())
            return 'E';

        return getState(seq, matchState, index);
    }

    private static char getState(String seq, boolean[] matchState, int index) {

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

    private static int baseToIndex(char base) {
        return charToIndex(BASES, base);
    }

    private static int stateToIndex(char state) {
        return charToIndex(STATES, state);
    }

    private static int transitionToIndex(String transition) {
        for (int i = 0, transitionsLength = TRANSITIONS.length; i < transitionsLength; i++) {
            if (TRANSITIONS[i].equals(transition))
                return i;
        }
        throw new IllegalArgumentException("Transition " + transition + " not found");
    }

    private static int initTransitionToIndex(String transition) {
        for (int i = 0, transitionsLength = INIT_TRANSITIONS.length; i < transitionsLength; i++) {
            if (INIT_TRANSITIONS[i].equals(transition))
                return i;
        }
        throw new IllegalArgumentException("Transition " + transition + " not found");
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
