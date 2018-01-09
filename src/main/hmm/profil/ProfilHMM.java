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

    public static final int PSEUDO_COUNT = 1;

    private static final char GAP = '-';
    private static final char[] BASES = {'A', 'C', 'G', 'U'};
    private static final char[] STATES = {'M', 'I', 'D'};

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
            int seqCountHalf = sequencesTrain.size() / 2;
            for (int i = 0; i < gapCounts.length; i++) {
                int gapCount = gapCounts[i];
                if (gapCount <= seqCountHalf) {
                    lengthModel++;
                    matchState[i] = true;
                }
            }
        }

        System.out.println(length);
        System.out.println("Model length = " + lengthModel);

        System.out.println(GAP + " " + Arrays.toString(gapCounts));

/*
        for (int i = 0; i < baseCounts.length; i++) {
            int[] a = baseCounts[i];
            System.out.println(BASES[i] + " " + Arrays.toString(a));
        }
        */

        double[][][] transitionProb = new double[lengthModel][3][3];
        double[][] emissionProb = new double[lengthModel][BASES.length];

        for (int i = 0, indexModel = 0; i < length; i++) {
            if (matchState[i]) { // if at this column is match-state
                double[] emissionProbVector = emissionProb[indexModel];
                int[] baseCountVector = baseCounts[i];

                int sumEmissionCount = 0;

                for (int baseCount : baseCountVector) {
                    sumEmissionCount += baseCount + PSEUDO_COUNT;
                }

                for (int j = 0; j < emissionProbVector.length; j++) {
                    emissionProbVector[j] = ((double) (PSEUDO_COUNT + baseCountVector[j])) / sumEmissionCount;
                }

                indexModel++;
            }
        }

        int sumTransitionFreq;


        // output
        for (int i = 0, j = 0; i < baseCounts.length; i++) {
            System.out.print(Arrays.toString(baseCounts[i]));
            if (matchState[i]) {
                System.out.print("  " + Arrays.toString(emissionProb[j]));
                j++;
            }
            System.out.println();
        }

    }

    private static int baseToIndex(char base) {
        for (int i = 0, basesLength = BASES.length; i < basesLength; i++) {
            char c = BASES[i];
            if (c == base)
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
