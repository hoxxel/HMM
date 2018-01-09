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

    private static final char[] BASES = {'-', 'A', 'C', 'G', 'U'};

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
        int[][] freq = new int[BASES.length][length];

        for (String seq : sequencesTrain) {
            for (int i = 0; i < length; i++) {
                char base = seq.charAt(i);
                freq[baseToIndex(base)][i]++;
            }
        }

        for (int[] a : freq) {
            System.out.println(Arrays.toString(a));
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
