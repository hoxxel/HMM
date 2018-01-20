package main.hmm.profil;

import main.argparser.*;
import main.fastaparser.FastaParser;
import main.fastaparser.FastaParserException;
import main.fastaparser.Sequence;
import main.logger.Log;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

public class ProfilHMMMain {


    public static void main(String[] args) {
        // set up Parameter
        ParameterSet parameterSet = new ParameterSet();
        Setting paramFileTrain = new Setting("filetrain", true);
        Setting paramFileTest = new Setting("filetest", true);
        Flag paramDebug = new Flag("debug", false);
        parameterSet.addSetting(paramFileTrain);
        parameterSet.addSetting(paramFileTest);
        parameterSet.addFlag(paramDebug);

        try {
            ArgumentParser parser = new ArgumentParser(parameterSet);
            parser.parseArgs(args);
        } catch (ArgumentParserException e) { // if parameter is missing or not intended
            Log.eLine(e.getMessage());
            System.exit(1);
        }

        if (paramDebug.isSet())
            Log.setPrintDebug(true);

        List<Sequence> sequencesTrain = readFile(paramFileTrain.getValue());

        ProfilHMM model = null;
        try {
            model = new ProfilHMM(sequencesTrain);
        } catch (IllegalArgumentException e) {
            Log.eLine("ERROR: Building ProfilHMM failed! " + e.getMessage());
            System.exit(1);
        }

        Log.iLine();

        List<Sequence> sequencesTest = readFile(paramFileTest.getValue());

        for (int i = 0; i < sequencesTest.size(); i++) {

            char[] observ = sequencesTest.get(i).getSequence().toCharArray();

            char[] statePath = null;
            try {
                statePath = model.viterbi(observ);
            } catch (IllegalArgumentException e) {
                Log.eLine("ERROR: Viterbi ProfilHMM failed! " + e.getMessage());
                System.exit(1);
            }

            Log.iLine(String.valueOf(observ));
            Log.iLine(String.valueOf(statePath));
            Log.iLine();
        }
    }

    private static List<Sequence> readFile(String filePath) {

        List<Sequence> ret = null;

        Log.iLine("reading " + filePath);
        try {
            ret = FastaParser.parseFile(filePath);
        } catch (FileNotFoundException e) {
            Log.eLine("ERROR: file " + filePath + " not found");
            System.exit(1);
        } catch (IOException e) {
            Log.eLine("ERROR: while reading file " + filePath);
            System.exit(1);
        } catch (FastaParserException e) {
            Log.eLine("ERROR: while parsing file " + filePath + ": " + e.getMessage());
            System.exit(1);
        }

        Log.iLine("successfully finished reading file");

        return ret;
    }
}
