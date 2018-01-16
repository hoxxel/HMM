package main.hmm.profil;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;
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
        parameterSet.addSetting(paramFileTrain);
        parameterSet.addSetting(paramFileTest);

        try {
            ArgumentParser parser = new ArgumentParser(parameterSet);
            parser.parseArgs(args);
        } catch (ArgumentParserException e) { // if parameter is missing or not intended
            Log.e(e.getMessage());
            System.exit(1);
        }

        List<Sequence> sequencesTrain = readFile(paramFileTrain.getValue());

        ProfilHMM model = new ProfilHMM(sequencesTrain);

        Log.iLine();

        List<Sequence> sequencesTest = readFile(paramFileTest.getValue());

        for (int i = 0; i < sequencesTest.size(); i++) {

            char[] observ = sequencesTest.get(i).getSequence().toCharArray();
            char[] statePath = model.viterbi(observ);

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
            Log.e("ERROR: file " + filePath + " not found");
            System.exit(1);
        } catch (IOException e) {
            Log.e("ERROR: while reading file " + filePath);
            System.exit(1);
        } catch (FastaParserException e) {
            Log.e("ERROR: while parsing file " + filePath + ": " + e.getMessage());
        }

        Log.iLine("successfully finished reading file");

        return ret;
    }
}
