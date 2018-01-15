package main.hmm.casino;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;

import java.io.*;
import java.util.Arrays;

/**
 * Ausfuehrbare Klasse, den Dateipfad als Parameter (-file <Path>) uebergeben bekommen muss.
 * Gibt Daten der Datei aus.
 * <p>
 * Generiert anhand der eingelesenen beobachteten Sequenz einen Zustands-Pfad mittels des Viterbi-Algorithmus in {@link CasinoHMM} und gibt diesen aus. F = Fair, L = Loaded.
 */
public class CasinoHMMMain {

    /**
     * ausfuehrbare Methode
     *
     * @param args Argumente
     */
    public static void main(String[] args) {

        // set up Parameter
        ParameterSet parameterSet = new ParameterSet();
        Setting filePath = new Setting("file", true);
        parameterSet.addSetting(filePath);

        try {
            ArgumentParser parser = new ArgumentParser(parameterSet);
            parser.parseArgs(args);
        } catch (ArgumentParserException e) { // if parameter is missing or not intended
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


        // print file data
        System.out.println(inSeqRolls);
        System.out.println(inSeqDice);
        System.out.println(inSeqViterbi);


        CasinoHMM model = new CasinoHMM();

        // generate state-path with viterbi
        char[] observations = inSeqRolls.toCharArray();

        System.out.println("Generated Viterbi Path: ");

        char[] statePath = model.viterbi(observations);
        System.out.println(String.valueOf(statePath));

        boolean equal = Arrays.equals(statePath, inSeqViterbi.toCharArray());
        String out = "loaded and generated Viterbi-State-Paths are " + (equal ? "" : "NOT ") + "equal";
        if (equal)
            System.out.println(out);
        else
            System.err.println(out);
    }
}
