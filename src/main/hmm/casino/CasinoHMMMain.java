package main.hmm.casino;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;
import main.logger.Log;

import java.io.*;
import java.util.Arrays;

/**
 * <p>
 * Ausfuehrbare Klasse, die den Dateipfad als Parameter (-file <Path>) uebergeben bekommen muss.
 * Gibt Daten der Datei aus.
 * </p>
 * <p>
 * Generiert anhand der eingelesenen beobachteten Sequenz einen Zustands-Pfad mittels des Viterbi-Algorithmus in {@link CasinoHMM}
 * und gibt diesen aus. F = Fair, L = Loaded.
 * </p>
 *
 * <p>
 * Aufgabe: Implementieren Sie den Viterbi-Algorithmus fuer das HMM zu dem Beispiel des
 * unehrlichen Casinos, wie in der Vorlesung vorgestellt (siehe auch im Buch von
 * Durbin et al. "Biological sequence analysis" Seite 54-57).
 * </p>
 *
 * @author Soeren Metje
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
            Log.eLine(e.getMessage());
            System.exit(1);
        }

        // reading file
        File file = new File(filePath.getValue());
        BufferedReader bufferedReader = null;
        try {
            bufferedReader = new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException e) {
            Log.eLine("ERROR: file " + file + " not found");
            System.exit(1);
        }

        Log.iLine("reading " + file);

        String inSeqRolls = null, inSeqDice = null, inSeqViterbi = null;
        try {
            inSeqRolls = bufferedReader.readLine(); // first line is dice sequence
            inSeqDice = bufferedReader.readLine();
            inSeqViterbi = bufferedReader.readLine();
        } catch (IOException e) {
            Log.eLine("ERROR: while reading file " + file);
            System.exit(1);
        } finally {
            try {
                bufferedReader.close();
            } catch (IOException e) {
                Log.eLine("ERROR: while closing reader");
                System.exit(1);
            }
        }
        Log.iLine("successfully finished reading file");


        // print file data
        Log.iLine(inSeqRolls);
        Log.iLine(inSeqDice);
        Log.iLine(inSeqViterbi);


        CasinoHMM model = new CasinoHMM();

        // generate state-path with viterbi
        char[] observations = inSeqRolls.toCharArray();

        Log.iLine("Generated Viterbi Path: ");

        char[] statePath = null;
        try {
            statePath = model.viterbi(observations);
        } catch (IllegalArgumentException e) {
            Log.eLine("ERROR: Viterbi Casino failed! " + e.getMessage());
            System.exit(1);
        }
        Log.iLine(String.valueOf(statePath));

        boolean equal = Arrays.equals(statePath, inSeqViterbi.toCharArray());
        String out = "loaded and generated Viterbi-State-Paths are " + (equal ? "" : "NOT ") + "equal";
        if (equal)
            Log.iLine(out);
        else
            Log.eLine(out);
    }
}
