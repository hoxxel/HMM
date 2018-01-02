package main.hmm;

import main.argparser.ArgumentParser;
import main.argparser.ArgumentParserException;
import main.argparser.ParameterSet;
import main.argparser.Setting;

import java.io.*;
import java.util.Arrays;

/**
 * Aufgabe: Implementieren Sie den Viterbi-Algorithmus fuer das HMM zu dem Beispiel des
 * unehrlichen Casinos, wie in der Vorlesung vorgestellt (siehe auch im Buch von
 * Durbin et al. "Biological sequence analysis" Seite 54-57).
 * <p>
 * Anhand  https://en.wikipedia.org/wiki/Viterbi_algorithm  implementiert.
 * <p>
 * Ausfuehrbare Klasse, den Dateipfad als Parameter (-file <Path>) uebergeben bekommen muss.
 * Gibt Daten der Datei aus.
 * Generiert anhand der eingelesenen beobachteten Sequenz einen Zustands-Pfad mittels des Viterbi-Algorithmus und gibt diesen aus. F = Fair, L = Loaded.
 */
public class HMM {

    /**
     * beobachtbare Ereignisse
     */
    private static final char[] OBSERVATION_SPACE = {'1', '2', '3', '4', '5', '6'};

    /**
     * Zustaende in denen sich das Modell befindet
     */
    private static final char[] STATE_CHAR = {'F', 'L'};

    /**
     * Anzahl der Zustaende
     */
    private static final int STATE_COUNT = STATE_CHAR.length;

    /**
     * Uebergangswahrscheinlichen aus dem Startzustand in die jeweiligen Zustaende
     */
    private static final double[] INIT_PROBABILITIES = new double[]{.5d, .5d};

    /**
     * Uebergangswahrscheinlichen zwischen den Zustaenden
     */
    private static final double[][] TRANSITION_MATRIX = new double[][]{{.95d, .05d}, {.1d, .9d}};
    /**
     * Beobachtungswahrscheinlichketen der Ereignisse in den jeweiligen Zustaenden
     */
    private static final double[][] EMISSION_MATRIX = new double[][]{{1d / 6, 1d / 6, 1d / 6, 1d / 6, 1d / 6, 1d / 6}, {.1d, .1d, .1d, .1d, .1d, .5d}};

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

        // generate state-path with viterbi
        char[] observations = inSeqRolls.toCharArray();

        System.out.println("generated Viterbi Path: ");

        char[] statePath = viterbi(observations);
        System.out.println(String.valueOf(statePath));

        System.out.println("loaded and generated Viterbi-State-Paths are " + (Arrays.equals(statePath, inSeqViterbi.toCharArray()) ? "" : "NOT ") + "equal");
    }

    /**
     * Implementation des Viterbi-Algorithmus
     *
     * @param observations Beobachtungsfolge
     * @return Zustands-Pfad
     */
    private static char[] viterbi(final char[] observations) {

        // init
        int[] observationIndices = observationsToIndices(observations);
        int length = observationIndices.length;

        double[][] tOne = new double[STATE_COUNT][length];
        int[][] tTwo = new int[STATE_COUNT][length];

        for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

            tOne[stateIndex][0] = INIT_PROBABILITIES[stateIndex] * EMISSION_MATRIX[stateIndex][observationIndices[0]];
            tTwo[stateIndex][0] = -1;
        }

        // iterate
        for (int i = 1; i < length; i++) {

            for (int j = 0; j < STATE_COUNT; j++) {


                //find max
                double maxProbability = -1d;
                int maxArg = -1;

                for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                    double prob = tOne[stateIndex][i - 1] * TRANSITION_MATRIX[stateIndex][j] * EMISSION_MATRIX[j][observationIndices[i]];
                    if (prob > maxProbability) {
                        maxProbability = prob;
                        maxArg = stateIndex;
                    }
                }

                tOne[j][i] = maxProbability;
                tTwo[j][i] = maxArg;
            }
        }

        // backtrace init
        int zLast = -1;
        {
            double probLast = 0d;
            for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                double prob = tOne[stateIndex][length - 1];
                if (prob > probLast) {
                    zLast = stateIndex;
                    probLast = prob;
                }
            }
        }

        int[] z = new int[length]; // stateIndexPath
        char[] x = new char[length]; // statePath

        z[length - 1] = zLast;
        x[length - 1] = STATE_CHAR[zLast];

        // backtrace iterate
        for (int i = length - 1; i > 0; i--) {
            int m = tTwo[z[i]][i];
            z[i - 1] = m;
            x[i - 1] = STATE_CHAR[m];
        }

        return x;
    }

    /**
     * mappt Beaobachtung-Folge auf entsprechende Index-Folge
     *
     * @param observations Beobachtungs-Folge
     * @return entsprechende Index-Folge
     */
    private static int[] observationsToIndices(final char[] observations) {
        int length = observations.length;
        int[] ret = new int[length];

        for (int i = 0; i < length; i++) {
            ret[i] = obesrvationToIndex(observations[i]);
        }

        return ret;
    }

    /**
     * mappt Beaobachtung auf den entsprechenden Index
     *
     * @param observation Beobachtung
     * @return entsprechender Index
     */
    private static int obesrvationToIndex(final char observation) {
        for (int i = 0, observation_spaceLength = OBSERVATION_SPACE.length; i < observation_spaceLength; i++) {
            char c = OBSERVATION_SPACE[i];
            if (c == observation)
                return i;
        }
        return -1;
    }
}
