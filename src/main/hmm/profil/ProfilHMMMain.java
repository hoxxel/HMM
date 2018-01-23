package main.hmm.profil;

import main.argparser.*;
import main.fastaparser.FastaParser;
import main.fastaparser.FastaParserException;
import main.fastaparser.Sequence;
import main.logger.Log;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Ausfuehrbare Klasse, die den Dateipfad der Traings-Sequencen als Parameter (-filetrain <Path>)
 * sowie der Test-Sequencen als Parameter (-filetest <Path>) uebergeben bekommen muss.
 * <p>
 * Erstellt anhand der Trainings-Sequencen ein {@link ProfilHMM}.
 * Anschliessend wird mittels des Viterbi-Algorithmus fuer jede Test-Sequenz ein Zustands-Pfad ermittelt.
 *
 * @author Soeren Metje
 */
public class ProfilHMMMain {

    /**
     * ausfuehrbare Methode. Erstellt anhand der Traings-Sequnzen ein {@link ProfilHMM} und
     * berechnet mittels des Viterbi-Algorithmus den wahrscheinlichsten Zustands-Pfad.
     *
     * @param args Argumente
     */
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

        // Test-Sequences --------------------------------------------------------
        List<Sequence> sequencesTest = readFile(paramFileTest.getValue());

        int sequenceCount = sequencesTest.size();
        int coreCount = Runtime.getRuntime().availableProcessors(); // returns count of logical cores available to JVM
        Log.dLine("available Cores = " + coreCount);
        // Init queues
        List<Queue<Sequence>> queues = new ArrayList<>(coreCount);
        for (int i = 0; i < coreCount; i++) {
            queues.add(new LinkedList<>());
        }
        // Distribute sequences to queues
        {
            int i = 0;
            for (Sequence sequence : sequencesTest) {
                queues.get(i % coreCount).add(sequence);
                i++;
            }
        }
        // Create and Start Threads
        int threadCount = 0;
        List<ViterbiPath> viterbiPaths = Collections.synchronizedList(new ArrayList<>(sequenceCount)); // list to hold results created in Threads
        Queue<Thread> threads = new LinkedList<>();
        for (Queue<Sequence> queue : queues) {
            if (!queue.isEmpty()) {
                Thread thread = new ThreadViterbi(model, queue, viterbiPaths);
                threads.add(thread);
                threadCount++;
                thread.start();
            }
        }
        Log.iLine(threadCount + " Threads created and started running Viterbi-Algo for " + sequenceCount + " Test-Sequences");
        Log.iLine("Waiting for async Output...");

        // Waiting for threads to finish
        while (!threads.isEmpty()) {
            Thread thread = threads.poll();
            try {
                thread.join();
                Log.dLine(thread.getName() + " finished");
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        // all Threads finished

        double threshold = calcThreshold(viterbiPaths);
        Log.iLine(String.format("Threshold = %.2f", threshold));

        Log.iLine("Table log Score and rRNA-Decision:");
        {
            StringBuilder out = new StringBuilder();
            for (ViterbiPath path : viterbiPaths) {
                double score = path.getScore();
                out.append(String.format("%.3f;%c\n", score, (score >= threshold ? '1' : '0')));
            }
            Log.iLine(out.toString());
        }
    }

    /**
     * berechnet den Score-Schwellwert der Sequenzen {@link Sequence} bzw. Zusatnds-Pfade {@link ViterbiPath}, mit dem zwischen rRNA und NonrRNA unterschieden werden soll.
     * Liefert diesen abschliessend zurueck.
     *
     * @param viterbiPaths zu betrachtende Sequenzen {@link Sequence} bzw. Zusatnds-Pfade {@link ViterbiPath}
     * @return Score-Schwellwert
     */
    private static double calcThreshold(final List<ViterbiPath> viterbiPaths) { // TODO improve
        int stateCount = 2;
        final ArrayList<List<ViterbiPath>> listrRNA = new ArrayList<>(stateCount);
        listrRNA.add(new LinkedList<>()); // 0 = rRNA
        listrRNA.add(new LinkedList<>()); // 1 = NonrRNA

        {
            // average score of all statepaths
            double avgScorePerState = 0d;
            for (ViterbiPath path : viterbiPaths) {
                avgScorePerState += path.getScore() / path.getStatePath().length;
            }
            avgScorePerState /= viterbiPaths.size();

            // use avarage as threshold and classify the sequences
            for (ViterbiPath path : viterbiPaths) {
                int listIndex = isrRNA(path, avgScorePerState) ? 1 : 0;
                listrRNA.get(listIndex).add(path);
            }
        }

        // calc averages of scores from classified sequences
        double[] avgScore = new double[stateCount];
        {
            int i = 0;
            for (List<ViterbiPath> list : listrRNA) {
                for (ViterbiPath path : list) {
                    avgScore[i] += path.getScore();
                }
                avgScore[i] /= list.size();
                i++;
            }
        }

        // place threshold in between average-rRNA and average-NonrRNA
        double threshold = 0;
        for (double avg : avgScore) {
            threshold += avg;
        }
        threshold /= stateCount;
        return threshold;
    }

    /**
     * Klassifiziert Sequnze {@link Sequence} bzw Zustands-Pfad {@link ViterbiPath} als rRNA oder NonrRNA
     * anhand uebergebenen Schwellwertes bzgl. durchschnittlichem Score pro Zustand des Pfades.
     * <p>
     * Liefert true zurueck, falls rRNA. Ansonsten false.
     *
     * @param path                      Sequnze {@link Sequence} bzw Zustands-Pfad {@link ViterbiPath}
     * @param thresholdAvgScorePerState Schwellwertes bzgl. durchschnittlichem Score pro Zustand
     * @return true, falls rRNA. Ansonsten false.
     */
    private static boolean isrRNA(final ViterbiPath path, double thresholdAvgScorePerState) { // TODO improve
        final char[] statePath = path.getStatePath();
        final double score = path.getScore();

        double scoreAveragePerState = score / statePath.length;

        if (scoreAveragePerState >= thresholdAvgScorePerState)
            return false;
        return true;
    }

    /**
     * Liesst Sequenzen aus Datei an uebergebenem Pfad mittels {@link FastaParser} ein und liefert sie zurueck.
     *
     * @param filePath Pfad zu Datei
     * @return Liste mit Sequenzen
     */
    private static List<Sequence> readFile(final String filePath) {

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
