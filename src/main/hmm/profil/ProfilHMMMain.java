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
     * ausfuehrbare Methode
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
        List<ViterbiPath> viterbiPaths = Collections.synchronizedList(new ArrayList<>(sequenceCount));
        int coreCount = Runtime.getRuntime().availableProcessors();
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

        // TODO implement rRNA detection

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
