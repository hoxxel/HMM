package main.hmm.profil.viterbi.parallel;

import main.fastaparser.Sequence;
import main.hmm.profil.ProfilHMM;
import main.hmm.profil.RNAProfilHMM;
import main.hmm.profil.viterbi.ViterbiPath;
import main.logger.Log;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

/**
 * Enthaelt Methode zur parallelisierten Ausfuehrung des Viterbi-Algorithmus fuer mehrere Sequenzen.
 *
 * @author Soeren Metje
 */
public class ParallelizationSupporter {

    /**
     * Fuehrt Viterbi-Algorithmus parallelisiert aus und liefert die berechneten Zustands-Pfade {@link ViterbiPath} zurueck.
     * <p>
     * Es werden so viele Threads erstellt, wie logische Kerne der JVM zur verfuegung stehen.
     * Die Threads berechnen anhand des uebergebenen Models fuer jede Sequenz den Zustands-Pfad.
     * Abschlissend wird auf die Threads gewartet und eine Liste mit den Zustands-Pfaden {@link ViterbiPath} zurueck geliefert.
     *
     * @param model     {@link RNAProfilHMM} Modell
     * @param sequences {@link Sequence} Sequenz
     * @return Liste mit den Zustands-Pfaden {@link ViterbiPath}
     */
    public static List<ViterbiPath> viterbiParallelized(ProfilHMM model, List<Sequence> sequences) {
        int sequenceCount = sequences.size();
        int coreCount = Runtime.getRuntime().availableProcessors(); // returns count of logical cores available to JVM
        Log.dLine("available Cores = " + coreCount);
        // Init sequenceQueue
        Queue<Sequence> sequenceQueue = new LinkedList<>(sequences); // synchronisation in ThreadViterbi

        // Create and Start Threads
        int threadCount = Math.min(coreCount, sequenceQueue.size());
        Log.iLine("Creating and starting " + threadCount + " Threads running Viterbi-Algo for " + sequenceCount + " Test-Sequences");
        Log.iLine("Waiting for async Output...");
        ViterbiPath[] viterbiPaths = new ViterbiPath[sequenceCount];// list to hold results created in Threads
        Queue<Thread> threads = new LinkedList<>();
        for (int i = 0; i < threadCount; i++) {
            Thread thread = new ThreadViterbi(model, sequenceQueue, sequenceCount, viterbiPaths);
            threads.add(thread);
            thread.start();
        }

        // Waiting for threads to finish
        while (!threads.isEmpty()) {
            Thread thread = threads.poll();
            try {
                thread.join();
                Log.dLine(thread.getName() + " finished");
            } catch (InterruptedException e) {
                Log.eLine("ERROR: " + thread.getName() + " got interrupted");
            }
        }
        // all Threads finished

        return Arrays.asList(viterbiPaths);
    }
}
