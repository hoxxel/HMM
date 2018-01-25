package main.hmm.profil;

import main.fastaparser.Sequence;
import main.logger.Log;

import java.util.*;

/**
 * Enthaelt Methode zur parallelisierten Ausfuehrung des Viterbi-Algorithmus.
 */
public class ParallelizationSupporter {

    /**
     * Fuehrt Viterbi-Algorithmus parallelisiert aus und liefert die berechneten Zustands-Pfade {@link ViterbiPath} zurueck.
     * <p>
     * Es werden so viele Threads erstellt, wie logische Kerne der JVM zur verfuegung stehen.
     * Danach werden die Sequenzen auf die Threads moeglichst gleichmaessig verteilt.
     * Die Threads berechnen anhand des uebergebenen Models fuer jede Sequenz den Zustands-Pfad.
     * Abschlissend wird eine Liste mit den Zustands-Pfaden {@link ViterbiPath} zurueck geliefert.
     *
     * @param model     {@link ProfilHMM} Modell
     * @param sequences {@link Sequence} Sequenz
     * @return Liste mit den Zustands-Pfaden {@link ViterbiPath}
     */
    public static List<ViterbiPath> viterbiParallelized(ProfilHMM model, List<Sequence> sequences) {
        int sequenceCount = sequences.size();
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
            for (Sequence sequence : sequences) {
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
                Log.eLine("ERROR: " + thread.getName() + " got interrupted");
            }
        }
        // all Threads finished
        return viterbiPaths;
    }
}
