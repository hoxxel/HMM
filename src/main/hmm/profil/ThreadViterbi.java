package main.hmm.profil;

import main.fastaparser.Sequence;
import main.logger.Log;

import java.util.List;
import java.util.Queue;

/**
 * {@link Thread}, der Sequnzen {@link Sequence} aus uebergebener Schlange abarbeitet.
 * Dabei wird fuer jede Sequenz anhand des uebergebenen Modells mittels des Viterbi-Algorithmus der maximierende Zustands-Pfad berechnet.
 * Der Score und der Zustands-Pfad der Sequenz wird dann mittles der Wrapper-Klasse {@link ViterbiPath} zur uebergebenen threadsicheren Liste hinzugefuegt.
 *
 * @author Soeren Metje
 */
public class ThreadViterbi extends Thread {
    /**
     * Monitor, um Ausgabe zu synchronisieren
     */
    private static final Object outputMonitor = new Object();

    /**
     * ProfilHMM, welches zur Berechnung verwendet wird
     */
    private final ProfilHMM model;

    /**
     * Schlange abzuarbeitender Sequenzen
     */
    private final Queue<Sequence> sequenceQueue;

    /**
     * threadsichere Liste zu der Ergebnisse hinzugefuegt werden
     */
    private final List<ViterbiPath> finishedPaths;

    /**
     * Konstruktor
     *
     * @param model                   zu verwendenes ProfilHMM
     * @param sequenceQueue           abzuarbeitende Sequenzen
     * @param threadSaveFinishedPaths threadsichere Liste fuer Ergebnisse
     */
    public ThreadViterbi(ProfilHMM model, Queue<Sequence> sequenceQueue, List<ViterbiPath> threadSaveFinishedPaths) {
        this.model = model;
        this.sequenceQueue = sequenceQueue;
        this.finishedPaths = threadSaveFinishedPaths;
    }

    /**
     * Arbeitet Sequnzen {@link Sequence} aus uebergebener Schlange ab.
     * Dabei wird fuer jede Sequenz anhand des uebergebenen Modells mittels des Viterbi-Algorithmus der maximierende Zustands-Pfad berechnet.
     * Der Score und der Zustands-Pfad der Sequenz wird dann mittles der Wrapper-Klasse {@link ViterbiPath} zur uebergebenen threadsicheren Liste hinzugefuegt.
     */
    @Override
    public void run() {
        Log.dLine(getName() + " started");
        while (!sequenceQueue.isEmpty()) {

            Sequence sequence = sequenceQueue.poll();

            ViterbiPath viterbiPath = null;

            long millis = System.currentTimeMillis();
            try {
                viterbiPath = model.viterbi(sequence);
            } catch (IllegalArgumentException e) {
                Log.eLine("ERROR: Viterbi ProfilHMM failed! " + e.getMessage());
                System.exit(1);
            } catch (OutOfMemoryError e) {
                Log.eLine("ERROR: Out of Memory " + e.getMessage() + ". Start with more Memory. (Argument -Xmx<Size>)");
                System.exit(1);
            }
            millis = (System.currentTimeMillis() - millis); // calc time of viterbi
            float time = (float) millis / 1000; // in sec

            synchronized (outputMonitor) {
                Log.iLine(String.format("(%.2fsec) %s -----------------------------", time, sequence.getDescription()));
                Log.iLine(sequence.getNucleotideSequence());
                Log.iLine(viterbiPath.toString());
                Log.iLine();
            }

            finishedPaths.add(viterbiPath);
        }
    }
}
