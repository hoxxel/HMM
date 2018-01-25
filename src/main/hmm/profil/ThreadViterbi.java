package main.hmm.profil;

import main.fastaparser.Sequence;
import main.logger.Log;

import java.util.Queue;

/**
 * Thread {@link Thread}, der Sequnzen {@link Sequence} aus uebergebener Schlange abarbeitet.
 * Dabei wird fuer jede Sequenz anhand des uebergebenen Modells mittels des Viterbi-Algorithmus der maximierende Zustands-Pfad berechnet.
 * Der Score und der Zustands-Pfad der Sequenz wird dann mittles der Wrapper-Klasse {@link ViterbiPath}
 * zum uebergebenen Array an entsprechender Position hinzugefuegt.
 *
 * @author Soeren Metje
 */
public class ThreadViterbi extends Thread {
    /**
     * Monitor, um Ausgabe zu synchronisieren
     */
    private static final Object outputMonitor = new Object();

    /**
     * Monitor, um die Bestimmung des Index der naechsten Sequenz zu synchronisieren
     */
    private static final Object indexPollMonitor = new Object();

    /**
     * ProfilHMM, welches zur Berechnung verwendet wird
     */
    private final ProfilHMM model;

    /**
     * Schlange abzuarbeitender Sequenzen
     */
    private final Queue<Sequence> sequenceQueue;

    /**
     * Initiale Anzahl der Sequnzen in Schlange
     */
    private final int sequenzeQueueInitSize;

    /**
     * Threadsichere Liste zu der Ergebnisse hinzugefuegt werden
     */
    private final ViterbiPath[] finishedPaths;

    /**
     * Konstruktor
     *
     * @param model         zu verwendenes ProfilHMM
     * @param sequenceQueue abzuarbeitende Sequenzen
     * @param finishedPaths threadsichere Liste fuer Ergebnisse
     */
    public ThreadViterbi(ProfilHMM model, Queue<Sequence> sequenceQueue, int sequenzeCount, ViterbiPath[] finishedPaths) {
        this.model = model;
        this.sequenceQueue = sequenceQueue;
        this.finishedPaths = finishedPaths;
        this.sequenzeQueueInitSize = sequenzeCount;
    }

    /**
     * Arbeitet Sequnzen {@link Sequence} aus uebergebener Schlange ab.
     * Dabei wird fuer jede Sequenz anhand des uebergebenen Modells mittels des Viterbi-Algorithmus der maximierende Zustands-Pfad berechnet.
     * Der Score und der Zustands-Pfad der Sequenz wird dann mittles der Wrapper-Klasse {@link ViterbiPath}
     * zum uebergebenen Array an entsprechender Position hinzugefuegt.
     */
    @Override
    public void run() {
        Log.dLine(getName() + " started");
        boolean running = true;
        while (running) {

            Sequence sequence = null;
            int index = 0;
            synchronized (indexPollMonitor) {
                running = !sequenceQueue.isEmpty();
                if (running) {
                    index = sequenzeQueueInitSize - sequenceQueue.size();
                    sequence = sequenceQueue.poll();
                }
            }

            if (sequence != null) {
                ViterbiPath viterbiPath = null;

                long millis = System.currentTimeMillis(); // measure calc time
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
                    Log.iLine(String.valueOf(viterbiPath.getStatePath()));
                    Log.iLine();
                }

                finishedPaths[index] = viterbiPath;
            }
        }
    }
}
