package main.hmm.profil;

import main.fastaparser.Sequence;
import main.logger.Log;

import java.util.List;
import java.util.Queue;

public class ThreadViterbi extends Thread {
    private final ProfilHMM model;
    private final Queue<Sequence> sequenceQueue;
    private final List<ViterbiPath> finishedPaths;

    public ThreadViterbi(ProfilHMM model, Queue<Sequence> sequenceQueue, List<ViterbiPath> addPathsTo) {
        this.model = model;
        this.sequenceQueue = sequenceQueue;
        this.finishedPaths = addPathsTo;
    }

    @Override
    public void run() {
        Log.dLine(getName() + " started");
        while (!sequenceQueue.isEmpty()) {

            Sequence sequence = sequenceQueue.poll();

            ViterbiPath viterbiPath = null;
            try {
                viterbiPath = model.viterbi(sequence);
            } catch (IllegalArgumentException e) {
                Log.eLine("ERROR: Viterbi ProfilHMM failed! " + e.getMessage());
                System.exit(1);
            } catch (OutOfMemoryError e) {
                Log.eLine("ERROR: Out of Memory " + e.getMessage() + ". Start with more Memory. (Argument -Xmx<Size>)");
                System.exit(1);
            }

            Log.iLine("Sequece " + sequence.getDescription() + " -----------------------------");
            Log.iLine(sequence.getSequence());
            Log.iLine(viterbiPath.toString());
            Log.iLine();

            finishedPaths.add(viterbiPath);
        }
    }
}
