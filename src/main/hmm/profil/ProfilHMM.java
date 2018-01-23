package main.hmm.profil;

import main.fastaparser.Sequence;
import main.hmm.HMM;
import main.logger.Log;

import java.util.LinkedList;
import java.util.List;

/**
 * <p><p>Aufgabe: </p>
 * Zum Einlesen des multiplen Sequenzalignments (MSA) der ribosomalen RNA
 * (rRNA) Trainingsequenzen (Datei: LSU train.fasta) koennen Sie sich an dem
 * C-Beispiel (im gleichen Daten Verzeichnis) orientieren oder dieses direkt in ihrem Programm verwenden.
 * Beim Bestimmen der Modellstruktur (Anzahl der Match-Zustaende) wenden Sie
 * die in der Vorlesung besprochene 50% Regel an: wenn in einer Spalte des MSA
 * weniger als 50% der Sequenzen Gaps aufweisen, dann wird diese Spalte einem
 * Match-Zustand zugeordnet. Danach koennen Sie die Emissionswahrscheinlichkeiten und Uebergangswahrscheinlichkeiten schaetzen.
 * Verwenden Sie zum Schaetzen aller Wahrscheinlichkeiten einen globalen Pseudocount-Parameter r = 1 (Laplace-Regel).
 * Lassen sie sich fuer die Match-Zustaende der Trainingssequenzen die Emissionswahrscheinlichkeiten ausgeben.
 * </p>
 * Anhand verschiedener Quellen implementiert.
 * <p>
 * Enthaelt Methode buildModel, die aus uebergebenen Trainings-Sequenzen ein ProfilHMM erstellt.
 * Enthaelt die Implementation des Viterbi-Algorithmus.
 * Dieser generiert aus einer uebergebenen Sequenz einen Zustands-Pfad. M = Match, D = Delete, I = Insert.
 *
 * @author Soeren Metje
 */
public class ProfilHMM {
    /*
    PROFIL HIDDEN MARKOV MODEL
               __              __                __                __                __                                  __                __
              (  )            (  )              (  )              (  )              (  )                                (  )              (  )
              /---\           /---\             /---\             /---\             /---\                               /---\             /---\
             /  I  \         /  I  \           /  I  \           /  I  \           /  I  \                             /  I  \           /  I  \
             \     /         \     /           \     /           \     /           \     /                             \     /           \     /
            />\---/\        />\---/\          />\---/\          />\---/\          />\---/\                            />\---/\          />\---/\
           /        \      /         \       /         \       /         \       /         \       /                           \       /         \
    |-----|         |------|          |------|          |------|          |------|          |------|                            |------|          |-----|
    |  B  |-------->|   M  |--------->|   M  |--------->|   M  |--------->|   M  |--------->|   M  |    ##   ##  ##   --------->|   M  |--------->|  E  |
    |-----|\        |------|\       />|------|\       />|------|\       />|------|\       />|------|    ##   ##  ##   \       />|------|        />|-----|
             \                \    /            \    /            \    /            \    /                              \    /                 /
               \                \/                \/                \/                \/                                  \/                 /
                 \     __     /    \     __     /    \     __     /    \     __     /    \     __                       /    \     __     /
                   \>/    \ /        \>/    \ /        \>/    \ /        \>/    \ /        \>/    \                   /        \>/    \ /
                    (  D   )          (  D   )          (  D   )          (  D   )          (  D   )                            (  D   )
                     \ __ /            \ __ /            \ __ /            \ __ /            \ __ /                              \ __ /

     */

    /**
     * Pseudo-Count fuer Berechnung der Emissions-Wahrscheinlichkeiten
     */
    private static final int PSEUDO_COUNT_EMISSION = 1;

    /**
     * Pseudo-Count fuer Berechnung der Uebergangs-Wahrscheinlichkeiten
     */
    private static final int PSEUDO_COUNT_TRANSITION = 1;

    /**
     * Anteil and Nukleotiden (also keine gaps), ab dem die Spalte als Match-State gezaehlt wird
     */
    private static final double THRESHOLD_MATCHSTATE = .5d; // min amount of nucleotides (no gap) to count column as match-state

    /**
     * Zeichen fuer Gap
     */
    private static final char GAP = '-';

    /**
     * Zeichen fuer Nukleotide
     */
    private static final char[] BASES = {'A', 'C', 'G', 'U'};

    /**
     * Match-Zustand
     */
    private static final char STATE_MATCH = 'M';

    /**
     * Insert-Zustand
     */
    private static final char STATE_INSERT = 'I';

    /**
     * Delete-Zustand
     */
    private static final char STATE_DELETE = 'D';

    /**
     * Zustand, wenn sich Gap in Insert-Spalte befindet
     */
    private static final char STATE_IGNORE = ' ';

    /**
     * Zustaende
     */
    private static final char[] STATES = {STATE_MATCH, STATE_INSERT, STATE_DELETE};

    /**
     * Anzahl der Zustaende
     */
    private static final int STATE_COUNT = STATES.length;

    /**
     * Beobachtungswahrscheinlichketen der Nukleotide im Match-Zustand an Position im Modell
     */
    private double[][] emissionProbMatch;

    /**
     * Beobachtungswahrscheinlichketen der Nukleotide im Insert-Zustand an Position im Modell
     */
    private double[][] emissionProbInsert;

    /**
     * Uebergangswahrscheinlichen zwischen den Zustaenden an Position im Modell
     */
    private double[][][] transitionProb;

    /**
     * Laenge des Modells bzw Anzahl der Match-Zustaende im Modell.
     * Der Start-Zustand wird auch als Match-Zustand interpretiert
     */
    private int lengthModel; // interpreting beginning-state als match-state

    /**
     * Konstruktor. Erstellt Modell und fuehrt die Methode buildModel aus.
     * Anschliessend wird in logarithmischen Raum konvertiert
     *
     * @param sequencesTrain Trainings-Sequenzen
     * @throws IllegalArgumentException falls in buildModel ein Fehler auftritt
     */
    public ProfilHMM(List<Sequence> sequencesTrain) throws IllegalArgumentException {

        buildModel(sequencesTrain);
        // convert into logspace (can be done before Viterbi-Algo is running)
        HMM.convertToLogspace(transitionProb);
        HMM.convertToLogspace(emissionProbMatch);
        HMM.convertToLogspace(emissionProbInsert);
    }

    /**
     * Extrahiert Daten fuer das ProfilHMM aus Trainings-Sequenzen.
     * Setzt Laenge des Modells.
     * Erstellt Felder fuer Beobachtungswahrscheinlichketen und Uebergangswahrscheinlichen zwischen den Zustaenden.
     *
     * @param sequencesTrain Trainings-Sequenzen
     * @throws IllegalArgumentException falls uebergebenes Feld == null
     *                                  oder das uebergebene Feld leer ist
     *                                  oder Sequenzen unterschiedlich lang
     *                                  oder Beobachtung nicht im Feld gefunden wird
     */
    private synchronized void buildModel(List<Sequence> sequencesTrain) throws IllegalArgumentException {
        Log.iLine("Building ProfilHMM -----------------------------");
        // checks
        if (sequencesTrain == null)
            throw new IllegalArgumentException("sequencesTrain is null");


        int seqenceCount = sequencesTrain.size();
        if (seqenceCount <= 0) {
            throw new IllegalArgumentException("sequencesTrain is empty");
        }

        int length = sequencesTrain.get(0).getNucleotideSequence().length();
        for (Sequence s : sequencesTrain) {
            int sLength = s.getNucleotideSequence().length();
            if (sLength != length) {
                throw new IllegalArgumentException("Sequence '" + s.getDescription()
                        + "' has different lenght (" + sLength + ") then the first Sequence (" + length + ")");
            }
        }

        Log.iLine("Sequence count = " + seqenceCount);
        Log.iLine("Sequence length = " + length);

        // find Match or Insertion-States and Model length --------------------------------------------------------
        // also count nucleotides and gaps ------------------------------------------------------------------------
        int[] gapCounts = new int[length];
        int[][] baseCounts = new int[length][BASES.length];

        for (Sequence seq : sequencesTrain) {
            for (int i = 0; i < length; i++) {
                char base = seq.getNucleotideSequence().charAt(i);
                if (base == GAP) {
                    gapCounts[i]++;
                } else {
                    baseCounts[i][observationToIndex(base)]++;
                }
            }
        }

        boolean[] matchState = new boolean[length];
        lengthModel = 1; // start-state is first match-state

        {
            int countThreshold = (int) (seqenceCount * (1d - THRESHOLD_MATCHSTATE));
            for (int i = 0; i < gapCounts.length; i++) {
                if (gapCounts[i] <= countThreshold) {
                    lengthModel++;
                    matchState[i] = true;
                }
            }
        }

        Log.iLine("Model length = " + lengthModel + " (interpreting start-state als match-state)");


        if (Log.isPrintDebug()) {
            // debug output gaps and nucleotide counts
            StringBuilder out = new StringBuilder(("Gap-Counts (m=Match-State \u001B[32mI\u001B[0m=Insert-State):\n"));
            for (int i = 0; i < length; i++) {
                out.append(matchState[i] ? "   m " : "\u001B[32m   I \u001B[0m");
            }
            out.append('\n');

            for (int i = 0; i < gapCounts.length; i++) {
                out.append(String.format("%4d ", gapCounts[i]));
            }
            out.append('\n');
            Log.dLine(out.toString());
            out = new StringBuilder(("Nucleotide-Counts: \n"));
            for (int i = 0; i < baseCounts[0].length; i++) {
                for (int j = 0; j < baseCounts.length; j++) {

                    out.append(String.format("%4d ", baseCounts[j][i]));
                }
                out.append('\n');
            }
            Log.dLine(out.toString());
        }


        // calc Emission Prob ---------------------------------------------------------------------------

        emissionProbMatch = new double[lengthModel][BASES.length]; // emission-probabilities for match-states
        emissionProbInsert = new double[lengthModel][BASES.length]; // emission-probabilities for insert-states


        int[] baseCountsInsert = new int[BASES.length];
        for (int i = 0, iModel = 1; i < length + 1; i++) {

            boolean end = i >= length;

            if (end || matchState[i]) {
                // set Insert-State Prob
                {
                    double[] emissionProbVector = emissionProbInsert[iModel - 1];

                    int sumEmissionCount = 0;

                    for (int baseCount : baseCountsInsert) {
                        sumEmissionCount += baseCount + PSEUDO_COUNT_EMISSION;
                    }

                    for (int j = 0; j < emissionProbVector.length; j++) {
                        emissionProbVector[j] = ((double) (PSEUDO_COUNT_EMISSION + baseCountsInsert[j])) / sumEmissionCount;
                    }
                }

            }

            if (!end) {
                if (matchState[i]) {
                    // set Match-State Prob
                    {
                        double[] emissionProbVector = emissionProbMatch[iModel];
                        int[] baseCountVector = baseCounts[i];

                        int sumEmissionCount = 0;

                        for (int baseCount : baseCountVector) {
                            sumEmissionCount += baseCount + PSEUDO_COUNT_EMISSION;
                        }

                        for (int j = 0; j < emissionProbVector.length; j++) {
                            emissionProbVector[j] = ((double) (PSEUDO_COUNT_EMISSION + baseCountVector[j])) / sumEmissionCount;
                        }
                    }

                    // reset Inset-State-Counts
                    for (int j = 0; j < baseCountsInsert.length; j++) {
                        baseCountsInsert[j] = 0;
                    }

                    iModel++;
                }
                // Insert-State
                else {
                    for (int j = 0; j < baseCountsInsert.length; j++) {
                        baseCountsInsert[j] += baseCounts[i][j];
                    }
                }
            }
        }


        // output Emission-Prob
        {
            StringBuilder outEimissionProb = new StringBuilder("\n");
            outEimissionProb.append("Emission Prob Match: (Pseudo-Count = " + PSEUDO_COUNT_EMISSION + ")\n");
            for (int i = 0; i < emissionProbMatch[0].length; i++) {
                outEimissionProb.append(String.format("%c:  ", BASES[i]));
                for (int j = 0; j < emissionProbMatch.length; j++) {
                    outEimissionProb.append(String.format("%.2f ", emissionProbMatch[j][i]));
                }
                outEimissionProb.append('\n');
            }
            Log.iLine(outEimissionProb.toString());

            if (Log.isPrintDebug()) {
                outEimissionProb = new StringBuilder();
                outEimissionProb.append("Emission Prob Insert: (Pseudo-Count = " + PSEUDO_COUNT_EMISSION + ")\n");
                for (int i = 0; i < emissionProbInsert[0].length; i++) {
                    outEimissionProb.append(String.format("%c:  ", BASES[i]));
                    for (int j = 0; j < emissionProbInsert.length; j++) {
                        outEimissionProb.append(String.format("%.2f ", emissionProbInsert[j][i]));
                    }
                    outEimissionProb.append('\n');
                }
                Log.dLine(outEimissionProb.toString());
            }
        }

        // calc Transition Count ----------------------------------------------------------------------

        if (Log.isPrintDebug()) {
            StringBuilder out = new StringBuilder("States and Transition Counts:\n");
            out.append("    ");
            for (int i = 0; i < length; i++) {
                out.append(matchState[i] ? "  m" : "\u001B[32m  I\u001B[0m");
            }
            out.append('\n');
            Log.d(out.toString());
        }

        int[][][] transitionCount = new int[STATES.length][STATES.length][lengthModel];
        {
            StringBuilder out = null;
            if (Log.isPrintDebug()) {
                out = new StringBuilder();
            }
            for (Sequence sequence : sequencesTrain) {
                String sequenceString = sequence.getNucleotideSequence();

                // output sequence
                if (out != null) {
                    out.append("\n\u001B[37m");
                    out.append(String.format("%.4s", sequence.getDescription()));
                    for (int i = 0; i < length; i++) {
                        out.append("  ").append(sequenceString.charAt(i));
                    }
                    out.append("\u001B[0m\n");
                }


                // get States and Transitions
                char lastState = STATE_MATCH; // interpreting Start-state as Match-state
                int insertCount = 0;

                if (out != null) {
                    out.append("    ");
                }
                for (int i = 0, iModel = 0; i < length + 1; i++) {
                    char state = getState(sequenceString, matchState, i);
                    if (state != STATE_IGNORE) {

                        if (state == STATE_INSERT && lastState == STATE_INSERT) {
                            insertCount++;
                        }
                        // no Insert-Insert (interpreting start and End as Match-state)
                        else {
                            transitionCount[stateToIndex(lastState)][stateToIndex(state)][iModel]++;

                            if (state != STATE_INSERT) {
                                if (lastState == STATE_INSERT) { // last Insert-Insert has ended
                                    transitionCount[stateToIndex(STATE_INSERT)][stateToIndex(lastState)][iModel] += insertCount;
                                    insertCount = 0;
                                }
                                iModel++;
                            }
                        }


                        lastState = state;
                    }
                    if (out != null) {
                        out.append("  ").append(state);
                    }
                }
                if (out != null) {
                    out.append('\n');
                }
            }
            if (out != null) {
                Log.dLine(out.toString());
            }
        }


        if (Log.isPrintDebug()) {
            // Debug output Transition Counts
            StringBuilder outTransCounts = new StringBuilder("\n");
            for (int i = 0; i < STATES.length; i++) {
                for (int j = 0; j < STATES.length; j++) {
                    outTransCounts.append(STATES[i]).append(STATES[j]).append(":  ");
                    for (int k = 0; k < transitionCount[i][j].length; k++) {
                        outTransCounts.append(String.format("%4d ", transitionCount[i][j][k]));
                    }
                    outTransCounts.append('\n');
                }
            }
            Log.dLine(outTransCounts.toString());
        }


        // calc Transition Prob ----------------------------------------------------------------------
        transitionProb = new double[STATES.length][STATES.length][lengthModel];
        for (int i = 0; i < lengthModel; i++) {

            int[] sumTrans = new int[STATES.length]; // vergleichsgroessen
            for (int j = 0; j < STATES.length; j++) {
                for (int k = 0; k < STATES.length; k++) {
                    sumTrans[j] += transitionCount[j][k][i] + PSEUDO_COUNT_TRANSITION;
                }
            }

            for (int j = 0; j < STATES.length; j++) {
                for (int s = 0; s < STATES.length; s++) {
                    transitionProb[j][s][i] = ((double) transitionCount[j][s][i] + PSEUDO_COUNT_TRANSITION) / sumTrans[j];
                }
            }
        }

        // output Transition Prob
        if (Log.isPrintDebug()) {
            StringBuilder outTransProb = new StringBuilder("\nTransition Prob: (Pseudo-Count = " + PSEUDO_COUNT_TRANSITION + ")\n");
            for (int i = 0; i < transitionProb.length; i++) {
                double[][] transProbMatrix = transitionProb[i];
                for (int j = 0; j < transProbMatrix.length; j++) {
                    double[] transProbVector = transProbMatrix[j];
                    outTransProb.append(String.format("%s:  ", STATES[i] + "" + STATES[j]));
                    for (double prob : transProbVector) {
                        outTransProb.append(String.format("%.2f ", prob));
                    }
                    outTransProb.append('\n');
                }
            }
            Log.dLine(outTransProb.toString());
        }
    }

    /**
     * Implementation des Viterbi-Algorithmus fuer den logarithmischen Raum
     *
     * @param sequence Beobachtungsfolge
     * @return Zustands-Pfad
     * @throws IllegalArgumentException falls uebergebenes Feld == null
     *                                  oder Beobachtung nicht im Feld gefunden wird
     */
    public ViterbiPath viterbi(final Sequence sequence) throws IllegalArgumentException {
        if (sequence == null)
            throw new IllegalArgumentException("sequence is null");

        // init
        char[] observations = sequence.getNucleotideSequence().toCharArray();
        int[] observationIndices = observationsToIndices(observations);
        int length = observationIndices.length + 1;

        // FILL MATRIX ----------------------------------------------------------------------------------
        double[][][] viterbiVar = new double[STATE_COUNT][length][lengthModel];
        int[][][] viterbiArg = new int[STATE_COUNT][length][lengthModel];

        // init
        {
            double initValue = Double.NEGATIVE_INFINITY;

            int stateIndex = stateToIndex(STATE_MATCH);
            viterbiVar[stateIndex][0][0] = 0d;
            viterbiArg[stateIndex][0][0] = -1;
            for (int j = 1; j < lengthModel; j++) {
                viterbiVar[stateIndex][0][j] = initValue;
            }
            for (int i = 1; i < length; i++) {
                viterbiVar[stateIndex][i][0] = initValue;
            }

            stateIndex = stateToIndex(STATE_INSERT);
            for (int j = 0; j < lengthModel; j++) {
                viterbiVar[stateIndex][0][j] = initValue;
            }

            stateIndex = stateToIndex(STATE_DELETE);
            for (int i = 0; i < length; i++) {
                viterbiVar[stateIndex][i][0] = initValue;
            }
        }

        // iterate observations indices
        for (int i = 0; i < length; i++) {
            // iterate model indices
            for (int j = 0; j < lengthModel; j++) {
                // iterate states indices
                for (int s = 0; s < STATE_COUNT; s++) { // order of iteration-loops is relevant!
                    char state = STATES[s];

                    int iShift = i, jShift = j;
                    double[][] emissionProbMatrix = null;
                    if (state == STATE_MATCH) {
                        iShift -= 1;
                        jShift -= 1;
                        emissionProbMatrix = emissionProbMatch;
                    } else if (state == STATE_INSERT) {
                        iShift -= 1;
                        emissionProbMatrix = emissionProbInsert;
                    } else if (state == STATE_DELETE) {
                        jShift -= 1;
                    } else
                        throw new RuntimeException("no valid state");

                    // calc must be possible for Delete-State in first column ans Insert-State in first row
                    if (iShift >= 0 && jShift >= 0) {
                        //find max
                        double maxProb = Double.NEGATIVE_INFINITY;
                        int maxArg = -1; // maximizing argument

                        for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                            double prob = viterbiVar[stateIndex][iShift][jShift] + transitionProb[stateIndex][s][jShift]; // log-space
                            if (prob > maxProb) {
                                maxProb = prob;
                                maxArg = stateIndex;
                            }
                        }

                        double emissionProb = 0d; // 0 is neutral element of addition (log-space)

                        if (emissionProbMatrix != null) {
                            emissionProb = emissionProbMatrix[j][observationIndices[i - 1]];
                        }
                        viterbiVar[s][i][j] = emissionProb + maxProb;
                        viterbiArg[s][i][j] = maxArg;
                    }
                }
            }
        }

        if (Log.isPrintDebug()) {
            // Debug output viterbi 3d-matrix
            StringBuilder outViterbiVar = new StringBuilder("ViterbiVar: \n");
            // [STATE_COUNT][length][lengthModel]
            for (int j = 0; j < lengthModel; j++) {
                for (int k = 0; k < STATES.length; k++) {
                    outViterbiVar.append("\u001B[37m").append(k == 0 ? String.format("j%3d%s ", j, STATES[k]) : "    " + STATES[k] + " ").append("\u001B[0m");
                    for (int i = 0; i < length; i++) {
                        outViterbiVar.append(String.format("%.5s ", String.format("%f", viterbiVar[k][i][j])));
                    }
                    outViterbiVar.append('\n');
                }
                outViterbiVar.append('\n');
            }
            Log.dLine(outViterbiVar.toString());

            StringBuilder outViterbiArg = new StringBuilder("ViterbiArg: \n");
            // [STATE_COUNT][length][lengthModel]
            for (int j = 0; j < lengthModel; j++) {
                for (int k = 0; k < STATES.length; k++) {
                    outViterbiArg.append("\u001B[37m").append(k == 0 ? String.format("j%3d%s ", j, STATES[k]) : "    " + STATES[k] + " ").append("\u001B[0m");
                    for (int i = 0; i < length; i++) {
                        int a = viterbiArg[k][i][j];
                        outViterbiArg.append(String.format("%5s ", (a >= 0 ? String.valueOf(STATES[a]) : a)));
                    }
                    outViterbiArg.append('\n');
                }
                outViterbiArg.append('\n');
            }
            Log.dLine(outViterbiArg.toString());
        }

        // BACKTRACE -------------------------------------------------------------------------------------
        List<Character> listStatePath = new LinkedList<>();
        double score = Double.NEGATIVE_INFINITY;
        {
            // backtrace init / Find path with max prob
            int i = length - 1, j = lengthModel - 1;
            int stateIndexEnd = -1;
            {
                for (int stateIndex = 0; stateIndex < STATE_COUNT; stateIndex++) {

                    double prob = viterbiVar[stateIndex][i][j] + transitionProb[stateIndex][stateToIndex(STATE_MATCH)][j]; // log-space
                    if (prob > score) {
                        stateIndexEnd = stateIndex;
                        score = prob;
                    }
                }
            }
            viterbiVar = null; // no reference left -> allow GC to trash

            listStatePath.add(STATES[stateIndexEnd]);

            // backtrace iterate
            try {
                while (i > 0 && j > 0 && (i != 1 || j != 1)) { // FIXME richtig?
                    int stateIndex = viterbiArg[stateIndexEnd][i][j];
                    char state = STATES[stateIndex];

                    listStatePath.add(0, state);

                    if (state == STATE_MATCH) {
                        i--;
                        j--;
                    } else if (state == STATE_INSERT) {
                        i--;
                    } else if (state == STATE_DELETE) {
                        j--;
                    }
                }

            } catch (ArrayIndexOutOfBoundsException e) {
                throw new ArrayIndexOutOfBoundsException(e.getMessage() + " i=" + i + " j=" + j + " seq=" + sequence.getDescription());
            }
        }

        // copy into char array
        char[] statePath = new char[listStatePath.size()];
        {
            int i = 0;
            for (Character c : listStatePath) {
                statePath[i] = c;
                i++;
            }
        }

        return new ViterbiPath(sequence, score, statePath);
    }

    /**
     * gibt den Zustand zurueck, in dem sich das HMM an uebergebenem index in uebergebener Sequenz befindet
     * oder ein Leerzeichen, falls sich das Modell an der Stelle im Insert-Zustand befindet, aber kein Zeichen in der Sequenz vorhanden ist.
     * Das uebergebene Feld matchState muss auskunft darueber geben, ob sich das Modell an uebergebenem index im Insert- oder Match-Zustand befindet (true falls Match-Zustand).
     *
     * @param seq        Sequenz
     * @param matchState Insert- oder Match-Zustaende
     * @param index      position Sequenz
     * @return Zustand
     */
    private static char getState(final String seq, final boolean[] matchState, final int index) {

        if (index >= seq.length())
            return STATE_MATCH; // ende

        boolean match = matchState[index];

        if (!match) {
            if (seq.charAt(index) != GAP) {
                return STATE_INSERT;
            } else {
                return STATE_IGNORE;
            }
        }

        // at i is a Match-state
        if (seq.charAt(index) == GAP)
            return STATE_DELETE;
        return STATE_MATCH;
    }

    /**
     * mappt Beaobachtung-Folge auf entsprechende Index-Folge
     *
     * @param observations Beobachtungs-Folge
     * @return entsprechende Index-Folge
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    private static int[] observationsToIndices(final char[] observations) throws IllegalArgumentException {
        return HMM.charsToIndices(BASES, observations);
    }

    /**
     * mappt Beaobachtung auf entsprechenden Index
     *
     * @param observation Beobachtungs
     * @return entsprechender Index
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    private static int observationToIndex(final char observation) throws IllegalArgumentException {
        return HMM.charToIndex(BASES, observation);
    }

    /**
     * mappt Zustand auf entsprechenden Index
     *
     * @param state Zustand
     * @return entsprechender Index
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    private static int stateToIndex(final char state) throws IllegalArgumentException {
        return HMM.charToIndex(STATES, state);
    }
}