package main.hmm.profil;

import main.fastaparser.Sequence;
import main.hmm.HMMFunc;
import main.logger.Log;

import java.util.List;

/**
 * <p>
 * Profil Hidden Markov Model.
 * </p>
 * Enthaelt Methode buildModel, die aus uebergebenen Trainings-Sequenzen ein RNAProfilHMM erstellt.
 *
 * @author Soeren Metje
 */
public abstract class ProfilHMM {
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
     * Match-Zustand
     */
    public static final char STATE_MATCH = 'M';

    /**
     * Insert-Zustand
     */
    public static final char STATE_INSERT = 'I';

    /**
     * Delete-Zustand
     */
    public static final char STATE_DELETE = 'D';

    /**
     * Zustand, wenn sich Gap in Insert-Spalte befindet
     */
    public static final char STATE_IGNORE = ' ';

    /**
     * Zustaende
     */
    public static final char[] STATES = {STATE_MATCH, STATE_INSERT, STATE_DELETE};

    /**
     * Anzahl der Zustaende
     */
    public static final int STATE_COUNT = STATES.length;

    /**
     * Index des Match-Zustand
     */
    public static final int STATE_MATCH_INDEX = stateToIndex(STATE_MATCH);

    /**
     * Index des Insert-Zustand
     */
    public static final int STATE_INSERT_INDEX = stateToIndex(STATE_INSERT);

    /**
     * Index des Delete-Zustand
     */
    public static final int STATE_DELETE_INDEX = stateToIndex(STATE_DELETE);

    // Instanz-Variablen ###############################################################################################

    /**
     * Pseudo-Count fuer Berechnung der Emissions-Wahrscheinlichkeiten
     */
    private final int pseudoCountEmission;

    /**
     * Pseudo-Count fuer Berechnung der Uebergangs-Wahrscheinlichkeiten
     */
    private final int pseudoCountTransition;

    /**
     * Anteil and Nukleotiden (also keine gaps), ab dem die Spalte als Match-State gezaehlt wird
     */
    private final double thresholdMatchState; // min amount of nucleotides (no gap) to count column as match-state

    /**
     * Zeichen fuer Gap
     */
    private final char gap;

    /**
     * Zeichen fuer Nukleotide
     */
    private final char[] bases;

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
     * Laenge des Modells bzw. Anzahl der Match-Zustaende im Modell.
     * Der Start-Zustand wird auch als Match-Zustand interpretiert
     */
    private int lengthModel; // interpreting beginning-state als match-state

    /**
     * Konstruktor. Erstellt Modell und fuehrt die Methode buildModel aus.
     * Anschliessend werden die logarithmierten Wahrscheinlichkeiten berechnet.
     *
     * @param sequencesTrain        Trainings-Sequenzen
     * @param gap                   Zeichen fuer Gap
     * @param bases                 Zeichen fuer Nukleotide
     * @param pseudoCountEmission   Pseudo-Count fuer Berechnung der Emissions-Wahrscheinlichkeiten
     * @param pseudoCountTransition Pseudo-Count fuer Berechnung der Uebergangs-Wahrscheinlichkeiten
     * @param thresholdMatchState   Anteil and Nukleotiden (also keine gaps), ab dem die Spalte als Match-State gezaehlt wird
     * @throws IllegalArgumentException falls in buildModel ein Fehler auftritt
     */
    public ProfilHMM(List<Sequence> sequencesTrain, char gap, char[] bases, int pseudoCountEmission, int pseudoCountTransition, double thresholdMatchState) throws IllegalArgumentException {
        this.gap = gap;
        this.bases = bases;
        this.pseudoCountEmission = pseudoCountEmission;
        this.pseudoCountTransition = pseudoCountTransition;
        this.thresholdMatchState = thresholdMatchState;
        buildModel(sequencesTrain);
        // calc log for each element in all matrices (can be done before Viterbi-Algo is running)
        HMMFunc.logspace(transitionProb);
        HMMFunc.logspace(emissionProbMatch);
        HMMFunc.logspace(emissionProbInsert);
    }

    /**
     * Extrahiert Daten fuer das RNAProfilHMM aus Trainings-Sequenzen.
     * Setzt Laenge des Modells.
     * Erstellt Felder fuer Beobachtungswahrscheinlichketen und Uebergangswahrscheinlichen zwischen den Zustaenden.
     *
     * @param sequencesTrain Trainings-Sequenzen
     * @throws IllegalArgumentException falls uebergebenes Feld == null
     *                                  oder das uebergebene Feld leer ist
     *                                  oder Sequenzen unterschiedlich lang
     *                                  oder Beobachtung nicht im Feld gefunden wird
     */
    private void buildModel(List<Sequence> sequencesTrain) throws IllegalArgumentException {
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
        int[][] baseCounts = new int[length][bases.length];

        for (Sequence seq : sequencesTrain) {
            for (int i = 0; i < length; i++) {
                char base = seq.getNucleotideSequence().charAt(i);
                if (base == gap) {
                    gapCounts[i]++;
                } else {
                    baseCounts[i][observationToIndex(base)]++;
                }
            }
        }

        boolean[] matchState = new boolean[length];
        lengthModel = 1; // start-state is first match-state

        {
            int countThreshold = (int) (seqenceCount * (1d - this.thresholdMatchState));
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

        emissionProbMatch = new double[lengthModel][bases.length]; // emission-probabilities for match-states
        emissionProbInsert = new double[lengthModel][bases.length]; // emission-probabilities for insert-states


        int[] baseCountsInsert = new int[bases.length];
        for (int i = 0, iModel = 1; i < length + 1; i++) {

            boolean end = i >= length;

            if (end || matchState[i]) {
                // set Insert-State Prob
                {
                    double[] emissionProbVector = emissionProbInsert[iModel - 1];

                    int sumEmissionCount = 0;

                    for (int baseCount : baseCountsInsert) {
                        sumEmissionCount += baseCount + pseudoCountEmission;
                    }

                    for (int j = 0; j < emissionProbVector.length; j++) {
                        emissionProbVector[j] = ((double) (pseudoCountEmission + baseCountsInsert[j])) / sumEmissionCount;
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
                            sumEmissionCount += baseCount + pseudoCountEmission;
                        }

                        for (int j = 0; j < emissionProbVector.length; j++) {
                            emissionProbVector[j] = ((double) (pseudoCountEmission + baseCountVector[j])) / sumEmissionCount;
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
            outEimissionProb.append("Emission Prob Match: (Pseudo-Count = " + pseudoCountEmission + ")\n");
            for (int i = 0; i < emissionProbMatch[0].length; i++) {
                outEimissionProb.append(String.format("%c:  ", bases[i]));
                for (int j = 0; j < emissionProbMatch.length; j++) {
                    outEimissionProb.append(String.format("%.2f ", emissionProbMatch[j][i]));
                }
                outEimissionProb.append('\n');
            }
            Log.iLine(outEimissionProb.toString());

            if (Log.isPrintDebug()) {
                outEimissionProb = new StringBuilder();
                outEimissionProb.append("Emission Prob Insert: (Pseudo-Count = " + pseudoCountEmission + ")\n");
                for (int i = 0; i < emissionProbInsert[0].length; i++) {
                    outEimissionProb.append(String.format("%c:  ", bases[i]));
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
                    sumTrans[j] += transitionCount[j][k][i] + pseudoCountTransition;
                }
            }

            for (int j = 0; j < STATES.length; j++) {
                for (int s = 0; s < STATES.length; s++) {
                    transitionProb[j][s][i] = ((double) transitionCount[j][s][i] + pseudoCountTransition) / sumTrans[j];
                }
            }
        }

        // output Transition Prob
        if (Log.isPrintDebug()) {
            StringBuilder outTransProb = new StringBuilder("\nTransition Prob: (Pseudo-Count = " + pseudoCountTransition + ")\n");
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
     * Gibt den Zustand zurueck, in dem sich das HMM an uebergebenem index in uebergebener Sequenz befindet
     * oder ein Leerzeichen, falls sich das Modell an der Stelle im Insert-Zustand befindet, aber kein Zeichen in der Sequenz vorhanden ist.
     * Das uebergebene Feld matchState muss auskunft darueber geben, ob sich das Modell an uebergebenem index im Insert- oder Match-Zustand befindet (true falls Match-Zustand).
     *
     * @param seq        Sequenz
     * @param matchState Insert- oder Match-Zustaende
     * @param index      position Sequenz
     * @return Zustand
     */
    private char getState(final String seq, final boolean[] matchState, final int index) {

        if (index >= seq.length())
            return STATE_MATCH; // ende

        boolean match = matchState[index];

        if (!match) {
            if (seq.charAt(index) != gap) {
                return STATE_INSERT;
            } else {
                return STATE_IGNORE;
            }
        }

        // at i is a Match-state
        if (seq.charAt(index) == gap)
            return STATE_DELETE;
        return STATE_MATCH;
    }

    /**
     * Mappt Beaobachtung-Folge auf entsprechende Index-Folge
     *
     * @param observations Beobachtungs-Folge
     * @return entsprechende Index-Folge
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    public int[] observationsToIndices(final char[] observations) throws IllegalArgumentException {
        return HMMFunc.charsToIndices(bases, observations);
    }

    /**
     * Mappt Beaobachtung auf entsprechenden Index
     *
     * @param observation Beobachtungs
     * @return entsprechender Index
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    public int observationToIndex(final char observation) throws IllegalArgumentException {
        return HMMFunc.charToIndex(bases, observation);
    }

    public int getPseudoCountEmission() {
        return pseudoCountEmission;
    }

    public int getPseudoCountTransition() {
        return pseudoCountTransition;
    }

    public double getThresholdMatchState() {
        return thresholdMatchState;
    }

    public char getGap() {
        return gap;
    }

    public char[] getBases() {
        return bases;
    }

    public double[][] getEmissionProbMatch() {
        return emissionProbMatch;
    }

    public double[][] getEmissionProbInsert() {
        return emissionProbInsert;
    }

    public double[][][] getTransitionProb() {
        return transitionProb;
    }

    public int getLengthModel() {
        return lengthModel;
    }

    /**
     * Mappt Zustand auf entsprechenden Index
     *
     * @param state Zustand
     * @return entsprechender Index
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    private static int stateToIndex(final char state) throws IllegalArgumentException {
        return HMMFunc.charToIndex(STATES, state);
    }
}