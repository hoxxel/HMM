package main.hmm.profil;

import main.fastaparser.Sequence;

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
 * Enthaelt Methode buildModel, die aus uebergebenen Trainings-Sequenzen ein RNAProfilHMM erstellt.
 * Enthaelt die Implementation des Viterbi-Algorithmus fuer logarithmische Werte.
 * Dieser generiert aus einer uebergebenen Sequenz einen Zustands-Pfad. M = Match, D = Delete, I = Insert.
 * </p>
 *
 * @author Soeren Metje
 */
public class RNAProfilHMM extends ProfilHMM {
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
     * Konstruktor. Erstellt Modell und fuehrt die Methode buildModel aus.
     * Anschliessend werden die logarithmierten Wahrscheinlichkeiten berechnet.
     *
     * @param sequencesTrain Trainings-Sequenzen
     * @throws IllegalArgumentException falls in buildModel ein Fehler auftritt
     */
    public RNAProfilHMM(List<Sequence> sequencesTrain) throws IllegalArgumentException {
        super(sequencesTrain, GAP, BASES, PSEUDO_COUNT_EMISSION, PSEUDO_COUNT_TRANSITION, THRESHOLD_MATCHSTATE);
    }
}