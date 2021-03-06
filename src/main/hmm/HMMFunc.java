package main.hmm;

import main.hmm.profil.RNAProfilHMM;

/**
 * Beinhaltet Funktionen, die sowohl in {@link main.hmm.casino.CasinoHMM}
 * als auch in {@link RNAProfilHMM} benoetigt werden
 *
 * @author Soeren Metje
 */
public class HMMFunc {

    /**
     * Berechne und speichere den Logarithmus jedes Elements der uebergebenen 3d-Matrix als Nebeneffekt
     *
     * @param matrix Matrix
     */
    public static void logspace(final double[][][] matrix) {
        for (double[][] matrixTwo : matrix) {
            logspace(matrixTwo);
        }
    }

    /**
     * Berechne und speichere den Logarithmus jedes Elements der uebergebenen Matrix als Nebeneffekt
     *
     * @param matrix Matrix
     */
    public static void logspace(final double[][] matrix) {
        for (double[] vector : matrix) {
            logspace(vector);
        }
    }

    /**
     * Berechne und speichere den Logarithmus jedes Elements des uebergebenen Vektors als Nebeneffekt
     *
     * @param vector Vektor
     */
    public static void logspace(final double[] vector) {
        for (int j = 0; j < vector.length; j++) {
            vector[j] = Math.log(vector[j]);
        }
    }

    /**
     * Mappt Beaobachtung-Folge auf entsprechende Index-Folge
     *
     * @param space        Feld aller Beobachtungen
     * @param observations Beobachtungs-Folge
     * @return entsprechende Index-Folge
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    public static int[] charsToIndices(final char[] space, final char[] observations) throws IllegalArgumentException {
        int length = observations.length;
        int[] ret = new int[length];

        for (int i = 0; i < length; i++) {
            ret[i] = charToIndex(space, observations[i]);
        }
        return ret;
    }

    /**
     * Mappt Beaobachtung auf entsprechenden Index
     *
     * @param space       Feld aller Beobachtungen
     * @param observation Beobachtung
     * @return entsprechender Index
     * @throws IllegalArgumentException falls Beobachtung nicht im Feld gefunden wird
     */
    public static int charToIndex(final char[] space, final char observation) throws IllegalArgumentException {
        for (int i = 0, basesLength = space.length; i < basesLength; i++) {
            if (space[i] == observation)
                return i;
        }
        throw new IllegalArgumentException("Character " + observation + " not found");
    }
}
