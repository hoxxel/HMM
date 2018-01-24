package main.fastaparser;

/**
 * Eingelesene Nukleotid-Sequenz aus {@link FastaParser}.
 * <p>
 * Enthaelt:
 * <p>
 * - Beschreibung, durch &gt; kenntlich gemacht
 * <p>
 * - Kommentar (optional, wenn in Datei angegeben), durch ; kenntlich gemacht
 * <p>
 * - Nukleotid-Sequenz, folgende Zeile ohne initiierende Zeichen
 *
 * @author Soeren Metje
 */
public class Sequence {

    /**
     * Beschreibung, durch &gt; kenntlich gemacht
     */
    private final String description;

    /**
     * Kommentar (optional, wenn in Datei angegeben), durch ; kenntlich gemacht.
     */
    private final String comments;

    /**
     * Nukleotid-Sequenz, folgende Zeile ohne initiierende Zeichen
     */
    private final String nucleotideSequence;

    /**
     * Konstruktor
     *
     * @param description        Beschreibung
     * @param comments           Kommentar
     * @param nucleotideSequence Nukleotid-Sequenz
     */
    public Sequence(String description, String comments, String nucleotideSequence) {
        this.description = description;
        this.comments = comments;
        this.nucleotideSequence = nucleotideSequence;
    }

    /**
     * liefert Beschreibung zurueck
     *
     * @return Beschreibung
     */
    public String getDescription() {
        return description;
    }

    /**
     * liefert Kommentar zurueck
     *
     * @return Kommentar
     */
    public String getComments() {
        return comments;
    }

    /**
     * liefert Nukleotid-Sequence zurueck
     *
     * @return Nukleotid-Sequence
     */
    public String getNucleotideSequence() {
        return nucleotideSequence;
    }

    /**
     * liefert true zurueck, falls Kommentar vorhanden ist. Ansonsten false
     *
     * @return true, falls Kommentar vorhanden ist. Ansonsten false
     */
    public boolean hasComments() {
        return comments != null;
    }
}
