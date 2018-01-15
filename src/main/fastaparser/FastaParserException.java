package main.fastaparser;

/**
 * exception, die bei Parse-fehlern durch {@link FastaParser} ausgeloesst wird.
 *
 * @author Soeren Metje
 */
public class FastaParserException extends Exception {
    /**
     * erstellt exception
     */
    public FastaParserException() {

    }

    /**
     * erstellt exception
     *
     * @param msg info
     */
    public FastaParserException(String msg) {
        super(msg);
    }

    /**
     * erstellt exception
     *
     * @param msg   info
     * @param cause grund
     */
    public FastaParserException(String msg, Throwable cause) {
        super(msg, cause);
    }
}