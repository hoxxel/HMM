package main.argparser;

/**
 * exception, die bei Parse-Fehlern durch {@link ArgumentParser} ausgeloesst wird.
 *
 * @author Dominick Leppich
 * @author Soeren Metje
 */
public class ArgumentParserException extends Exception {
    /**
     * erstellt exception
     */
    public ArgumentParserException() {

    }

    /**
     * erstellt exception
     *
     * @param msg info
     */
    public ArgumentParserException(String msg) {
        super(msg);
    }

    /**
     * erstellt exception
     *
     * @param msg   info
     * @param cause grund
     */
    public ArgumentParserException(String msg, Throwable cause) {
        super(msg, cause);
    }
}