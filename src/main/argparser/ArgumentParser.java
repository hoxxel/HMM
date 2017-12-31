package main.argparser;

/**
 * Wird zum parsen der Argumente verwendet.
 * Ist an den parser aus dem letzten Allgemeinen-Programmier-Praktikum (TowerWarsPP) angelehnt.
 * Urspruenglicher autor dabei ist Dominick Leppich, der mir die Erlaubnis erteilt hat, dies hier zu verwenden.
 *
 * @author Dominick Leppich
 * @author Soeren Metje
 */
public class ArgumentParser {

    private ParameterSet params;

    // ------------------------------------------------------------

    /**
     * Erzeuge einen neuen ArgumentParser.
     *
     * @param parameterSet Parameter
     * @throws IllegalArgumentException wenn uebergebenes {@link ParameterSet} null oder leer ist
     */
    public ArgumentParser(ParameterSet parameterSet) throws IllegalArgumentException {
        if (parameterSet == null) {
            throw new IllegalArgumentException("passed ParameterSet is null");
        }
        if (parameterSet.size() < 1)
            throw new IllegalArgumentException("passed ParameterSet is empty");

        this.params = parameterSet;
    }

    /**
     * Erzeuge einen neuen ArgumentParser und parse die Argumente
     *
     * @param parameterSet Parameter
     * @param args         Argumente
     * @throws IllegalArgumentException wenn uebergebenes {@link ParameterSet} null oder leer ist
     * @throws ArgumentParserException  wenn das Parsen der Argumente fehlschlaegt
     */
    public ArgumentParser(ParameterSet parameterSet, String[] args) throws ArgumentParserException {
        this(parameterSet);
        parseArgs(args);
    }
    // ------------------------------------------------------------

    /**
     * Parse die Argumente.
     *
     * @param args Argumente
     * @throws ArgumentParserException wenn das Parsen der Argumente fehlschlaegt
     */
    public void parseArgs(String[] args) throws ArgumentParserException {

        // Index to parse
        int index = 0;

        try {
            while (index < args.length) {
                // Check if argument is a flag or setting
                if (args[index].startsWith("--")) {
                    String parameter = args[index].substring(2);
                    Flag flag = params.getFlag(parameter);
                    if (flag == null) {
                        throw new ArgumentParserException("Flag " + parameter + " not found");
                    }

                    flag.setValue(true);
                    index += 1;
                } else if (args[index].startsWith("-")) {
                    String parameter = args[index].substring(1);
                    Setting setting = params.getSetting(parameter);
                    if (setting == null) {
                        throw new ArgumentParserException("Flag " + parameter + " not found");
                    }

                    String value = args[index + 1];
                    if (value.startsWith("-"))
                        throw new ArgumentParserException("Setting value wrong format: "
                                + value);
                    setting.setValue(value);
                    index += 2;
                } else
                    throw new ArgumentParserException("Error parsing: " + args[index]);
            }
        } catch (IndexOutOfBoundsException e) {
            throw new ArgumentParserException("Missing value");
        }

        AbstractParameter parameter = params.firstRequiredNotSetParameter();
        if (parameter != null) {
            throw new ArgumentParserException("missing parameter " + parameter);
        }
    }
}
