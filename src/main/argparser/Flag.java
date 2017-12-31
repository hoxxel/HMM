package main.argparser;

/**
 * Parameter, der den wert true und false annehmen kann.
 *
 *
 * @author Soeren Metje
 */
public class Flag extends AbstractParameter<Boolean> {

    /**
     * erstellt Flag
     *
     * @param name     bezeichnung
     * @param required flag zwingend notwendig
     * @throws IllegalArgumentException name == null
     */
    public Flag(String name, boolean required) throws IllegalArgumentException {
        super(name, required);
    }

    /**
     * gibt string-representation zurueck
     *
     * @return string-representation
     */
    @Override
    public String toString() {
        return super.toString() + "=" + getValue();
    }
}
