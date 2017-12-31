package main.argparser;

/**
 * Parameter, der einen String als Wert besitzt.
 *
 *
 * @author Soeren Metje
 */
public class Setting extends AbstractParameter<String> {

    /**
     * erstellt Setting
     *
     * @param name     bezeichnung
     * @param required flag zwingend notwendig
     * @throws IllegalArgumentException name == null
     */
    public Setting(String name, boolean required) {
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
