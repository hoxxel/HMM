package main.argparser;

/**
 * abstrakte Parametergrundlage
 *
 *
 * @author Soeren Metje
 */
abstract class AbstractParameter<T> {
    /**
     * bezeichnung des Parameters
     */
    private final String name;
    /**
     * flag, ob der Parameter zwingend notwendig ist.
     * Loesst exception beim parsen aus, falls required == true und nicht vorhanden
     */
    private final boolean required;
    /**
     * wert des Parameters
     */
    private T value;

    /**
     * erstellt einen Parameter
     *
     * @param name     bezeichnung
     * @param required flag zwingend notwendig
     * @throws IllegalArgumentException name == null
     */
    public AbstractParameter(String name, boolean required) throws IllegalArgumentException {
        if (name == null)
            throw new IllegalArgumentException("passed name is null");
        this.name = name;
        this.required = required;
    }

    /**
     * gibt bezeichnung zurueck
     *
     * @return bezeichnung
     */
    public String getName() {
        return name;
    }

    /**
     * gibt flag zwingend notwendig zurueck
     *
     * @return flag zwingend notwendig
     */
    public boolean isRequired() {
        return required;
    }

    /**
     * gibt wert zurueck
     *
     * @return wert
     */
    public T getValue() {
        return value;
    }

    /**
     * setzt wert
     *
     * @param value zu setzender wert
     * @throws ArgumentParserException  falls aktueller wert schon gesetzt wurde
     * @throws IllegalArgumentException falls wert == null
     */
    void setValue(T value) throws ArgumentParserException, IllegalArgumentException {
        if (value == null)
            throw new IllegalArgumentException("passed value is null");
        if (this.isSet())
            throw new ArgumentParserException("Parameter " + name + " is already set");
        this.value = value;
    }

    /**
     * gibt true zurueck, falls wert gesetzt ist. Ansonsten false.
     *
     * @return true, falls wert gesetzt ist. Ansonsten false.
     */
    public boolean isSet() {
        return value != null;
    }

    /**
     * gibt true zurueck, falls die bezeichnungen der parameter uebereinstimmen. Ansonsten false.
     *
     * @param o zu vergleichendes Objekt
     * @return true, falls die bezeichnungen der parameter uebereinstimmen. Ansonsten false.
     */
    @Override
    public boolean equals(Object o) {
        if (o == null)
            return false;
        if (!(o instanceof AbstractParameter))
            return false;
        if (!name.equals(((AbstractParameter) o).getName()))
            return false;
        return true;
    }

    /**
     * gibt string-representation zurueck
     *
     * @return string-representation
     */
    @Override
    public String toString() {
        return name;
    }
}
