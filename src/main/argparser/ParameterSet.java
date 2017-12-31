package main.argparser;

import java.util.ArrayList;
import java.util.List;

/**
 * Set von Parametern
 *
 *
 * @author Soeren Metje
 */
public class ParameterSet {

    /**
     * enthaelt hinzugefuegte Parameter
     */
    private List<AbstractParameter> params;

    /**
     * erstellt einen neues {@link ParameterSet}
     */
    public ParameterSet() {
        params = new ArrayList<>();
    }

    /**
     * fuegt ein {@link Setting} hinzu und gibt referenz auf {@link ParameterSet} zurueck um Konkatenation zu ermoeglichen
     *
     * @param setting Setting
     * @return referenz auf {@link ParameterSet}
     */
    public ParameterSet addSetting(Setting setting) {
        return addParameter(setting);
    }

    /**
     * fuegt ein {@link Flag} hinzu und gibt referenz auf {@link ParameterSet} zurueck um Konkatenation zu ermoeglichen
     *
     * @param flag Flag
     * @return referenz auf {@link ParameterSet}
     */
    public ParameterSet addFlag(Flag flag) {
        return addParameter(flag);
    }

    /**
     * fuegt ein {@link AbstractParameter} hinzu und gibt referenz auf {@link ParameterSet} zurueck um Konkatenation zu ermoeglichen
     *
     * @param parameter AbstractParameter
     * @return referenz auf {@link ParameterSet}
     */
    private ParameterSet addParameter(AbstractParameter parameter) {
        if (parameter == null)
            throw new IllegalArgumentException("passed parameter == null");
        if (params.contains(parameter))
            throw new IllegalArgumentException("passed parameter is already added");

        params.add(parameter);
        return this;
    }

    /**
     * gibt Falg entsprechend dem uebergebenem Parameter oder null zurueck
     *
     * @param parameter Parameter
     * @return Falg oder null
     */
    Flag getFlag(String parameter) {
        try {
            return (Flag) getParameter(parameter);
        } catch (ClassCastException e) {
            throw new IllegalArgumentException("parameter matching to passed String is not a Flag");
        }
    }

    /**
     * gibt Setting entsprechend dem uebergebenem Parameter oder null zurueck
     *
     * @param parameter Parameter
     * @return Setting oder null
     */
    Setting getSetting(String parameter) {
        try {
            return (Setting) getParameter(parameter);
        } catch (ClassCastException e) {
            throw new IllegalArgumentException("parameter matching to passed String is not a Setting");
        }
    }

    /**
     * gibt {@link AbstractParameter} entsprechend dem uebergebenem Parameter oder null zurueck
     *
     * @param parameter Parameter
     * @return AbstractParameter oder null
     */
    private AbstractParameter getParameter(String parameter) {
        if (parameter == null)
            throw new IllegalArgumentException("passed String is null");
        for (AbstractParameter flag :
                params) {
            if (parameter.equals(flag.getName()))
                return flag;
        }
        return null;
    }

    /**
     * gibt das erste {@link AbstractParameter} zurueck, das required und nicht gesetzt ist.
     *
     * @return {@link AbstractParameter}, das required und nicht gesetzt ist.
     */
    AbstractParameter firstRequiredNotSetParameter() {
        for (AbstractParameter parameter :
                params) {
            if (parameter.isRequired() && !parameter.isSet())
                return parameter;
        }
        return null;
    }

    /**
     * gibt anzahl der Parameter zurueck
     *
     * @return anzahl der Parameter
     */
    int size() {
        return params.size();
    }

    /**
     * gibt true zurueck, falls das uebergebene {@link Setting} schon vorhanden ist. Ansonsten false.
     *
     * @param setting zu ueberpruefendes {@link Setting}
     * @return true, falls das uebergebene {@link Setting} schon vorhanden ist. Ansonsten false.
     */
    public boolean contains(Setting setting) {
        return contains((AbstractParameter) setting);
    }


    /**
     * gibt true zurueck, falls das uebergebene {@link Flag} schon vorhanden ist. Ansonsten false.
     *
     * @param flag zu ueberpruefendes {@link Flag}
     * @return true zurueck, falls das uebergebene {@link Flag} schon vorhanden ist. Ansonsten false.
     */
    public boolean contains(Flag flag) {
        return contains((AbstractParameter) flag);
    }

    /**
     * gibt true zurueck, falls der uebergebene {@link AbstractParameter} schon vorhanden ist. Ansonsten false.
     *
     * @param abstractParameter zu ueberpruefende {@link AbstractParameter}
     * @return true, falls der uebergebene {@link AbstractParameter} schon vorhanden ist. Ansonsten false.
     */
    private boolean contains(AbstractParameter abstractParameter) {
        return params.contains(abstractParameter);
    }


    /**
     * gibt true zurueck, falls das {@link Setting} entsprechend dem uebergebenem parameter schon vorhanden ist. Ansonsten false.
     *
     * @param parameter zu ueberpruefender parameter
     * @return true, falls das {@link Setting} entsprechend dem uebergebenem parameter schon vorhanden ist. Ansonsten false.
     */
    public boolean contains(String parameter) {
        for (AbstractParameter mParameter :
                params) {
            if (mParameter.getName().equals(parameter)) {
                return true;
            }
        }
        return false;
    }
}
