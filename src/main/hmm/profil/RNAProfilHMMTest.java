package main.hmm.profil;

import main.fastaparser.Sequence;
import main.hmm.profil.viterbi.Viterbi;
import main.hmm.profil.viterbi.ViterbiPath;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import static org.junit.runners.Parameterized.Parameter;
import static org.junit.runners.Parameterized.Parameters;

/**
 * Test-Klasse fuer {@link RNAProfilHMM}.
 *
 * @author Soeren Metje
 */
@RunWith(Parameterized.class)
public class RNAProfilHMMTest {

    /**
     * Trainings-Sequenzen
     */
    @Parameter(0)
    public String[] seqTrain;

    /**
     * Test-Sequenzen
     */
    @Parameter(1)
    public String seqTest;

    /**
     * Erwarteter Zustandspfad
     */
    @Parameter(2)
    public String result;

    /**
     * Liefert List mit Parametern der Testfaelle zurueck
     *
     * @return List mit Parametern der Testfaelle
     */
    @Parameters
    public static Collection<Object[]> data() {
        Object[][] data = new Object[][]{
                {new String[]{"UACAAUCAAGG", "UA-AAUCAAGG", "U--AAUCAAGG", "U-CAAUCAAGG", "U--AAUCAAGG", "U--AAUCAAGG", "UACAAUCAAGG", "U--AAUCAAGG", "U--AAUCAAGG", "UACAAUCAAGG"}, "UACAAUCAAGG", "MIIMMMMMMMM"},
                {new String[]{"UACAAUCAAGG", "UA-AAUCAAGG", "U--AAUCAAGG", "U-CAAUCAAGG", "U--AAUCAAGG", "U--AAUCAAGG", "UACAAUCAAGG", "U--AAUCAAGG", "U--AAUCAAGG", "UACAAUCAAGG"}, "", "DDDDDDDDD"}};
        return Arrays.asList(data);
    }

    /**
     * Test von {@link RNAProfilHMM}.
     * Sowohl die Erstellung des Modells, als auch der Viterbi-Algorithmus
     */
    @Test
    public void testViterbiOne() {
        ArrayList<Sequence> sequences = new ArrayList<>(seqTrain.length);
        for (int i = 0; i < seqTrain.length; i++) {
            sequences.add(new Sequence(String.valueOf(i), null, seqTrain[i]));
        }
        ProfilHMM model = new RNAProfilHMM(sequences);

        ViterbiPath path = Viterbi.viterbi(model, new Sequence("test", null, seqTest));
        String stringPath = String.valueOf(path.getStatePath());
        Assert.assertEquals(result, stringPath);
    }
}
