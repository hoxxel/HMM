package main.fastaparser;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Parser fuer das .fasta Dateiformat.
 * Parset die Datei zu einer {@link Sequence}-Liste.
 */
public class FastaParser {

    /**
     * Liesst die Datei am uebergebenen Dateipfad ein und parset sie anschliessend.
     * Liefert eine Liste mit den geparseten {@link Sequence} zurueck.
     *
     * @param filePath Dateipfad
     * @return Liste mit den geparseten {@link Sequence}
     * @throws FileNotFoundException falls Dateipfad ungueltig
     * @throws IOException           falls beim einlesen Fehler auftritt
     * @throws FastaParserException  falls der Inhalt der Datei nicht dem fasta Format entspricht
     */
    public static List<Sequence> parseFile(String filePath) throws FileNotFoundException, IOException, FastaParserException {

        File file = new File(filePath);

        List<Sequence> ret = new LinkedList<>();

        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(file))) { // closes reader
            String line;
            String description = null, comment = null;
            while ((line = bufferedReader.readLine()) != null) {
                line = line.trim();
                char firstChar = line.charAt(0);

                // Description -----------------------------------------------
                if (firstChar == '>') {
                    if (description == null)
                        description = line.substring(1);
                    else
                        throw new FastaParserException("Missing sequence!");
                }
                // Comment -----------------------------------------------
                else if (firstChar == ';') {
                    if (description == null)
                        throw new FastaParserException("Comment at wrong position or missing description!");
                    if (comment == null)
                        comment = line.substring(1); // create comment
                    else
                        comment += '\n' + line.substring(1); // add comment line
                }
                // Sequence -----------------------------------------------
                else {
                    if (description == null) {
                        throw new FastaParserException("Missing description! (line starting with >)");
                    }
                    ret.add(new Sequence(description, comment, line));
                    description = null;
                    comment = null;
                }
            }
        }

        ret = new ArrayList<>(ret);

        return ret;
    }
}
