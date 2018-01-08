package main.hmm.profil;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class FastaParser {

    public static List<String> parseFile(String filePath) throws FileNotFoundException, IOException {

        File file = new File(filePath);

        List<String> ret = new LinkedList<>();

        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                line = line.trim();
                if (line.charAt(0) != '>')
                    ret.add(line);
            }
        }

        ret = new ArrayList<>(ret);

        return ret;
    }
}
