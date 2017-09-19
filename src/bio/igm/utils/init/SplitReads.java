package bio.igm.utils.init;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class SplitReads {

    String path;
    String out;
    static StringBuilder left = new StringBuilder();
    static StringBuilder right = new StringBuilder();

    public SplitReads(String path, String out) throws FileNotFoundException, IOException {
        this.out = out;
        this.path = path;
        readFastQ(path);

    }

    private void readFastQ(String s) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(s));
        String p = new StringBuilder().append(this.out).append("left.fastq").toString();
        String q = new StringBuilder().append(this.out).append("right.fastq").toString();
        BufferedWriter bwl = new BufferedWriter(new FileWriter(p));
        BufferedWriter bwr = new BufferedWriter(new FileWriter(q));
        String line = "";
        String l = "";
        String r = "";

        int counter = 4;
        int index = 0;
        while ((line = br.readLine()) != null) {
            index = counter % 4;
            switch (index) {
                case 0:
                    l = line;
                    r = line;
                    break;
                case 1:
                    l = line.substring(0, 20);
                    r = line.substring(line.length() - 20, line.length());
                    break;
                case 2:
                    l = line;
                    r = line;
                    break;
                case 3:
                    l = line.substring(0, 20);
                    r = line.substring(line.length() - 20, line.length());
            }
            bwl.write(l + "\n");
            bwr.write(r + "\n");
            

            counter++;
        }
        bwl.close();
        bwr.close();
       
    }

   

    public static void main(String[] args) {
        try {
            String path = args[0];
            String out = args[1];
            new SplitReads(path, out);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SplitReads.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SplitReads.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}