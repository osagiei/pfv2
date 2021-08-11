package bio.igm.utils.annotate;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author osagie
 */
public class GenerateFlankingSequenceBeds {
    String ptes_bed, path, id;
    Map<String, String> left_flank, left_internal, right_flank, right_internal, annotated;
    int window_size;
    
    public GenerateFlankingSequenceBeds(String _ptes_bed, String _path, String _id, int _window_size) throws IOException{
        this.path = _path;
        this.id = _id;
        this.ptes_bed = _ptes_bed;
        this.window_size = _window_size;
        
        left_flank = new HashMap<>();
        left_internal = new HashMap<>();
        right_flank = new HashMap<>();
        right_internal = new HashMap<>();
        annotated = new HashMap<>();
        
        getPTESStructures();
        extractSegments();
        writeToFile();
        
    }

    private void getPTESStructures() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.ptes_bed));
        String line = "";
        
        br.readLine(); //removes header
        
        while((line = br.readLine()) != null){
            annotated.put(line.split("\t")[3], line);
        }
        br.close();
    }

    private void extractSegments() {
        for(String s : annotated.keySet()){
            String [] structure = annotated.get(s).split("\t");
            String lf = structure[0] + " \t"+ (Integer.parseInt(structure[1]) - window_size) + "\t" +
                    structure[1] + "\t" + s + "\t.\t" + structure[5];
            String li = structure[0] + "\t" + structure[1] + " \t"+ (Integer.parseInt(structure[1]) + window_size) + "\t" +
                     s + "\t.\t" + structure[5];
            String rf = structure[0] + "\t" + structure[2] + " \t"+ (Integer.parseInt(structure[2]) + window_size) + "\t" +
                     s + "\t.\t" + structure[5];
            String ri = structure[0] + " \t"+ (Integer.parseInt(structure[2]) - window_size) + "\t" +
                    structure[2] + "\t" + s + "\t.\t" + structure[5];
            
            left_flank.put(s, lf);
            left_internal.put(s, li);
            right_flank.put(s, rf);
            right_internal.put(s, ri);
        }
    }

    private void writeToFile() throws IOException {
        BufferedWriter blf = new BufferedWriter(new FileWriter(this.path + "/" + this.id + "_left_flanking.bed"));
        BufferedWriter bli = new BufferedWriter(new FileWriter(this.path + "/" + this.id + "_left_internal.bed"));
        BufferedWriter brf = new BufferedWriter(new FileWriter(this.path + "/" + this.id + "_right_flanking.bed"));
        BufferedWriter bri = new BufferedWriter(new FileWriter(this.path + "/" + this.id + "_right_internal.bed"));
        
        for(String s : annotated.keySet()){
            blf.write(left_flank.get(s) + "\n");
            bli.write(left_internal.get(s) + "\n");
            brf.write(right_flank.get(s) + "\n");
            bri.write(right_internal.get(s) + "\n");
        }
        
        blf.close();
        bli.close();
        brf.close();
        bri.close();
    }
    
    public static void main(String[] args){
        try {
            new GenerateFlankingSequenceBeds(args[0], args[1], args[2], Integer.parseInt(args[3]));
        } catch (IOException ex) {
            System.out.println("Class requires as input:\nPTESFinder.bed, working_directory, id and "
                    + "window_size");
            Logger.getLogger(GenerateFlankingSequenceBeds.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
}
