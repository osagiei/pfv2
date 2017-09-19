/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.init;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author osagie izuogu - 05/2013 Class converts transcriptome bed obtained from UCSC
 * to coordinates file required by Pipeline 
 * input: path to transcriptome bed 
 * output: coordinates written to path.coords
 */
public class ConvertBedToCoordinates {

    String path;
    Map<String, String> genes = new HashMap<String, String>();
    Map<String, String> coords = new HashMap<String, String>();
    private static Logger LOG;

    public ConvertBedToCoordinates(String _path, String _wd) throws IOException {
        this.path = _path;
        File f = new File(_wd);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_wd, ConvertBedToCoordinates.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), ConvertBedToCoordinates.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(ConvertBedToCoordinates.class.getName()).log(Level.SEVERE, null, ex);
        }



        read_input_bed();
        convert_to_coords();

        writeToFile();

    }

    private void read_input_bed() throws IOException {
        LOG.entering(getClass().getName(), "read_input_bed() - Preparing Coordinates File...");
        BufferedReader br = new BufferedReader(new FileReader(this.path));
        String line = "";

        while ((line = br.readLine()) != null) {
            String id = line.split("\t")[3];    //extracts the refseq id on bed entry

            if ((id.startsWith("NM")) || (id.startsWith("NR"))) {
                genes.put(id, line);
            } else {
                String _id = id.split("_")[0].replace(".", ":"); //assumes knowngene id
                String _line = "";
                String[] content = line.split("\t");
                content[3] = _id;
                for(int i = 0; i < content.length; i++){
                    _line += content[i] + "\t";
                }
                
                
                genes.put(_id, _line.trim()); 
                
            }
            
        }
    }

    private void convert_to_coords() {
        LOG.entering(getClass().getName(), "convert_to_coords() - Converting input BED to coords file");
        for (String s : genes.keySet()) {
            String new_coord = "";
            String[] sizes = null;
            if (genes.get(s).split("\t")[5].equalsIgnoreCase("-")) {
                sizes = negative_strand_reverse(genes.get(s).split("\t")[10].split(","));
            } else {
                sizes = genes.get(s).split("\t")[10].split(",");
            }

            String starts = "";
            String ends = "";
            int x = 0;
            int y = 0;
            new_coord = s + "\t" + s;
            for (int i = 0; i < sizes.length; i++) {
                starts += x + 1 + " ";
                y = x + Integer.parseInt(sizes[i]);
                ends += y + " ";
                x = y;
            }
            new_coord += "\t" + starts + "\t" + ends;
            coords.put(s, new_coord);

        }
    }

    private String[] negative_strand_reverse(String[] data) {

        int left = 0;
        int right = data.length - 1;

        while (left < right) {
            // swap the values at the left and right indices
            String temp = data[left];
            data[left] = data[right];
            data[right] = temp;

            // move the left and right index pointers in toward the center
            left++;
            right--;
        }
        return data;
    }

    private void writeToFile() throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.path + ".coords"));

        for (String s : coords.values()) {
            bw.write(s + "\n");
        }
        bw.close();
        LOG.info("Finished processing coordinates File...");
    }

    public static void main(String[] args) {
        String path = args[0];
        String wd = args[1];
        try {
            new ConvertBedToCoordinates(path, wd);
        } catch (IOException ex) {
            LOG.log(Level.SEVERE, null, ex);
        }
    }
}
