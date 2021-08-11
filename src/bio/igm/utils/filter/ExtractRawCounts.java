package bio.igm.utils.filter;

import bio.igm.entities.Reads;
import bio.igm.utils.init.Logging;
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
 * @author osagie
 */
public class ExtractRawCounts {

    String path;
    
    Map<String, Integer> ptes_counts = new HashMap<String, Integer>();
   

    private static Logger LOG;

    public ExtractRawCounts(String _path) throws IOException {
        this.path = _path;
    
           File f = new File(_path);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_path, ExtractRawCounts.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), ExtractRawCounts.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(ExtractRawCounts.class.getName()).log(Level.SEVERE, null, ex);
        }

        LOG.info("Generating raw counts from ptes.sam .. ");
        
        ptes_counts = extract_counts_from_sam("ptes.sam");
        LOG.info("Generating raw counts from canonical.sam .. ");
        
       // canonical_counts = extract_counts_from_sam("canonical.sam");

        writeToFile();
    }

    private void writeToFile() throws IOException {
        BufferedWriter bwP = new BufferedWriter(new FileWriter(this.path + "ptes-raw.counts"));
       

      
        for (String s : ptes_counts.keySet()) {
            bwP.write(s + "\t" + ptes_counts.get(s) + "\n");
        }
      

        bwP.close();
     
    }



    private Map<String, Integer> extract_counts_from_sam(String filename) throws IOException {
        Map<String, Integer> temp = new HashMap<String, Integer>();
        BufferedReader br = new BufferedReader(new FileReader(this.path + filename));

        String line = "";

        while ((line = br.readLine()) != null) {
            Reads r = new Reads(line);
            if (temp.containsKey(r.getTarget())) {
                int x = temp.get(r.getTarget()) + 1;
                temp.put(r.getTarget(), x);
            } else {
                temp.put(r.getTarget(), 1);
            }
        }
        return temp;

    }

    public static void main(String[] args) {
        String path = args[0];
        try {
            new ExtractRawCounts(path);
            
        } catch (IOException ex) {
            LOG.log(Level.SEVERE, null, ex);
        }
    }
}
