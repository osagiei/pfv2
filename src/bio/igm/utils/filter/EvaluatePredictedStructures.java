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
import java.util.ArrayList;
import java.util.List;
/**
 *
 * @author osagie
 */
public class EvaluatePredictedStructures {

    String path;
    
    Map<String, List<Reads>> structures = new HashMap<>();
   

    private static Logger LOG;

    public EvaluatePredictedStructures(String _path) throws IOException {
        this.path = _path;
    
           File f = new File(_path);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_path, EvaluatePredictedStructures.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), EvaluatePredictedStructures.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(EvaluatePredictedStructures.class.getName()).log(Level.SEVERE, null, ex);
        }

        LOG.info("Generating raw counts from ptes.sam .. ");

        extractStructureCounts();
        writeToFile();
    }

    private void writeToFile() throws IOException {
        BufferedWriter bwP = new BufferedWriter(new FileWriter(this.path + "/junctions.txt"));
        String header = "#chrom\tstart\tstop\tcoords\tread_count\tstrand\tdistinct_read_count\tNM_le_2\tNM_gt_2"; 

        for (String s : structures.keySet()) {
            List<Reads> reads = structures.get(s);

            String chrom = s.split(":")[0];
            String start = s.split(":")[1].split("-")[0];
            String stop = s.split(":")[1].split("-")[1];
            String strand = s.split("_")[1];

            int read_count = reads.size();

            int distinct_read_count = 0;
            int nm_gt_2 = 0;
            int nm_lt_2 = 0;
            List<Integer> starts = new ArrayList<>();

            for(Reads read : reads){
                String[] temp = read.getLine().split("\t");
                int _start = Integer.parseInt(temp[3]);
                int edit_distance = Integer.parseInt(temp[temp.length - 2].split(":")[2]);

                if(edit_distance > 2){
                    nm_gt_2 += 1;
                }else{
                    nm_lt_2 +=1;
                }

                if(!starts.contains(_start)){
                    starts.add(_start);
                }
               distinct_read_count = starts.size(); 
            }

            bwP.write(chrom + "\t" + start + "\t" + stop + "\t" + 
                s.replace("_", ":") + "\t" + read_count + "\t" + strand + "\t" + distinct_read_count +
                "\t" + nm_lt_2 + "\t" + nm_gt_2 + "\n");            
        }
      

        bwP.close();
     
    }

    private void extractStructureCounts() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "/ptes.sam"));
        String line = "";

        while (( line = br.readLine()) != null) {
            Reads read = new Reads(line);
            List<Reads> _reads = new ArrayList<>();
            if (structures.containsKey(read.getTarget())) {
                _reads = structures.get(read.getTarget());
            }
             _reads.add(read);
             structures.put(read.getTarget(), _reads);
        }

    }



    public static void main(String[] args) {
        String path = args[0];
        try {
            new EvaluatePredictedStructures(path);
            
        } catch (IOException ex) {
            LOG.log(Level.SEVERE, null, ex);
        }
    }
}
