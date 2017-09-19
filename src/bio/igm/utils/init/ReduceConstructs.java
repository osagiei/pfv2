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
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author osagie
 */
public class ReduceConstructs {

    String path;
    Map<String, String> constructs = new HashMap<String, String>();
    Map<String, String> can_constructs = new HashMap<String, String>();
    int segment_size = 65;
    private static Logger LOG;

    public ReduceConstructs(String _path, int _segment_size) throws IOException {
        this.path = _path;
        File f = new File(_path);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_path, ReduceConstructs.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), ReduceConstructs.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(ReduceConstructs.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.segment_size = _segment_size;
        LOG.info("Reading supplied reference FASTA files..");
        constructs = readConstructs("tempConstructs.fa");
        can_constructs = readConstructs("tempCan.fa");

        LOG.info("Shrinking constructs to specified segment length..");
        reduce_constructs("ptes");
        reduce_constructs("canonical");
    }

    public ReduceConstructs(String _path, String _ptes_path, String _can_path, int _segment_size) throws IOException {
        this.path = _path;
        File f = new File(_path);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_path, ReduceConstructs.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), ReduceConstructs.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(ReduceConstructs.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.segment_size = _segment_size;
        LOG.info("Reading supplied reference FASTA files..");
        File p = new File(_ptes_path);
        if (p.exists()) {
            constructs = readConstructs(_ptes_path);
        } else {
            LOG.info("No PTES FASTA reference supplied - running unguided analysis..");
        }
        File c = new File(_can_path);
        if (c.exists()) {
            can_constructs = readConstructs(_can_path);
        } else {
            LOG.info("No Canonical Junctions FASTA reference supplied");
        }


        LOG.info("Shrinking constructs to specified segment length..");
        reduce_constructs("ptes");
        reduce_constructs("canonical");
    }

    private Map<String, String> readConstructs(String filename) throws IOException {
        Map<String, String> temp = new HashMap<String, String>();

        BufferedReader br = new BufferedReader(new FileReader(filename));
        StringBuilder builder = new StringBuilder();

        String line = "";


        while ((line = br.readLine()) != null) {
            builder.append(line);
            builder.append(System.getProperty("line.separator"));

        }

        String[] seqs = builder.toString().split(">");
        String seq = "";

        for (String s : seqs) {
            if (StringUtils.isNotBlank(s)) {
                String[] oneSeq = s.trim().split("\n");
                String id = oneSeq[0];

                seq = StringUtils.join(oneSeq, "", 1, oneSeq.length);
                temp.put(id, seq);
            }
        }

        br.close();
        return temp;
    }

    private void reduce_constructs(String structure) throws IOException {
        String filename = "";
        Map<String, String> temp_in = new HashMap<String, String>();
        Map<String, String> temp_out = new HashMap<String, String>();

        if (structure.equalsIgnoreCase("ptes")) {
            filename = "tempConstructs.fasta";
            temp_in = constructs;
        } else {
            filename = "tempCan.fasta";
            temp_in = can_constructs;
        }
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.path + filename));
        for (String s : temp_in.keySet()) {
            int junction = Integer.parseInt(s.split("\\.")[3]);
            String[] seq = trim_seq(temp_in.get(s), junction);
            String id = s.split("\\.")[0] + "." + s.split("\\.")[1] + "." + s.split("\\.")[2] + "." + seq[1];
            bw.write(">" + id + "\n" + seq[0] + "\n");

        }
        bw.close();
    }

    private String[] trim_seq(String sequence, int junction) {
        String[] temp = new String[2];


        if (sequence.length() > junction) {
            String seq1 = sequence.substring(0, junction);
            String seq2 = sequence.substring(junction);
            int x1 = seq1.length() > segment_size ? segment_size : seq1.length();
            int x2 = seq2.length() > segment_size ? segment_size : seq2.length();
            String temp1 = seq1.substring(seq1.length() - x1);
            String temp2 = seq2.substring(0, x2);

            temp[0] = temp1 + temp2;
            temp[1] = x1 + "";
            return temp;
        }

        return new String[]{"", ""};
    }

    public static void main(String[] args) {
        String path = args[0];
        String ppath = args[2];
        String cpath = args[3];
        int segementSize = 65;
        try {

            segementSize = Integer.parseInt(args[1]);
        } catch (NumberFormatException nfe) {
            LOG.info("Error parsing input parameter, proceeding with default value ..");
        }
        try {

            new ReduceConstructs(path, ppath, cpath, segementSize);

        } catch (IOException ex) {
            LOG.log(Level.SEVERE, null, ex);
        }
    }
}
