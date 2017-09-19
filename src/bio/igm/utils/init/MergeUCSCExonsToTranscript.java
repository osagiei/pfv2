package bio.igm.utils.init;

import bio.igm.entities.Exons;
import bio.igm.entities.Transcript;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author osagie izuogu - 05/2013
 */
public class MergeUCSCExonsToTranscript {

    String path;
    Map<String, Exons> exons = new HashMap<String, Exons>();
    Map<String, Transcript> transcripts = new HashMap<String, Transcript>();
    private static Logger LOG;

    public MergeUCSCExonsToTranscript(String _path, String _wd) throws IOException {
        this.path = _path;
        File f = new File(_wd);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_wd, MergeUCSCExonsToTranscript.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), MergeUCSCExonsToTranscript.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(MergeUCSCExonsToTranscript.class.getName()).log(Level.SEVERE, null, ex);
        }
        readExons();
        buildTranscripts();
        writeToFile();
    }

    private void readExons() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        StringBuilder builder = new StringBuilder();

        String line = "";

        while ((line = br.readLine()) != null) {
            builder.append(line);
            builder.append(System.getProperty("line.separator"));
        }
        String[] seqs = builder.toString().split(">");

        for (String s : seqs) {
            if (StringUtils.isNotBlank(s)) {
                String[] oneSeq = s.trim().split("\n");
                String id = oneSeq[0];

                String seq = StringUtils.join(oneSeq, "", 1, oneSeq.length);
                Exons e = new Exons(id, seq);
                exons.put(id, e);
            }
        }
        LOG.info("Number of Exons: " + exons.size());
    }

    private void buildTranscripts() {
        for (Exons e : exons.values()) {
            if (transcripts.containsKey(e.getRefid())) {
                transcripts.get(e.getRefid()).addExon(e);
            } else {
                Transcript t = new Transcript(e.getRefid());
                t.addExon(e);
                transcripts.put(t.getRefseq(), t);
            }

        }
        LOG.info("Finished Adding Exons to Transcripts.. ");

        for (Transcript t : transcripts.values()) {

            List<Integer> temp = new ArrayList<Integer>(t.getExons().keySet());
            Collections.sort(temp);
            Integer[] index = new Integer[temp.size()];
            index = temp.toArray(index);
            String temp_seq = "";
            for (int i = 0; i < index.length; i++) {
                temp_seq += t.getExons().get(index[i]).getSequence();
            }
            t.setSequence(temp_seq);

        }
    }

    private void writeToFile() throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.path + ".mrna"));
        int counter = 0;
        for (Transcript t : transcripts.values()) {
            bw.write(">" + t.getRefseq() + "\n" + t.getSequence() + "\n");
            counter++;
        }
        bw.close();
        LOG.info("Number of mRNA Sequences: " + counter);
    }

    public static void main(String[] args) {
        String path = args[0];
        String wd = args[1];
        try {
            new MergeUCSCExonsToTranscript(path, wd);
        } catch (IOException ex) {
            Logger.getLogger(MergeUCSCExonsToTranscript.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
