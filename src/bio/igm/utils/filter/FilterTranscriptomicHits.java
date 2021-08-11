/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.filter;

import bio.igm.entities.Reads;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author osagie izuogu - 05/2013
 */
public class FilterTranscriptomicHits {

    private String path;
    private Map<String, Reads> reads;
    

    FilterTranscriptomicHits(Map<String, Reads> reads, String _path) throws IOException {
        this.path = _path;
        
        this.reads = reads;
        PipelineFilter.LOG.info("Comparing reads aligned to constructs to reads aligned to transcriptomic reference..");
        filter_refseq_reads();
    }

    private void filter_refseq_reads() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "transcriptomic.sam"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.path + "transcriptomic-filtered-out.sam"));
        BufferedWriter bw2 = new BufferedWriter(new FileWriter(this.path + "transcriptomic-better.sam"));
        BufferedWriter bp = new BufferedWriter(new FileWriter(this.path + "transcriptomic-unique.sam"));
        String line = "";
        String id = "";
        Reads r = null;

        while ((line = br.readLine()) != null) {
            if ((line.split("\t").length > 13)) {
                id = line.split("\t")[0];
                if (reads.containsKey(id)) {
                    r = reads.get(id);
                    if (!compareNMs(r.getLine(), line)) {
                        bw.write(r.getLine() + "\n");
                        bw2.write(line + "\n");
                        reads.remove(r.getId());
                    } 
                }
            }
        }
        for (Reads p : reads.values()) {

            bp.write(p.getLine() + "\n");
        }

        br.close();
        bw.close();
        bw2.close();
        bp.close();
    }

    private boolean compareNMs(String ptes_read, String refseq_read) {
        boolean better = false;

        int x = extract_aligned(ptes_read.split("\t")[ptes_read.split("\t").length - 2]);
        int y = extract_aligned(refseq_read.split("\t")[refseq_read.split("\t").length - 2]);
        int ptes_nm = Integer.parseInt(ptes_read.split("\t")[ptes_read.split("\t").length - 3].split(":")[2]);
        int refseq_nm = Integer.parseInt(refseq_read.split("\t")[refseq_read.split("\t").length - 3].split(":")[2]);

        if ((x >= y) && (ptes_nm < refseq_nm)) {
            better = true;
        }
        return better;
    }

    private int extract_aligned(String MD) {
        String pattern = "[0-9]+";
        int temp = 0;
        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(MD.split(":")[2]);

        while (m.find()) {
            temp += Integer.parseInt(m.group());
        }
        return temp;
    }

    public Map<String, Reads> getReads() {
        return reads;
    }

    public void setReads(Map<String, Reads> reads) {
        this.reads = reads;
    }

    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }
}
