package bio.igm.utils.filter;

import bio.igm.entities.Reads;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author osagie izuogu - 05/2013
 */
public class FilterGenomicHits {

    Map<String, Reads> reads = new HashMap<String, Reads>();
    String path;

    /**
     *
     * @param path
     * @throws IOException
     */
    public FilterGenomicHits(String _path) throws IOException {

        this.path = _path;

        read_ptes_sam();
        PipelineFilter.LOG.info("Comparing reads aligned to constructs to reads aligned to genomic reference..");
        //filter_genomic_reads();
    }
    /*
     * Read reads aligning to PTES constructs
     */

    private void read_ptes_sam() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "ptes.sam"));
        String line = "";
        int counter = 0;

        while ((line = br.readLine()) != null) {

            try {
                if ((line.split("\t").length > 13)) {
                    Reads r = new Reads(line);
                    counter++;
                   // if (Integer.parseInt(r.getEditDistance().split(":")[2]) < 10) { //allow maximum of 5 edits before filtering
                        this.reads.put(line.split("\t")[0], r);

                        if ((counter % 5000000) == 0) {
                            filter_genomic_reads();
                            reads = new HashMap<String, Reads>();
                            counter = 0;
                        }
                    //}
                }
            } catch (Exception e) {
            }
        }
        if (counter < 5000000) {
            filter_genomic_reads();
            reads = new HashMap<String, Reads>();
        }
        br.close();

    }
    /*
     * Read genomic.sam and compare alignment qualities to reads aligned to
     * ptes structures.
     */

    private void filter_genomic_reads() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "genomic.sam"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.path + "genomic-filtered-out.sam", true));
        BufferedWriter bw2 = new BufferedWriter(new FileWriter(this.path + "genomic-better.sam", true));
        BufferedWriter bp = new BufferedWriter(new FileWriter(this.path + "unique.sam", true));
        String line = "";
        String id = "";
        Reads r = null;

        while ((line = br.readLine()) != null) {
            if ((line.split("\t").length > 13)) {
                id = line.split("\t")[0];
                if (reads.containsKey(id)) {
                    r = reads.get(id);
                    if (!compareNMs(r.getLine(), line)) { //if ptes read alignment is worse
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

    private boolean compareNMs(String ptes_read, String genomic_read) {
        boolean better = false;

        int x = extract_aligned(ptes_read.split("\t")[ptes_read.split("\t").length - 2]); //number of aligned bp in ptes read
        int y = extract_aligned(genomic_read.split("\t")[genomic_read.split("\t").length - 2]); //number of aligned bp to genomic target
        int ptes_nm = Integer.parseInt(ptes_read.split("\t")[ptes_read.split("\t").length - 3].split(":")[2]); //ptes read edit distance
        int genomic_nm = Integer.parseInt(genomic_read.split("\t")[genomic_read.split("\t").length - 3].split(":")[2]); //genomic read edit distance

        //if extracted MD is higher in ptes read and NM is lower - better alignment to contruct
        if ((x >= y) && (ptes_nm < genomic_nm)) {
            better = true;
        }
        return better;
    }

    private int extract_aligned(String MD) {

        String pattern = "[0-9]+";      // checks for all numbers in MD
        int temp = 0;
        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(MD.split(":")[2]);

        while (m.find()) {
            temp += Integer.parseInt(m.group()); // calculates number of aligned nucleotides - from MD field
        }
        return temp;
    }

    /**
     *
     * @return
     */
    public Map<String, Reads> getReads() {
        return reads;
    }

    /**
     *
     * @param reads
     */
    public void setReads(Map<String, Reads> reads) {
        this.reads = reads;
    }

    /**
     *
     * @return
     */
    public String getPath() {
        return path;
    }

    /**
     *
     * @param path
     */
    public void setPath(String path) {
        this.path = path;
    }
}
