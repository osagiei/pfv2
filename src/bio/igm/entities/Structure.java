/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.entities;

import java.util.HashMap;
import java.util.Map;
import org.apache.commons.lang3.Range;

/**
 *
 * @author osagie
 */
public class Structure {

    String bed, id, coords, strand, chromosome, refseq, classification = "NA", splice_signal;
    int start, stop, read_count, exon5, exon3;
    String sequence;
    int size, As, Ts, Cs, Gs;
    double gc, at, jpm, mirna_bs_density;
    Range outer_left, outer_right, inner, ss5, ss3;
    int right_alus, left_alus, mirna_bs, ss_count;
    boolean exonic, exon_inton, intronic;

    public Structure(String _bed) {
        this.bed = _bed;
        String[] temp = _bed.split("\t");
        chromosome = temp[0];
        start = Integer.parseInt(temp[1]);
        stop = Integer.parseInt(temp[2]);
        size = stop - start;
        id = temp[3];
        coords = temp[0] + ":" + temp[1] + "-" + temp[2];
        read_count = Integer.parseInt(temp[4]);
        strand = temp[5];
        splice_signal = temp[6];

        outer_left = Range.between(start - 1000, start);
        outer_right = Range.between(stop, stop + 1000);
        inner = Range.between(start, stop);

        ss5 = Range.between(start - 1, start + 1);
        ss3 = Range.between(stop - 1, stop + 1);

    }

    public void feature_search(Map<String, String> bed_like, String feature) {
        //    if (!exonic) { TODO: Add chromosome check here
        for (String s : bed_like.values()) {
            if (chromosome.equalsIgnoreCase(s.split("\t")[0])) {
                int _start = Integer.parseInt(s.split("\t")[1]);

                if (feature.equalsIgnoreCase("alu")) {
                    if (outer_left.contains(_start)) {
                        left_alus++;
                    }
                    if (outer_right.contains(_start)) {
                        right_alus++;
                    }
                } else if (feature.equalsIgnoreCase("mirna")) {
                    if (inner.contains(_start)) {
                        mirna_bs++;
                    }
                }
            }
        }
        //    }else if(feature.equalsIgnoreCase("mirna")){

        //   }
        if (feature.equalsIgnoreCase("mirna")) {
            mirna_bs_density = (double) mirna_bs / size;
        }
    }

    public void resolve_coords_to_exons(Map<String, String> exons) {
        boolean left = false;
        boolean right = false;
        boolean annotated = false;
        exon5 = -1;
        exon3 = -1;
        Map<String, Integer> left_ids = new HashMap<>();
        Map<String, Integer> right_ids = new HashMap<>();

        for (String s : exons.values()) {
            if (chromosome.equalsIgnoreCase(s.split("\t")[0])) {

                Range _start = Range.between(Integer.parseInt(s.split("\t")[1]) - 1, Integer.parseInt(s.split("\t")[1]) + 1);
                Range _stop = Range.between(Integer.parseInt(s.split("\t")[2]) - 1, Integer.parseInt(s.split("\t")[2]) + 1);
                String _id = s.split("\t")[3];

                int exon_order = Integer.parseInt(s.split("\t")[4]);
                exon_order = s.split("\t")[5].contains("-") ? exon_order - 1 : exon_order;

                if (ss5.isOverlappedBy(_start)) {
                    //left = true;
                    exon5 = exon_order;
                    left_ids.put(_id, exon5);

                }

                if (ss3.isOverlappedBy(_stop)) {
                    //right = true;
                    exon3 = exon_order;
                    right_ids.put(_id, exon3);

                }
            }

        }
        //System.out.println(left_ids.size() + "\t" + right_ids.size());

        for (String s : left_ids.keySet()) {
            exon5 = left_ids.get(s);
            left = true;
            if (right_ids.containsKey(s)) {
                exon3 = right_ids.get(s);
                this.refseq = s;
                right = true;
                break;
            }
        }
        
        if (left && right) {
            ss_count = 2;
            classification = "exonic";

        } else if (left || right) {
            ss_count = 1;
            classification = "intronic";
        } else {
            ss_count = 0;
            classification = "intergenic_intronic";

        }
        exon5 = Math.max(exon5, exon3);
        exon3 = Math.min(exon5, exon3);

    }

    public void compute_jpm(int _count) {
        jpm = (double) read_count / _count;
        jpm *= 1000000;
        setJpm(jpm);
    }

    public String getBed() {
        return bed;
    }

    public void setBed(String bed) {
        this.bed = bed;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getCoords() {
        return coords;
    }

    public void setCoords(String coords) {
        this.coords = coords;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getRefseq() {
        return refseq;
    }

    public void setRefseq(String refseq) {
        this.refseq = refseq;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getStop() {
        return stop;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }

    public int getRead_count() {
        return read_count;
    }

    public void setRead_count(int read_count) {
        this.read_count = read_count;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
        this.size = size;
    }

    public int getAs() {
        return As;
    }

    public void setAs(int As) {
        this.As = As;
    }

    public int getTs() {
        return Ts;
    }

    public void setTs(int Ts) {
        this.Ts = Ts;
    }

    public int getCs() {
        return Cs;
    }

    public void setCs(int Cs) {
        this.Cs = Cs;
    }

    public int getGs() {
        return Gs;
    }

    public void setGs(int Gs) {
        this.Gs = Gs;
    }

    public double getGc() {
        return gc;
    }

    public void setGc(double gc) {
        this.gc = gc;
    }

    public double getAt() {
        return at;
    }

    public void setAt(double at) {
        this.at = at;
    }

    public double getJpm() {
        return jpm;
    }

    public void setJpm(double jpm) {
        this.jpm = jpm;
    }

    public double getMirna_bs_density() {
        return mirna_bs_density;
    }

    public void setMirna_bs_density(double mirna_bs_density) {
        this.mirna_bs_density = mirna_bs_density;
    }

    public int getRight_alus() {
        return right_alus;
    }

    public void setRight_alus(int right_alus) {
        this.right_alus = right_alus;
    }

    public int getLeft_alus() {
        return left_alus;
    }

    public void setLeft_alus(int left_alus) {
        this.left_alus = left_alus;
    }

    public int getMirna_bs() {
        return mirna_bs;
    }

    public void setMirna_bs(int mirna_bs) {
        this.mirna_bs = mirna_bs;
    }

    public int getSs_count() {
        return ss_count;
    }

    public void setSs_count(int ss_count) {
        this.ss_count = ss_count;
    }

    public boolean isExonic() {
        return exonic;
    }

    public void setExonic(boolean exonic) {
        this.exonic = exonic;
    }

    public boolean isExon_inton() {
        return exon_inton;
    }

    public void setExon_inton(boolean exon_inton) {
        this.exon_inton = exon_inton;
    }

    public boolean isIntronic() {
        return intronic;
    }

    public void setIntronic(boolean intronic) {
        this.intronic = intronic;
    }

    public String getClassification() {
        return classification;
    }

    public void setClassification(String classification) {
        this.classification = classification;
    }

    public int getExon5() {
        return exon5;
    }

    public void setExon5(int exon5) {
        this.exon5 = exon5;
    }

    public int getExon3() {
        return exon3;
    }

    public void setExon3(int exon3) {
        this.exon3 = exon3;
    }

    @Override
    public String toString() {
        String temp = chromosome + "\t" + start + "\t" + stop + "\t" + 
                id + "\t" + read_count + "\t" + strand + "\t" + splice_signal + "\t" +
                jpm + "\t" + ss_count + "\t" + classification + "\t" + left_alus + "\t" + right_alus
                + "\t" + mirna_bs + "\t" + mirna_bs_density;
        return temp;
    }

}
