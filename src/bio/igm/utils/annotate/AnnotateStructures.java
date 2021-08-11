package bio.igm.utils.annotate;

import bio.igm.entities.Structure;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.Range;

/**
 *
 * @author osagie
 */
public class AnnotateStructures {
    /*
     * read files
     * resolve coordinates to exons / classify
     * 
     */

    String ptes_bed, exons_bed, alus_bed, mirna_bed, transcripts_bed, name, path;
    Map<String, String> alus = new HashMap<String, String>();
    Map<String, String> mirnas = new HashMap<>();
    Map<String, String> exons = new HashMap<>();
    //Map<String, String> exons_w_ids = new HashMap<>();
    Map<String, Structure> structures = new HashMap<>();
    Structure structure;
    int junctions;

    public AnnotateStructures(String _ptes_bed, String _exons_bed, String _alus_bed, String _mirna_bed,
            String _transcripts_bed, String _name, String _path) throws IOException {
        this.ptes_bed = _ptes_bed;
        this.exons_bed = _exons_bed;
        this.alus_bed = _alus_bed;
        this.mirna_bed = _mirna_bed;

        this.name = _name;
        this.path = _path;

        get_canonical_junction_count();
        get_structures();
        get_exons();
        get_alus();
        get_mirna();

        writeToFile();

    }

    public AnnotateStructures(String _ptes_bed, String _exons_bed, String _transcripts_bed, String _name, String _path) throws IOException {
        this.ptes_bed = _ptes_bed;
        this.exons_bed = _exons_bed;
        this.transcripts_bed = _transcripts_bed;

        this.name = _name;
        this.path = _path;

        get_canonical_junction_count();
        get_structures();
        get_exons();

        get_alus();
        get_mirna();
        annotate_non_exonic_structures();

        writeToFile();

    }

    private void get_canonical_junction_count() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "/flanking-canonical-counts.tsv.bed"));

        String line = "";

        while ((line = br.readLine()) != null) {
            if (!line.split("\t")[0].contains("M")) {
                junctions += Integer.parseInt(line.split("\t")[4]);
            }
        }

        br.close();

    }

    private void get_structures() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.ptes_bed));
        String line = "";

        while ((line = br.readLine()) != null) {
            junctions += Integer.parseInt(line.split("\t")[4]);
            Structure s = new Structure(line);
            if(structures.containsKey(s.getCoords())){
                structure = structures.get(s.getCoords());
                if(structure.getStrand().equalsIgnoreCase(s.getStrand())){
                    int count = s.getRead_count() + structure.getRead_count();
                    structure.setRead_count(count);
                    structures.put(structure.getCoords(), structure);
                }
            }else{
                structures.put(s.getCoords(), s);
            }
            
        }
        br.close();

    }

    private void get_exons() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.exons_bed));
        String line = "";

        while ((line = br.readLine()) != null) {
            String id = line.split("\t")[3] + "_" + line.split("\t")[4];
            exons.put(id, line);
        }
        br.close();
        int counter = 0;
        for (Structure s : structures.values()) {
            structure = s;
            structure.resolve_coords_to_exons(exons);

            structures.put(structure.getCoords(), structure);
            counter++;
        }

        System.out.println("Finished reading exons and classifying structures..");
    }

    private void get_alus() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "/alus.bed"));
        String line = "";

        while ((line = br.readLine()) != null) {
            String id = line.split("\t")[3].split("_")[0];
            if (structures.containsKey(id)) {
                structure = structures.get(id);
                String[] temp = line.split("\t");
                structure.setLeft_alus(Integer.parseInt(temp[temp.length - 2]));
                structure.setRight_alus(Integer.parseInt(temp[temp.length - 1]));

                structures.put(id, structure);
            }
        }
        br.close();

        System.out.println("Finished annotating structures with alus..");
    }

    private void get_mirna() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(this.path + "/mirnas.bed"));
        String line = "";

        while ((line = br.readLine()) != null) {
            String id = line.split("\t")[0] + ":" + line.split("\t")[1] + "-" + line.split("\t")[2];
            //mirnas.put(id, line);

            if (structures.containsKey(id)) {
                structure = structures.get(id);
                if (structure.getSs_count() == 2) {
                    try {
                        structure = get_exonic_circrna_mirnabs_count(structure);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                } else {
                    String[] temp = line.split("\t");
                    structure.setMirna_bs(Integer.parseInt(temp[temp.length - 1]));
                    structure.setMirna_bs_density((double) structure.getMirna_bs() / structure.getSize());
                }
                structure.compute_jpm(junctions);

                structures.put(id, structure);
            }
        }
        br.close();

        System.out.println("Finished annotating structures with mirna_bs..");
    }

    public Structure get_exonic_circrna_mirnabs_count(Structure structure) {
        int prime5 = structure.getExon5();
        int prime3 = structure.getExon3();
        int size = 0;
        int mirna_bs = 0;
        boolean strand = structure.getStrand().contains("+") ? true : false;
        if (!strand) {
            prime3 += 1;
            prime5 += 1;
        }
        String refseq = structure.getRefseq();
        if (refseq.contains(":")) {
            refseq = refseq.replace(":", ".");
        }
        for (int i = prime3; i <= prime5; i++) {

            String temp = refseq + "_" + i;
            String exon = exons.get(temp);
            size += (Integer.parseInt(exon.split("\t")[2]) - Integer.parseInt(exon.split("\t")[1]));
            mirna_bs += Integer.parseInt(exon.split("\t")[6]);
        }

        structure.setMirna_bs(mirna_bs);
        structure.setSize(size);
        structure.setMirna_bs_density((double) mirna_bs / size);

        return structure;

    }

    private void annotate_non_exonic_structures() throws IOException {
        Map<String, Range> transcripts = new HashMap<>();
        BufferedReader br = new BufferedReader(new FileReader(transcripts_bed));
        String line = "";
        while ((line = br.readLine()) != null) {
            String[] temp = line.split("\t");
            String id = temp[0] + ":" + temp[1] + "-" + temp[2] + ":" + temp[5];
            Range r = Range.between(Integer.parseInt(temp[1]), Integer.parseInt(temp[2]));
            transcripts.put(id, r);
        }
        br.close();
        //annotate_structures
        for (Structure structure : structures.values()) {
            if (structure.getClassification().contains("intronic")) {
                String classification = "intronic";
                boolean found = false;
                for (String s : transcripts.keySet()) {
                    String strand = s.split(":")[2];
                    Range r = transcripts.get(s);
                    if (r.contains(structure.getStart())) {
                        found = true;
                        if (!strand.equalsIgnoreCase(structure.getStrand())) {
                            classification += "_antisense";
                        }

                        break;
                    }
                }
                if (found) {
                    structure.setClassification(classification);
                } else {
                    structure.setClassification("intergenic");
                }
                structures.put(structure.getCoords(), structure);
            }

        }
    }

    private void writeToFile() throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.path + "/" + name + ".bed"));

        String header = "chr\tstart\tstop\tid\tcount\tstrand\tsplice_motif\tjpm\tname\tsplice_site_count\t"
                + "region\tleft_alus\tright_alus\tmirna_bs\tmirna_bs_density";

        bw.write(header + "\n");

        for (Structure s : structures.values()) {
            bw.write(s.toString() + "\n");
        }
        bw.close();
    }

    public static void main(String[] args) {
        try {
            //new AnnotateStructures(args[0], args[1], args[2], args[3], args[4], args[5]);
            new AnnotateStructures(args[0], args[1], args[2], args[3], args[4]);
        } catch (IOException ex) {
            Logger.getLogger(AnnotateStructures.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
