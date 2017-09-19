/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.annotate;


import bio.igm.entities.Structure;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author osagie
 */
@Deprecated
public class AnnotateCircRNAs {

    static String _ptes_bed;
    static String _exons_bed;
    static String _alus_bed;
    static String _mirna_bed;
    static String _name;
    static String _exons_fa;
    static int kmer;
    static String _path;

    public static void mirnabs_search(String ptespath, String exonspath, String path, String id) throws IOException {
        Map<String, String> exons = get_exons(exonspath);
        Map<String, String> structures = get_bed_like(ptespath);
        //Map<String, String> mirna = get_bed_like(mirnapath);

        for (String s : structures.keySet()) {
            structures.put(s, get_mirnabs_count(structures.get(s), exons));
        }

        path += "/" + id + "_mirnabs.bed";
        writeToFile(structures, path);

        System.out.println("[ Finished searching for miRNA binding sites @ :"
                + new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime()) + " ]\n");
    }

    
    private static Map<String, String> get_bed_like(String path) throws IOException {
        Map<String, String> bed_like = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(path));

        String line = "";

        while ((line = br.readLine()) != null) {
            String id = line.split("\t")[0] + ":" + line.split("\t")[1] + "-" + line.split("\t")[2];
            bed_like.put(id.trim(), line);
        }
        br.close();

        return bed_like;
    }

    public static String get_mirnabs_count(String s, Map<String, String> exons) {

        int prime5 = Integer.parseInt(s.split("\t")[3].split("\\.")[1]);
        int prime3 = Integer.parseInt(s.split("\t")[3].split("\\.")[2]);
        int mirna_bs = 0;
        boolean strand = s.split("\t")[5].contains("+") ? true : false;
        if (!strand) {
            prime3 += 1;
            prime5 += 1;
        }
        String refseq = s.split("\t")[3].split("\\.")[0];
        if (refseq.contains(":")) {
            refseq = refseq.replace(":", ".");
        }
        for (int i = prime3; i <= prime5; i++) {

            String temp = refseq + "_" + i;
            String exon = exons.get(temp);
            mirna_bs += Integer.parseInt(exon.split("\t")[6]);
        }

        s += "\t" + mirna_bs;


        return s;

    }

    public static void merge_annotations(String path, String id) throws IOException {
        Map<String, String> structures = new HashMap<>();

        String apath = path + "/" + id + "_alus.bed";
        Map<String, String> alus = get_exons(apath);


        String mpath = path + "/" + id + "_mirnabs.bed";
        Map<String, String> mirnas = get_bed_like(mpath);


        for (String s : mirnas.keySet()) {
             String structure = mirnas.get(s);
           

            if (alus.containsKey(structure.split("\t")[3])) {
                String[] temp = alus.get(structure.split("\t")[3]).split("\t");
                structure += "\t" + temp[temp.length - 2] + "\t" + temp[temp.length - 1];
            } else {
                structure += "\tNA\tNA";
            }

            structures.put(s, structure);
        }
        path += "/" + id + "_annotated.bed";
        writeToFile(structures, path);

        System.out.println("[ Finished merging annotations  @ :"
                + new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime()) + " ]\n");

    }

    public static void post_annotate(String ptesbed, String ptessam, String path, String id) throws IOException {
        Map<String, String> structures = get_bed_like(ptesbed);
        BufferedWriter bw = new BufferedWriter(new FileWriter(path + "/" + id + ".bed"));

        //structures = extract_aligned_reads(ptessam, structures);
       

        String header = "chr\tstart\tstop\tid\t.\tstrand\tsplice_junction\tleft_alus\tright_alus";

        for (int i = 0; i < 100; i++) {
            header += "\tcov_" + (i + 1);
        }
        header += "\tcircseq_count";

        bw.write(header + "\n");


        for (String s : structures.keySet()) {

            bw.write(structures.get(s) + "\n");
        }

        bw.close();

    }

    private static void writeToFile(Map<String, String> structures, String path) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(path));

        for (String s : structures.values()) {
            bw.write(s + "\n");
        }
        bw.close();

    }

    private static Map<String, Structure> get_structures(String path) throws IOException {
        Map<String, Structure> structures = new HashMap<>();
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line = "";

        while ((line = br.readLine()) != null) {
            Structure structure = new Structure(line);
            structures.put(structure.getCoords(), structure);
        }
        br.close();

        return structures;
    }

    private static Map<String, String> get_exons(String exonspath) throws IOException {
        Map<String, String> bed_like = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(exonspath));

        String line = "";

        while ((line = br.readLine()) != null) {
            String id = line.split("\t")[3];
            bed_like.put(id, line);
        }
        br.close();

        return bed_like;
    }

    public static void main(String[] args) {
        if (args.length == 4) {
            try {
                post_annotate(args[0], args[1], args[2], args[3]);

            } catch (IOException ex) {
                System.out.println("Class requires as input: \nString _ptes_bed, String _exons_bed, String _alus_bed, String _mirna_bed,"
                        + "String _name, String _exons_fa, int kmer, String _path\n\nOR\n\nptes_tsv_path, annotated_ptes.bed file and sam file if just adding annotations with a reference");
                Logger.getLogger(AnnotateCircRNAs.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                System.exit(0);
            }
        }

        try {
            _ptes_bed = args[0];
            _exons_bed = args[1];
            _name = args[2];
            _exons_fa = args[3];
            kmer = Integer.parseInt(args[4]);
            _path = args[5];



            Thread t2 = new Thread(new Runnable() {
                @Override
                public void run() {
                    try {
                        System.out.println(">>> Running screen for miRNA binding sites...\n");
                        mirnabs_search(_ptes_bed, _exons_bed, _path, _name);
                    } catch (IOException ex) {
                        Logger.getLogger(AnnotateCircRNAs.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            });

          
            t2.start();
            

            t2.join();
            
            merge_annotations(_path, _name);

        } catch (Exception ex) {
            System.out.println("Class requires as input: \nString _ptes_bed, String _exons_bed, String _alus_bed, String _mirna_bed,"
                    + "String _name, String _exons_fa, int kmer, String _path\n\nOR\n\nptes_tsv_path and annotated_ptes.bed file if just adding annotations with a reference");
            Logger.getLogger(AnnotateCircRNAs.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    
}
