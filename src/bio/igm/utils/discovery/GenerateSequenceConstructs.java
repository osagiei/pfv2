/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.discovery;

import bio.igm.entities.PutativeStructure;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author osagie izuogu
 */
public class GenerateSequenceConstructs {

    String path, chroms_path;
    Map<String, String> chromosomes = new HashMap<String, String>();
   
    Map<String, Map<String, PutativeStructure>> structures = new HashMap<>();

    public GenerateSequenceConstructs(String _path, String _gen_fasta) throws IOException {
        this.path = _path;
        this.chroms_path = _gen_fasta;

        get_putative_structures();
        get_canonical_structures();

    }

    private void read_genomic_fasta_by_chromosome(String _type) throws IOException {

        for (String chr : structures.keySet()) {
            try {
                Map<String, PutativeStructure> temp = structures.get(chr);
                BufferedReader br = new BufferedReader(new FileReader(chroms_path + "/" + chr + ".fa"));
                String line = "";
                StringBuilder builder = new StringBuilder();
                while ((line = br.readLine()) != null) {
                    builder.append(line);
                    builder.append(System.getProperty("line.separator"));

                }
                br.close();
                String[] seqs = builder.toString().split(">");
                String seq = "";

                for (String s : seqs) {
                    if (StringUtils.isNotBlank(s)) {
                        String[] oneSeq = s.trim().split("\n");
                        String id = oneSeq[0];

                        seq = StringUtils.join(oneSeq, "", 1, oneSeq.length);
                        //chromosomes.put(id, seq.toUpperCase());
                    }
                }
                generate_sequences(seq, temp, _type);
                br.close();
            } catch (Exception e) {
                System.out.println(chr + " not in the directory provided");
            }

        }
    }

    private void get_putative_structures() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path + "/putative_structures.txt"));
        String line = "";
        int counter = 0;

        while ((line = br.readLine()) != null) {
            PutativeStructure structure = new PutativeStructure(line);
            //putative_structures.put(structure.getId(), structure);
            if (structures.containsKey(structure.getChromosome())) {
                structures.get(structure.getChromosome()).put(structure.getId(), structure);
                counter++;
            } else {
                Map<String, PutativeStructure> temp = new HashMap<>();
                temp.put(structure.getId(), structure);
                structures.put(structure.getChromosome(), temp);
                counter++;
            }
            if (counter % 5000000 == 0) {
                read_genomic_fasta_by_chromosome("ptes");
                structures = new HashMap<>();
            }

        }
        if (counter < 5000000) {
            read_genomic_fasta_by_chromosome("ptes");
            structures = new HashMap<>();
        }
        br.close();
        System.out.println("Finished reading structures.." + structures.size());

    }

    private void get_canonical_structures() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path + "/canonical_structures.txt"));
        structures = new HashMap<>();
        String line = "";
        int counter = 0;

        while ((line = br.readLine()) != null) {
            PutativeStructure structure = new PutativeStructure(line);

            if (structures.containsKey(structure.getChromosome())) {
                structures.get(structure.getChromosome()).put(structure.getId(), structure);
                counter++;
            } else {
                Map<String, PutativeStructure> temp = new HashMap<>();
                temp.put(structure.getId(), structure);
                structures.put(structure.getChromosome(), temp);
                counter++;
            }
            if (counter % 5000000 == 0) {
                read_genomic_fasta_by_chromosome("canonical");
                structures = new HashMap<>();
            }

        }
        if (counter < 5000000) {
            read_genomic_fasta_by_chromosome("canonical");
            structures = new HashMap<>();
        }
        br.close();
        System.out.println("Finished reading canonical structures.." + structures.size());
    }

    private void generate_sequences(String sequence, Map<String, PutativeStructure> temp_structures, String _type) throws IOException {
        if (_type.equalsIgnoreCase("ptes")) {
            BufferedWriter bw = new BufferedWriter(new FileWriter(path + "/Constructs.fa", true));

            for (PutativeStructure structure : temp_structures.values()) {
                try {
                    //removed 1 because of java's 0 index, different from bowtie's start pos
                    String seq1 = sequence.substring(structure.getStart1() - 1, structure.getStop1());
                    String seq2 = sequence.substring(structure.getStart2() - 1, structure.getStop2());

                    String signal = sequence.substring(structure.getStop1(), structure.getStop1() + 2)
                            + sequence.substring(structure.getStart2() - 3, structure.getStart2() - 1);



                    String seq = seq1 + seq2;

                    if (!structure.getId().contains("+")) {
                        seq = reverse_complement_sequence(seq2) + reverse_complement_sequence(seq1);
                        signal = reverse_complement_sequence(signal);
                    }
                    String id = structure.getId() + ":" + signal;



                    // if ((signal.equalsIgnoreCase("GTAG")) || (signal.equalsIgnoreCase("CTAC"))) {
                    bw.write(">" + id + "\n" + seq.toUpperCase() + "\n");
                    // }
                } catch (Exception e) {
                    System.out.println("Error generating sequence construct for : " + structure.getId() + 
                            " -- omitting");
                }

            }
            bw.close();
        } else {
            BufferedWriter bw = new BufferedWriter(new FileWriter(path + "/Can.fa", true));

            for (PutativeStructure structure : temp_structures.values()) {
                try {

                    String seq1 = sequence.substring(structure.getStart1() - 1, structure.getStop1());
                    String seq2 = sequence.substring(structure.getStart2() - 1, structure.getStop2());

                    String signal = sequence.substring(structure.getStop1(), structure.getStop1() + 2)
                            + sequence.substring(structure.getStart2() - 3, structure.getStart2() - 1);

                    String seq = seq1 + seq2;

                    if (!structure.getId().contains("+")) {
                        seq = reverse_complement_sequence(seq2) + reverse_complement_sequence(seq1);
                        signal = reverse_complement_sequence(signal);
                    }
                    String id = structure.getId() + ":" + signal;



                    if ((signal.equalsIgnoreCase("GTAG")) || (signal.equalsIgnoreCase("CTAC"))) {
                        bw.write(">" + id + "\n" + seq.toUpperCase() + "\n");
                    }
                } catch (Exception e) {
                    System.out.println("Error generating construct for " + structure.getId() + " -- omitting!");
                }

            }
            bw.close();
        }
    }

    private String reverse_complement_sequence(String sequence) {
        String rev_comp = "";
        Map<String, String> nucs = new HashMap<>();
        nucs.put("A", "T");
        nucs.put("T", "A");
        nucs.put("C", "G");
        nucs.put("G", "C");
        nucs.put("N", "N");

        for (int i = (sequence.length() - 1); i >= 0; i--) {
            String temp = "" + sequence.charAt(i);
            rev_comp += nucs.get(temp.toUpperCase());
        }
        return rev_comp;
    }

    public static void main(String[] args) {
        String path = args[0];
        String genomic_fasta = args[1];
        try {
            new GenerateSequenceConstructs(path, genomic_fasta);
        } catch (IOException ex) {
            System.out.println("Class takes as input: \nworking_dir,\tgenomic_fasta\n");
            Logger.getLogger(GenerateSequenceConstructs.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void read_genomic_fasta() throws IOException {
        FilenameFilter filter = new FilenameFilter() {
            String ext = ".fa";

            @Override
            public boolean accept(File file, String string) {
                return string.endsWith(ext);
            }
        };

        File f = new File(this.chroms_path);


        if (f.isDirectory()) {
            String[] list = f.list(filter);
            for (String p : list) {
                BufferedReader br = new BufferedReader(new FileReader(chroms_path + "/" + p));
                String line = "";
                StringBuilder builder = new StringBuilder();
                while ((line = br.readLine()) != null) {
                    builder.append(line);
                    builder.append(System.getProperty("line.separator"));

                }
                br.close();
                String[] seqs = builder.toString().split(">");
                String seq = "";

                for (String s : seqs) {
                    if (StringUtils.isNotBlank(s)) {
                        String[] oneSeq = s.trim().split("\n");
                        String id = oneSeq[0];

                        seq = StringUtils.join(oneSeq, "", 1, oneSeq.length);
                        chromosomes.put(id, seq.toUpperCase());
                    }
                }
                builder = null;
                System.out.println("Finished reading " + p + "\t" + chromosomes.size());
            }
        }

        //LOG.info("Finished loading exons\nSize: " + transcripts.size());
    }

    
}
