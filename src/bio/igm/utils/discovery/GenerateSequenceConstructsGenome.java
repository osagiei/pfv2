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
import java.util.List;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author osagie izuogu
 */
public class GenerateSequenceConstructsGenome {

    String path, genome_path;
    Map<String, String> genome = new HashMap<String, String>();

    Map<String, Map<String, PutativeStructure>> structures = new HashMap<>();

    public GenerateSequenceConstructsGenome(String _path, String _gen_fasta) throws IOException {
        this.path = _path;
        this.genome = get_genome(_gen_fasta);

        get_putative_structures();
        get_canonical_structures();

    }

    private Map<String, String> get_genome(String _path) throws IOException {
        Map<String, String> genome = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(_path));
        String line = br.readLine();;
        int counter = 0;

        StringBuilder builder = new StringBuilder();
        List<String> chromosomes = new ArrayList<String>();

        while ((line = br.readLine()) != null){
             if (line.startsWith(">")){
                //counter++;
                if (counter % 2 == 0){
                  chromosomes.add(builder.toString());
                  builder = new StringBuilder();
                }
                counter++;
             }
             builder.append(line);
             builder.append(System.getProperty("line.separator"));

        }

        br.close();

        for (String sequence : chromosomes){

          String[] seqs = sequence.split(">"); //builder.toString().split(">");

          String seq = "";

          for(String s : seqs){
            if (StringUtils.isNotBlank(s)) {
                String[] oneSeq = s.trim().split("\n");
                String id = oneSeq[0].split(" ")[0];

                seq = StringUtils.join(oneSeq, "", 1, oneSeq.length);
                genome.put(id, seq.toUpperCase());
            }
          }
        }

        return genome;
    }

    private void generate_sequence_jobs(String _type) throws IOException {

        for (String chr : this.structures.keySet()) {
            Map<String, PutativeStructure> temp = structures.get(chr);
            String sequence = this.genome.get(chr);

            generate_sequences(sequence, temp, _type);
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
                generate_sequence_jobs("ptes");
                structures = new HashMap<>();
            }

        }
        //if (counter < 5000000) {
            generate_sequence_jobs("ptes");
            structures = new HashMap<>();
        //}
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
                generate_sequence_jobs("canonical");
                structures = new HashMap<>();
            }

        }
        //if (counter < 5000000) {
            generate_sequence_jobs("canonical");
            structures = new HashMap<>();
        //}
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
            new GenerateSequenceConstructsGenome(path, genomic_fasta);
        } catch (IOException ex) {
            System.out.println("Class takes as input: \nworking_dir,\tgenomic_fasta\n");
            Logger.getLogger(GenerateSequenceConstructsGenome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}