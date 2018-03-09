/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.discovery;

import bio.igm.utils.init.Logging;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author osagie izuogu
 */
public class ProcessShuffledCoordinates {

    String path;
    int chrom_width, segment_size;
    private static Logger LOG;
   
    public ProcessShuffledCoordinates(String _path, int _chrom_width, int _segment_size) throws IOException {
        this.path = _path;

        //gsc = new GenerateSequenceConstructs(_path, _chrom_paths);
        File f = new File(_path);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_path, ProcessShuffledCoordinates.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), ProcessShuffledCoordinates.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(ProcessShuffledCoordinates.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.chrom_width = _chrom_width;
        this.segment_size = _segment_size;
        Thread t = new Thread(new Runnable() {

            @Override
            public void run() {
                try {
                    read_star_chimeric_junction();
                } catch (IOException ex) {
                    LOG.info("Error reading STAR Chimeric junction file.. ");
                }
            }
        });
        
        t.start();        
        
        read_star_canonical_junctions();
    }

    private void read_star_chimeric_junction() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path + "/star_Chimeric.out.junction"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(path + "/putative_structures.txt"));


        String line = "";


        LOG.info("Reading STAR Chimeric junction file.. ");
        while ((line = br.readLine()) != null) {

            try {
                String chrL = line.split("\t")[0];
                String chrR = line.split("\t")[3];
                String oL = line.split("\t")[2];
                String oR = line.split("\t")[5];
                String id = "";
                String to_print = "";

                int left = Integer.parseInt(line.split("\t")[1]);
                int right = Integer.parseInt(line.split("\t")[4]);
                int tl = 0;
                int tr = 0;

                if (!chrL.equalsIgnoreCase(chrR)) {
                    continue;
                }
                if (chrL.equalsIgnoreCase("chrM") || chrL.equalsIgnoreCase("chrMT")) {
                    continue;
                }
                if (!oL.equalsIgnoreCase(oR)) {
                    continue;
                }
                if (Math.abs(left - right) > chrom_width) {
                    continue;
                }
                if(Math.abs(left -right) > 1000000 || Math.abs(left -right) < 50) {
                  continue;
                }
                if (Integer.parseInt(line.split("\t")[6]) <= 0) {
                    continue;
                }

                if (oL.contains("-")) {
                    right -= 1;
                    left += 1;
                    if (left > right) {
                        continue;
                    } else {
                        id = chrL + ":" + (left - 1) + "-" + right;
                        tr = right - segment_size;
                        tl = left + segment_size;
                        to_print = chrL + "\t" + tr + "\t" + right + "\t" + left + "\t" + tl + "\t" + id + "_" + oL;
                    }
                } else {
                    left -= 1;
                    right += 1;
                    if (right > left) {
                        continue;
                    } else {
                        id = chrL + ":" + (right - 1) + "-" + left;
                        tr = right + segment_size;
                        tl = left - segment_size;
                        to_print = chrL + "\t" + tl + "\t" + left + "\t" + right + "\t" + tr + "\t" + id + "_" + oL;
                    }
                }

                bw.write(to_print + "\n");



            } catch (Exception e) {
                LOG.info("Error processing this read: " + line);
            }

        }


        br.close();
        bw.close();
        LOG.info("Finished reading processed SAM file .. ");

    }

    private void read_star_canonical_junctions() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path + "/star_SJ.out.tab"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(path + "/canonical_structures.txt"));


        String line = "";


        LOG.info("Reading STAR canonical junctions file.. ");
        while ((line = br.readLine()) != null) {
            String chrL = line.split("\t")[0];

            String id = "";
            String to_print = "";

            int left = Integer.parseInt(line.split("\t")[1]);
            int right = Integer.parseInt(line.split("\t")[2]);
            int tr = 0;
            int tl = 0;

            if ((Integer.parseInt(line.split("\t")[4])) > 2) {
                continue; //non-canonical splice site
            }

            int o = Integer.parseInt(line.split("\t")[3]);


            switch (o) {
                case 1:
                    left -= 1;
                    right += 1;

                    id = chrL + ":" + left + "-" + (right - 1);
                    tr = right + segment_size;
                    tl = left - segment_size;
                    to_print = chrL + "\t" + tl + "\t" + left + "\t" + right + "\t" + tr + "\t" + id + "_+";
                    break;
                case 2:
                    left -= 1;
                    right += 1;

                    id = chrL + ":" + left + "-" + (right - 1);
                    tr = right + segment_size;
                    tl = left - segment_size;
                    to_print = chrL + "\t" + tl + "\t" + left + "\t" + right + "\t" + tr + "\t" + id + "_-";
                    break;
                default:
                    continue;
            }
            bw.write(to_print + "\n");
        }
        br.close();
        bw.close();
    }

    public static void main(String[] args) {
        String path = args[0];
        int _chrom_width = Integer.parseInt(args[1]);

        int segment_size = Integer.parseInt(args[2]);
        try {
            new ProcessShuffledCoordinates(path, _chrom_width, segment_size);
        } catch (IOException ex) {
            System.out.println("Class requires as input:\nworking_directory,\tmax_width_between_anchors,\t"
                    + "segment_size\n");
            Logger
                    .getLogger(ProcessShuffledCoordinates.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
    }
}
