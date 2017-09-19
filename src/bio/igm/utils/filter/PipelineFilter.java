/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.filter;

import bio.igm.entities.PTES;
import bio.igm.entities.Reads;
import bio.igm.utils.init.Logging;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author osagie
 */
public class PipelineFilter {

    int jspan = 0;
    double pid = 0.00;
    Map<String, PTES> putative_ptes = new HashMap<String, PTES>();
     Map<String, PTES> canonical = new HashMap<String, PTES>();

    Map<String, Reads> processed_ptes_reads = new HashMap<String, Reads>();
    Map<String, String> constructs = new HashMap<String, String>();
    
     Map<String, Reads> raw_canonical_reads = new HashMap<String, Reads>();
    Map<String, Reads> processed_canonical_reads = new HashMap<String, Reads>();
    
    String path;
    boolean filters = true;
    boolean refseq = true;
    boolean genomic = true;
    public static Logger LOG;

    public PipelineFilter(String _path, int _jspan, double _pid, boolean _filters, boolean _genomic, boolean _refseq) throws IOException {
        this.path = _path;
        this.filters = _filters;
        this.pid = _pid;
        this.jspan = _jspan;
        this.filters = _filters;
        this.refseq = _refseq;
        this.genomic = _genomic;
        File f = new File(_path);

        try {
            if (f.isDirectory()) {
                LOG = new Logging(_path, PipelineFilter.class.getName()).setup();
            } else {
                LOG = new Logging(f.getParent(), PipelineFilter.class.getName()).setup();
            }
        } catch (IOException ex) {
            Logger.getLogger(PipelineFilter.class.getName()).log(Level.SEVERE, null, ex);
        }
        LOG.info("Reading ptes.sam and applying filters..\nInput Parameters:\n\tPID -> " + _pid + "\n\tJSpan -> "
                + _jspan);
      


        if (this.filters) { //default 
            this.refseq = false;
            this.genomic = false;
           
            processed_ptes_reads = apply_filters();

        } else if (this.genomic || this.refseq) { //using default settings, both booleans should be false here
            
            processed_ptes_reads = apply_filters();

        } else {

            processed_ptes_reads = get_sam_reads("ptes");
        }

      
        System.gc();

        LOG.info("Applying Minimum Junction Span and Construct Percent Identity Filters .. ");
        putative_ptes = new MDFilter(processed_ptes_reads, jspan, pid, this.path).getPtes();

        processed_ptes_reads = null;
        System.gc();

        LOG.info("Processing reads mapped to Canonical Junctions ..");
        canonical = new MDFilter(this.path, jspan, pid).getPtes();
        processed_canonical_reads = null;



    }

    private Map<String, Reads> get_sam_reads(String filename) throws IOException {
        Map<String, Reads> temp = new HashMap<String, Reads>();
        
        BufferedReader br = new BufferedReader(new FileReader(this.path + "unique.sam"));

        String line = "";

        while ((line = br.readLine()) != null) {
            Reads r = new Reads(line);
            temp.put(r.getId(), r);
        }

        return temp;
    }

    private Map<String, Reads> apply_filters(String type, Map<String, Reads> _structures) throws IOException {
        Map<String, Reads> temp = _structures;

        if (type.equalsIgnoreCase("ptes")) {
            if (this.genomic) { //default = false

                LOG.info("[ Applying Genomic Filter Only.. ");
                temp = new FilterGenomicHits(path).getReads();
            } else if (this.refseq) { // default = false
                LOG.info("Applying Transcriptomic Filter Only.. ");
                return new FilterTranscriptomicHits(temp, path).getReads(); //temp here contains either raw ptes reads or reads already filtered using genomic filter
            } else {
                LOG.info("Applying Genomic & Transcriptomic Filters ..");
                return new FilterTranscriptomicHits(new FilterGenomicHits(path).getReads(), this.path).getReads(); //used for both filters; default!!
            }

        }

        return temp; //expected to return canonical reads if line encountered
    }

    private Map<String, Reads> apply_filters() throws IOException {
        Map<String, Reads> temp = new HashMap<String, Reads>();


        if (this.genomic) { //default = false

            LOG.info("[ Applying Genomic Filter Only.. ");
            temp = new FilterGenomicHits(path).getReads();
        } else if (this.refseq) { // default = false
            LOG.info("Applying Transcriptomic Filter Only.. ");
            return new FilterTranscriptomicHits(temp, path).getReads(); //temp here contains either raw ptes reads or reads already filtered using genomic filter
        } else {
            LOG.info("Applying Genomic & Transcriptomic Filters ..");
            new FilterGenomicHits(path);
            System.out.print("Finished applying genomic filter...\n");
            //return new FilterTranscriptomicHits(new FilterGenomicHits(path).getReads(), this.path).getReads(); //used for both filters; default!!
            return new FilterTranscriptomicHits(get_sam_reads(path), this.path).getReads();
        }



        return temp; //expected to return canonical reads if line encountered
    }

    public static void main(String[] args) {
        String path = args[0];
        int _jspan = 8;
        double _pid = 0.85;
        boolean all_filters = true;
        boolean genomic = true;
        boolean refseq = true;

        try {

            _jspan = Integer.parseInt(args[1]);
            _pid = Double.parseDouble(args[2]);
            all_filters = args[3].equalsIgnoreCase("0") ? false : true;
            genomic = args[4].equalsIgnoreCase("0") ? false : true;
            refseq = args[5].equalsIgnoreCase("0") ? false : true;
        } catch (Exception e) {
            LOG.severe("\nThere are errors in your input parameters. Please re-check before running again.\nProceeding with default values:\n"
                    + ">>> Run all filters = Yes\n\t PID = 0.85 \n\t Junction Span = 8\n\n");
        }
        try {
            new PipelineFilter(path, _jspan, _pid, all_filters, genomic, refseq);
        } catch (IOException ex) {
            LOG.log(Level.SEVERE, null, ex);
        }

    }
}
