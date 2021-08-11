package bio.igm.entities;

/**
 *
 * @author osagie - 05/2013
 */
public class PutativeStructure {
    String chromosome, id, strand, sequence, splice_signal;
    int start1, stop1, start2, stop2, junction;
    
    
    public PutativeStructure(String line){
        chromosome = line.split("\t")[0];
        start1 = Integer.parseInt(line.split("\t")[1]);
        stop1 = Integer.parseInt(line.split("\t")[2]);
        start2 = Integer.parseInt(line.split("\t")[3]);
        stop2 = Integer.parseInt(line.split("\t")[4]);
        junction = Math.abs(stop1 - start1);
        id = line.split("\t")[5] + ":" + junction;
        strand = id.split("_")[1];
        
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public int getStart1() {
        return start1;
    }

    public void setStart1(int start1) {
        this.start1 = start1;
    }

    public int getStop1() {
        return stop1;
    }

    public void setStop1(int stop1) {
        this.stop1 = stop1;
    }

    public int getStart2() {
        return start2;
    }

    public void setStart2(int start2) {
        this.start2 = start2;
    }

    public int getStop2() {
        return stop2;
    }

    public void setStop2(int stop2) {
        this.stop2 = stop2;
    }

    public int getJunction() {
        return junction;
    }

    public void setJunction(int junction) {
        this.junction = junction;
    }

    public String getSplice_signal() {
        return splice_signal;
    }

    public void setSplice_signal(String splice_signal) {
        this.splice_signal = splice_signal;
    }
    
    
    
    
    
  
}
