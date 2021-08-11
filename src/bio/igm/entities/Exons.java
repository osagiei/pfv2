package bio.igm.entities;

import org.apache.commons.lang3.Range;

/**
 *
 * @author osagie izuogu - 05/2013
 */
public class Exons {
  int order = -1;
  int start = -1;
  int stop = -1;
  String refid;
  Transcript transcript;
  Range<Integer> range;
  String sequence = "";
  int size = 0;

  public Exons(Transcript rna, int start, int stop) {
    this.refid = rna.getRefseq();
    this.transcript = rna;
    this.start = start;
    this.stop = stop;
    this.size = (stop - start);
  }
  public Exons(String id, String seq){ //used when converting exons.bed (obtained from UCSC) to full transcriptome reference
      if((id.startsWith("NM")) || (id.startsWith("NR"))){
          this.refid = id.split("_")[0] + "_" + id.split("_")[1]; //assumes refseq NM_002111_X
          order = Integer.parseInt(id.split("_")[2]);
      }else{
          this.refid = id.split("_")[0].replace(".", ":"); //assumes knowngene id
          order = Integer.parseInt(id.split("_")[1]);
      }
      
      this.sequence = seq;
      this.size = seq.length();
      
  }

  public int getOrder() {
    return this.order;
  }

  public void setOrder(int order) {
    this.order = order;
  }

  public Range<Integer> getRange() {
    this.range = Range.between(Integer.valueOf(getStart()), Integer.valueOf(getStop()));
    return this.range;
  }

  public void setRange(Range<Integer> range) {
    this.range = range;
  }

  public int getSize() {
    return this.size;
  }

  public void setSize(int size) {
    this.size = size;
  }

  public int getStart() {
    return this.start;
  }

  public void setStart(int start) {
    this.start = start;
  }

  public int getStop() {
    return this.stop;
  }

  public void setStop(int stop) {
    this.stop = stop;
  }

  public String getRefid() {
    return this.refid;
  }

  public void setRefid(String refid) {
    this.refid = refid;
  }

  public Transcript getTranscript() {
    return this.transcript;
  }

  public void setTranscript(Transcript transcript) {
    this.transcript = transcript;
  }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
  
}
