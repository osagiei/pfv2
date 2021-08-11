package bio.igm.entities;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author osagie izuogu - 05/2013
 */
public class Transcript {
  String refseq;
  Map<Integer, Exons> exons = new HashMap();
  int numOfExons = 0;
  String sequence;

  public Transcript(String id)
  {
    this.refseq = id;
  }
  public void addExon(Exons e) {
    this.exons.put(Integer.valueOf(e.getOrder()), e);
  }

  public Map<Integer, Exons> getExons() {
    return this.exons;
  }

  public void setExons(Map<Integer, Exons> exons) {
    this.exons = exons;
  }

  public int getNumOfExons() {
    this.numOfExons = this.exons.size();
    return this.numOfExons;
  }

  public void setNumOfExons(int numOfExons) {
    this.numOfExons = numOfExons;
  }

  public String getRefseq() {
    return this.refseq;
  }

  public void setRefseq(String refseq) {
    this.refseq = refseq;
  }

  public String getSequence() {
    return this.sequence;
  }

  public void setSequence(String sequence) {
    this.sequence = sequence;
  }
}
