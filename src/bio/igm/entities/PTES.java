/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.entities;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author osagie
 */
public class PTES {

    Map<String, Reads> reads = new HashMap();
    String id, splice_signal;
    String locus;
    int count;
    boolean confirmed;
    boolean spanned;

    public PTES(String id, int count) {
        setId(id);
        splice_signal = id.split(":")[3];
        //setLocus(id);

    }

    public void addRead(Reads read) {
        if (read.getTarget().equalsIgnoreCase(this.id)) {
            this.reads.put(read.getId(), read);
            this.count += 1;
        }
    }

    public Map<String, Reads> getReads() {
        return reads;
    }

    public void setReads(Map<String, Reads> reads) {
        this.reads = reads;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getLocus() {
        return locus;
    }

    public void setLocus(String locus) {
        this.locus = locus;
    }

    public boolean isConfirmed() {
        return confirmed;
    }

    public void setConfirmed(boolean confirmed) {
        this.confirmed = confirmed;
    }

    public boolean isSpanned() {
        return spanned;
    }

    public void setSpanned(boolean spanned) {
        this.spanned = spanned;
    }

    public String getSplice_signal() {
        return splice_signal;
    }

    public void setSplice_signal(String splice_signal) {
        this.splice_signal = splice_signal;
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }
    
}
