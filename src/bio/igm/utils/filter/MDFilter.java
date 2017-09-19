/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bio.igm.utils.filter;

import bio.igm.entities.PTES;
import bio.igm.entities.Reads;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author Osagie
 */
public class MDFilter {

    Map<String, PTES> ptes = new HashMap<String, PTES>();
    Map<String, List<Reads>> canonical = new HashMap<String, List<Reads>>();

    /**
     *
     * @param putative
     * @param reads
     * @param span
     * @param pid
     */
    public MDFilter(Map<String, PTES> putative, Map<String, Reads> reads, int span, double pid) {
        Reads read = null;
        PTES p = null;

        /*
         * For each read (accepted after previous filters), check that read aligns to
         * minimum junction span using Cigar and MD fields
         */
        for (String s : reads.keySet()) {
            read = (Reads) reads.get(s);
            int junction = read.getRefJunction() - read.getStart();
            if ((read.getEditDistance().split(":")[2].equalsIgnoreCase("0")) && (putative.containsKey(read.getTarget()))) {
                p = putative.get(read.getTarget());
                p.addRead(read);
                p.setCount(p.getReads().size());
                p.setConfirmed(true);

                putative.put(p.getId(), p);
            } else {
                String[] parsedMd = parseMD(read.getMdfield(), read.getCigar(), junction);
                if ((checkJunctionSpan(read, parsedMd, span, pid))
                        && (putative.containsKey(read.getTarget()))) {
                    p = putative.get(read.getTarget());
                    p.addRead(read);
                    p.setCount(p.getReads().size());
                    p.setSpanned(true);

                    putative.put(p.getId(), p);

                }
            }
        }
        this.ptes = putative;

    }

    public MDFilter(Map<String, Reads> reads, int span, double pid, String path) throws IOException {
        Reads read = null;
        List<Reads> temp = null;
        Map<String, Integer> counts = new HashMap<String, Integer>();

        BufferedWriter bpR = new BufferedWriter(new FileWriter(path + "PTESReads"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(path + "ptescounts.tsv"));
        BufferedWriter bwB = new BufferedWriter(new FileWriter(path + "ptescounts.tsv.bed"));
        BufferedWriter bwJ = new BufferedWriter(new FileWriter(new StringBuilder().append(path).append("junctions.fa").toString()));
        BufferedWriter bwP = new BufferedWriter(new FileWriter(new StringBuilder().append(path).append("pid.tsv").toString()));

        BufferedWriter bf = new BufferedWriter(new FileWriter(path + "junctional-filtered.sam"));
        bwP.write("Read_ID\tPTES_ID\tEdit_Distance\tLeftPID\tRightPID\n");
        /*
         * For each read (accepted after previous filters), check that read aligns to
         * minimum junction span using Cigar and MD fields
         */
        for (String s : reads.keySet()) {
            read = (Reads) reads.get(s);
            int junction = read.getRefJunction() - read.getStart();
            String[] parsedMd = parseMD(read.getMdfield(), read.getCigar(), junction);
            if (((checkJunctionSpan(read, parsedMd, span, pid))) || (read.getEditDistance().equalsIgnoreCase("NM:i:0"))) {
            //if (checkJunctionSpan(read, parsedMd, span, pid)) {
                bpR.write(new StringBuilder().append(read.getLine()).append("\n").toString());
                bwP.write(new StringBuilder().append(read.getId()).append("\t").append(read.getTarget()).append("\t").append(read.getEditDistance()).append("\t").append(read.getLeftpid()).append("\t").append(read.getRightpid()).append("\n").toString());
                bwJ.write(new StringBuilder().append(">").append(read.getId()).append("\t").append(read.getTarget()).append("\tStart: ").append(read.getStart()).append("\tJunction:").append(read.getRefJunction() - read.getStart() + read.getJunctionShift()).append("\t").append(read.getJunctionSeq()).append("\t").append(read.getJunctionShift()).append("\t").append(read.getHex()).append("\t").append(read.getEditDistance()).append("\t").append(read.getMdfield()).append("\t").append(read.getCigar()).append("\n").append(read.getSequence()).append("\n").append(read.getMdTransformed()).append("\n").toString());

                if (counts.containsKey(read.getTarget())) {
                    int x = counts.get(read.getTarget()) + 1;
                    counts.put(read.getTarget(), x);
                } else {
                    counts.put(read.getTarget(), 1);
                }


            } else {
                bf.write(read.getLine() + "\n");
            }
        }

        for (String s : counts.keySet()) {
            bw.write(s + "\t" + counts.get(s) + "\n");
            try {
                String[] str = s.split(":"); //chr4:144464659-144465123_+:41:GGTC
                bwB.write(str[0] + "\t" + str[1].split("-")[0] + "\t" + str[1].split("-")[1].split("_")[0]
                        + "\t" + s + "\t" + counts.get(s) + "\t" + str[1].split("_")[1] + "\t" + str[3] + "\n");
            } catch (Exception e) {
                System.out.println(s);
            }
        }

        bpR.close();
        bwJ.close();
        bwP.close();
        bf.close();
        bw.close();
        bwB.close();

    }

    public MDFilter(String path, int span, double pid) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path + "canonical.sam"));
        Map<String, Integer> counts = new HashMap<String, Integer>();

        BufferedWriter bcR = new BufferedWriter(new FileWriter(path + "flanking-canonical-reads.sam"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(path + "flanking-canonical-counts.tsv"));
        BufferedWriter bwB = new BufferedWriter(new FileWriter(path + "flanking-canonical-counts.tsv.bed"));
        String line = "";

        while ((line = br.readLine()) != null) {
            Reads read = new Reads(line);
            int junction = read.getRefJunction() - read.getStart();
            String[] parsedMd = parseMD(read.getMdfield(), read.getCigar(), junction);
            if ((checkJunctionSpan(read, parsedMd, span, pid))) {
                bcR.write(new StringBuilder().append(read.getLine()).append("\n").toString());
                if (counts.containsKey(read.getTarget())) {
                    int x = counts.get(read.getTarget()) + 1;
                    counts.put(read.getTarget(), x);
                } else {
                    counts.put(read.getTarget(), 1);
                }
            }

        }
        for (String s : counts.keySet()) {
            bw.write(s + "\t" + counts.get(s) + "\n");
            try {
                String[] str = s.split(":"); //chr4:144464659-144465123_+:41:GGTC
                bwB.write(str[0] + "\t" + str[1].split("-")[0] + "\t" + str[1].split("-")[1].split("_")[0]
                        + "\t" + s + "\t" + counts.get(s) + "\t" + str[1].split("_")[1] + "\t" + str[3] + "\n");
            } catch (Exception e) {
                System.out.println(s);
            }
        }
        bw.close();
        bcR.close();
        br.close();
        bwB.close();
    }

    private String[] parseMD(String md, String cigar, int junction) {
        StringBuilder seq = new StringBuilder();
        String temp = "";

        String pattern = "[0-9]+";
        String[] contents = md.split(":");

        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(contents[2]);

        int counter = 0;
        List starts = new ArrayList();
        List ends = new ArrayList();
        List readGroups = new ArrayList();

        while (m.find()) {
            starts.add(Integer.valueOf(m.start()));
            ends.add(Integer.valueOf(m.end()));

            readGroups.add(Integer.valueOf(Integer.parseInt(m.group())));

            counter++;
        }

        contents[2] = new StringBuilder().append(contents[2]).append(" ").toString();
        for (int i = 0; i < readGroups.size(); i++) {
            for (int j = 0; j < ((Integer) readGroups.get(i)).intValue(); j++) {
                seq.append("M");
            }
            seq.append(contents[2].charAt(((Integer) ends.get(i)).intValue()));
        }
        temp = seq.toString().replaceAll("[^a-zA-Z0-9]", "");
        temp.trim();

        return parseIndels(temp, cigar, junction);
    }

    /**
     *
     * @param seq
     * @param cigar
     * @param junction
     * @return
     */
    public String[] parseIndels(String seq, String cigar, int junction) {
        List result = new ArrayList();
        List minShift = new ArrayList();
        int shift = 0;
        int unshift = 0;
        int counter2 = 0;
        Map index = new HashMap();

        List ind = new ArrayList();
        StringBuilder sb = new StringBuilder(seq);
        String pattern = "[A-Za-z]";

        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(cigar);

        while (m.find()) {
            char c = m.group().charAt(0);
            index.put(Integer.valueOf(m.end()), Character.valueOf(c));

            ind.add(Integer.valueOf(m.end()));
        }

        int counter = 0;
        char c;
        int a;
        for (Iterator i$ = ind.iterator(); i$.hasNext();) {
            int z = ((Integer) i$.next()).intValue();
            int pos = z;
            c = ((Character) index.get(Integer.valueOf(z))).charValue();
            switch (c) {
                case 'I':
                    result = processDetail(cigar, c, pos);
                    break;
                case 'D':
                    result = processDetail(cigar, c, pos);
                    break;
                case 'S':
                    result = processDetail(cigar, c, pos);
                    c = 's';
                    break;
                default:
                    result = processDetail(cigar, c, pos);
            }

            if (c != 'M') {
                for (a = 1; a < result.size();) {
                    for (int i = 0; i < ((Integer) result.get(a)).intValue(); i++) {
                        sb.insert(((Integer) result.get(a - 1)).intValue() + i, c);
                        if (c == 'D') {
                            counter--;
                        }
                        counter++;
                        counter2++;
                    }

                    if (((Integer) result.get(a - 1)).intValue() < junction) {
                        shift = counter;
                        unshift = counter2;
                    }
                    a += 2;
                }
            }
        }

        String[] res = {sb.toString(), new StringBuilder().append(shift).append("").toString(), new StringBuilder().append(unshift).append("").toString()};

        return res;
    }

    private List<Integer> processDetail(String cigar, char c, int pos) {
        List result = new ArrayList();

        int start = 0;
        int width = 0;

        Pattern p = Pattern.compile("[0-9]+");
        Matcher m = p.matcher(cigar.substring(0, pos));
        start = 0;

        while (m.find()) {
            start += Integer.parseInt(m.group());
            width = Integer.parseInt(m.group());
        }
        result.add(Integer.valueOf(start - width));
        result.add(Integer.valueOf(width));

        return result;
    }

    private boolean checkJunctionSpan(Reads read, String[] parsedMD, int span, double pid) {
        boolean junctionSpan = false;
        int shift = Integer.parseInt(parsedMD[1]);
        int unshift = Integer.parseInt(parsedMD[2]);
        int oldJunction = read.getRefJunction();
        String leftMD = "";
        String rightMD = "";
        double leftpid = 0.0D;
        double rightpid = 0.0D;
        read.setJunctionShift(shift);
        read.setHex(shift);

        int position = oldJunction - read.getStart() + shift + 1;
        junctionSpan = (position > (read.getStart() + (span / 2))) && ((position + (span / 2)) < (read.getSequence().length() - 2));
        int counter = 0;
        if (junctionSpan) {
            int startpos = position - (span / 2);
            int lastpos = position + (span / 2);

            read.setJunctionSeq(parsedMD[0].substring(startpos - 1, lastpos + 1));
            for (int i = startpos; i <= lastpos; i++) {
                if (parsedMD[0].charAt(i) != 'M') {
                    junctionSpan = false;
                }

                counter++;
            }
            StringBuilder parsed = new StringBuilder(parsedMD[0]);
            leftMD = parsedMD[0].substring(0, startpos);
            rightMD = parsedMD[0].substring(lastpos);

            int matchCount = StringUtils.countMatches(leftMD, "M");
            leftpid = (double) matchCount / leftMD.length();
            matchCount = StringUtils.countMatches(rightMD, "M");
            rightpid = (double) matchCount / rightMD.length();

            parsed.insert(startpos - 1, "<");
            parsed.insert(lastpos + 1, ">");
            read.setLeftpid(String.format("%.3f", new Object[]{Double.valueOf(leftpid)}));
            read.setRightpid(String.format("%.3f", new Object[]{Double.valueOf(rightpid)}));
            read.setMdTransformed(new StringBuilder().append(parsed.toString()).append("|").append(counter).append("|").append(position).append("|").append(startpos).append("|").append(lastpos).append("|").append(String.format("%.3f", new Object[]{Double.valueOf(leftpid)})).append("|").append(String.format("%.3f", new Object[]{Double.valueOf(rightpid)})).toString());

            if ((rightpid < pid) || (leftpid < pid)) {
                junctionSpan = false;
            }
        }
        return junctionSpan;
    }

    /**
     *
     * @return
     */
    public Map<String, PTES> getPtes() {
        return this.ptes;
    }

    public Map<String, List<Reads>> getCanonical() {
        return canonical;
    }

    /*
     public static void main(String[] args) {
     String s = "chr4:144464659-144465123_+:41:GGTC";
     String[] str = s.split(":"); //chr4:144464659-144465123_+:41:GGTC
     System.out.println(str[0] + "\t" + str[1].split("-")[0] + "\t" + str[1].split("-")[1].split("_")[0]
     + "\t" + s + "\t" + str[1].split("-")[1].split("_")[1] + "\t"+ str[3] + "\n");
     }
     * */
}
