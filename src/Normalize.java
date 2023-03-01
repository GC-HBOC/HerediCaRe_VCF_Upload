import java.io.*;
import java.util.*;
import java.util.zip.DataFormatException;


class Normalize{

    static String referenceFastaPath = "resources\\ChromFa";
    public static String getStackTrace(final Throwable throwable) {
         final StringWriter sw = new StringWriter();
         final PrintWriter pw = new PrintWriter(sw, true);
         throwable.printStackTrace(pw);
         return sw.getBuffer().toString().replace("\n", "").replace("\r", "");
    }
    
    //input an arrayList of variants to split the lines with more than one alternative allele
        //First Copy the Arraylist and clear the original ArrayList. Then add the processed lines to the cleared ArrayList

	 public static String[] split_lines(ArrayList<String> lines, ArrayList<String> normlines){
		try{
			for (String line: lines) {

                if(!line.startsWith("#")){
                    String[] colms = line.split("\t");

					// if lines have to be break down
					if(colms[4].contains(",")){
						String[] alts = colms[4].split(",");
						String start= colms[0]+"\t"+colms[1]+"\t"+colms[2]+"\t"+colms[3]+"\t"; //(chr, pos, id, ref)
						int anzahl =alts.length;
						String[] outputs = new String[anzahl]; // contains lines to output
                        for(int i = 0; i<anzahl;i++) outputs[i]=start; //Add beginning ((chr, pos, id, ref)) 
						for(int i = 0; i<anzahl;i++) outputs[i]+=alts[i]+"\t"+colms[5]+"\t"+colms[6]+"\t"; //then alt, qual, filter
						// after that Infos (splitted corresponding to the nimber of output lines)
						String[] infos = colms[7].split(";");
						String classmt = (colms[7].contains("CLASS=")) ? "CLASS=" : "MT=" ;
						for(String inf: infos){
							String[] infs = get_all_comb(inf,anzahl,classmt);
							for(int i = 0; i<anzahl;i++) outputs[i]+= infs[i] + ";"; 
						}
                        //ending with genotype
						colms[9] = repair_gt(colms);
						for(int i = 0; i<anzahl;i++) outputs[i]=outputs[i].substring(0,outputs[i].length()-1)+"	";
						// add the rest
						for(int i = 8;i<colms.length;i++)
							for(int j = 0; j<anzahl;j++) outputs[j]+= colms[i]+"\t"; 
						for(int i = 0; i<anzahl;i++) outputs[i]=outputs[i].substring(0,outputs[i].length()-1);
						// write into the file 
						for(int i = 0; i<anzahl;i++) normlines.add(outputs[i]+"\n");
					}else
						normlines.add(line);
				}else
					normlines.add(line);
			}
        }catch (DataFormatException e){
            //e.printStackTrace(System.out);
            return new String[]{"2","Error: Number of CLASS/MT in vcf wrong: "+getStackTrace(e)};
        }catch (IndexOutOfBoundsException e){
            //e.printStackTrace(System.out);
            return new String[]{"4","ERROR: index out of bounds: "+getStackTrace(e)};
        }catch (Exception e) {
            //e.printStackTrace(System.out);
            //System.out.println(e.getMessage());
            
        }
        return new String[]{"0",null};
    }
    
    public static String repair_gt(String[] colms){
        if (colms[8].startsWith("GT")){
            String[] parts = colms[9].split(":",2);
            String p1 = (parts.length == 2)? parts[1]:"";
            String s = (parts[0].contains("|")) ? "|":"/";
            if(parts[0].contains(s)){
                String[] ar = parts[0].split(s);
                String a = "";
                a+= (ar[0].equals("0")) ? "0" : ".";
                a+=s;
                //a+= (Integer.parseInt(ar[1])>1) ? "1" : ar[1];
                a+= (ar[1].equals(".") || ar[1].equals("0")) ? ar[1] : "1";
                return  (!p1.equals("")) ? a+":"+ p1 : a;
            }else
                if(Integer.parseInt(parts[0])>1)
                    return  (!p1.equals("")) ? "1:"+ p1 : "1";
        }
        return colms[9];
    }
    
    // for every line in the input file normalize ref and alt
    // output to be written in a file
    public static int genenotation( ArrayList<String> header, ArrayList<String> normalized, String outfname, String reference){
        ArrayList<String> tempgetrennt = new ArrayList<String>(normalized);
		normalized.clear();
        BufferedWriter output=null;
        
        try {
            File outfile = new File(outfname);
            output = new BufferedWriter(new FileWriter(outfile));
            for(String line: header){
                output.write(line+"\n");
            }
            for(String line: tempgetrennt){
                if(!line.startsWith("#")){
                    String[] colms = line.split("\t");
                    String chrom = colms[0];
                    List<String> valid = Arrays.asList("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","M","X","Y");
					if(chrom.startsWith("chr") || chrom.startsWith("Chr")){
                        chrom = chrom.substring(3);
                    }
					if(chrom.equals("MT")){
						chrom = "M";
					}
                    colms[0] = chrom;
					
                    if(!valid.contains(chrom)){
						System.out.println("ERROR: invalid chrom value in vcf: "+chrom);
						output.close();
                        return 1;
                    }
                    
                    long pos = Long.valueOf(colms[1]);
                    //use algorithm, get surrounding from the method
                    String[] alleles=normalize(colms[3],colms[4],chrom, pos, reference); 
                    colms[1] = alleles[2];
                    colms[3] = alleles[0];
                    colms[4] = alleles[1];
                    
                    line = "";
                    for(String s: colms)
                        line += s+"	";
                    line = line.substring(0,line.length()-1);
                    output.write(line+"\n");
                    normalized.add(line);
                }
                else 
                    output.write(line+"\n");
            }
            output.close();
        }catch (Exception e) {
            try{
                output.close();
            }catch(Exception ex){}
			System.out.println("Error while closing file");
			e.printStackTrace();
            return 1;
        }
        return 0;
    }

    /*
     * help function for split_lines for infos(Column 7).
     * if number of class erroneous output an exception 
     */
    private static String[] get_all_comb(String inf,int anzahl, String classmt) throws DataFormatException{
        String[] out=new String[anzahl];
        if(inf.startsWith(classmt))
            if(inf.split("=")[1].split(",").length != anzahl)
                throw new DataFormatException();
        
        
        if(!inf.contains("=")){
            for(int i = 0; i<anzahl;i++) out[i] = inf;
            return out;
        }
        String[] tmp = inf.split("=");
        String[] re = tmp[1].split(",");
        if(re.length == 1){
            String retmp=re[0];
            re = new String[anzahl];
            for(int i =0; i<anzahl;i++) re[i]= retmp;
        }
        if(re.length != anzahl) {System.out.println("#### Erroneous number \n "+re.length+"  "+tmp[0]);return null;}
        for(int i =0; i<anzahl;i++)
            out[i]=tmp[0] + "=" + re[i];
        return out;
    }
    
    // delete all "N" right and left from ref and alt & replace "." to "" 
    private static Object[] prenormalize_allel(String ref, String alt, long pos){
        if(ref.equals(".")) ref = "";
        if(alt.equals(".")) alt = "";
        
        while(true){
            if(ref.equals("") || alt.equals("")) {
                break;
            }else if(ref.charAt(0) == 'N' && alt.charAt(0) == 'N'){
                ref = ref.substring(1);
                alt = alt.substring(1);
                pos += 1;
            }else if(ref.charAt(ref.length()-1) == 'N' && alt.charAt(alt.length()-1) == 'N'){
                ref = ref.substring(0, ref.length()-1);
                alt = alt.substring(0, alt.length()-1);
            }else{
                break;
            }
        }
        return new Object[]{ref, alt, pos};
    }

    // use algorithm to normalize 
    // algorithm: https://genome.sph.umich.edu/wiki/Variant_Normalization
    private static String[] normalize(String re,String al,String chrom, long p, String reference){
        if(re.equals(al)){ 
            return new String[]{re,al,p+""};
        }
        
        
        Object[] temp = prenormalize_allel(re, al, p);
        String ref = (String) temp[0];  
        String alt = (String) temp[1];
        long pos = (Long) temp[2];
        
        
        // get the surroundings
        long len = 100;
        String surr=null;//surrounding(chrom,pos,len).toUpperCase();
                
        // delete same chars at the end & add on the left if one is missing 
        boolean changes = true;

        while(changes){
            changes = false;
            if(ref.length() > 0 && alt.length() > 0 && 
              ref.charAt(ref.length() - 1) == alt.charAt(alt.length() - 1)){
                ref = ref.substring(0, ref.length()-1);
                alt = alt.substring(0, alt.length()-1);
                changes = true;
            }

            if(ref.equals("") || alt.equals("")){
                if(surr == null){
                    surr = surrounding(chrom,pos,len, reference).toUpperCase();
                }
                ref = surr.charAt(surr.length() - 1) + ref;
                alt = surr.charAt(surr.length() - 1) + alt;
                surr = surr.substring(0, surr.length()-1);
                if(surr.length()<2){
                    surr = surrounding(chrom,pos,len, reference).toUpperCase();
                }
                changes = true;
                pos--;
            }
        }
        // delete chars from the left so far as possible
        while(ref.charAt(0) == alt.charAt(0) && ref.length()>=2 && alt.length()>=2 ){
            ref = ref.substring(1);
            alt = alt.substring(1);
            pos++;
        }
        return new String[]{ref,alt,pos+""};
    }
    
    // look for base before the given position 
    private static String surrounding(String chrom,long postmp,long len, String reference){
        long pos = postmp-(len+1);
        pos += pos/50;
        if(pos < 0){ 
            pos =0;
            len=postmp-1;
        }
        try {
            File file = new File( referenceFastaPath + "\\" + reference + "\\chr"+chrom+".fa");
            BufferedReader reader = new BufferedReader(new FileReader(file));
            
            reader.readLine();
            reader.skip(pos); 
            String line = "";
            int i =0;
            while(i<len){
                char c= (char) reader.read();
                if(c == '\n') continue;
                line += c ;
                i++;
            }
            reader.close();
            return line;
        } catch (IOException e) {
            //e.printStackTrace();
            return null;
        }
    }
    
}
    