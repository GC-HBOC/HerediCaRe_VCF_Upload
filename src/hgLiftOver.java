import java.io.*;
import java.util.*;

//chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 

/*
 *`hgLiftOver: contains static methods to read a vcf file and LiftOver to the corresponding reference genome
 * 
 * Workflow:
 *  1. encoding_and_tohg19() is called by ParseVCF main method
 *  2. read_and_encoding() is called (returns list of variant lines with new coordinates)
 *  3. if the vcf file is hg38 it is transformed to hg19 in liftOver()
 *     (line by line using the corresponding chain file)
*/
public class hgLiftOver{

    /*
     * is called by main method in ParseVCF
    */
    static String chainFileHg19 = "resources\\hg19ToHg38.over.chain";
    static String chainFileHg38 = "resources\\hg38ToHg19.over.chain";

    public static ArrayList<String> encoding_and_liftover(ArrayList<String> lines, String genome) throws Exception{
        

        switch(genome.toUpperCase()){
            case "HG38":
			    System.out.println("LiftOver from hg38 to hg19");
			    lines = liftOver(lines, chainFileHg38);
                break;
            case "HG19":
                System.out.println("LiftOver from hg19 to hg38");
			    lines = liftOver(lines, chainFileHg19);
                break;
        }
        return lines;
    }    
	
	/*
     * change each line to hg19/hg38 if needed
     * chrom and pos are changed
    */
	public static ArrayList<String> liftOver(ArrayList<String> lines, String chainFileName) throws Exception{
		ArrayList<String> out = new ArrayList<String>();
        for(String line: lines) {
			if(line.startsWith("#")){
				out.add(line);
				continue;
			}
			String[] columns= line.split("\t",3);
			String chrom = columns[0];
			int pos = Integer.parseInt(columns[1]);
			ArrayList<String> chain = findChain(chainFileName, chrom, pos);
            try{
                String newchrom = fixChr(chain.get(0).split("\\s+")[7]);
                int newpos = searchInChain(chain,pos);
					
                out.add(newchrom + "\t" + newpos + "\t" + columns[2]);
            }catch(IndexOutOfBoundsException e){
                throw new Exception("Error with the following Variant\n"+line+"\nPosition not available in the other reference genome.");
            } 
			

		}
		return out;
	}


    /*
     * for a chrom and pos: finds the relevant chain entry of teh corresponding chain file
     * returns the entry as an ArrayList<String> of lines
    */
	public static ArrayList<String> findChain(String chainFileName, String chrom, int chrompos){
        ArrayList<String> outchain = new ArrayList<String>();
        BufferedReader br = null;
        try{
			br = new BufferedReader(new FileReader(chainFileName));
            String line;
            while ((line = br.readLine()) != null) {
                if(!line.contains("chain"))
                    continue;

                String[] infos = line.split("\\s+");
                
                if(infos[2].endsWith(chrom)){
                    int startpos = Integer.parseInt(infos[5]);
                    int endpos = Integer.parseInt(infos[6]);
                    if(startpos - 1 <= chrompos && chrompos <= endpos){
                        //System.out.println(line);
                        outchain.add(line);
                        while ((line = br.readLine()) != null) {
                            if(line.equals("") || line.contains("chain"))
                                break;
                            outchain.add(line);
                        }
                    }
                }
            }
			
        }catch(IOException e){
			e.printStackTrace(System.out);
		} finally{
            if(br!=null){
                try {
                    br.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }
        return outchain;
    }

    /*
     * Is given a chain block and a position
     * finds the hg19/hg38 position using the information of the chain block
    */
	public static int searchInChain(ArrayList<String> chain, int chrompos){
        String[] infos = chain.get(0).split("\\s+");
        int tpos = Integer.parseInt(infos[5]);
        int qpos = Integer.parseInt(infos[10]);
        for(int i = 1; i< chain.size() -1; i++){
        
            String[] blockinfo = chain.get(i).split("\t");
            int size = Integer.parseInt(blockinfo[0]);
            int dt = Integer.parseInt(blockinfo[1]);
            int dq = Integer.parseInt(blockinfo[2]);
            tpos += size; // matches
            qpos += size;
            
            if(tpos >= chrompos)
                return qpos - (tpos-chrompos);
            
            tpos += dt; // gaps
            qpos += dq;
            
            if(tpos >= chrompos)
                return -1;
            
            
        }

        //last line
		String blockinfo = chain.get(chain.size() -1);
		int size = Integer.parseInt(blockinfo);
		
		tpos += size; // matches
        qpos += size;
            
        if(tpos >= chrompos)
            return qpos - (tpos-chrompos);
            
        return -2; //should not happen
    }

    private static String fixChr(String chrom){
        
        if (chrom.contains("chr")){
            return chrom.substring(3);
        } else{ return chrom;}
    }

}

