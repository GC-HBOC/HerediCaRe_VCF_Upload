import java.io.*;
import java.util.*;
import java.nio.charset.*;

//chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 

/*
 *`hgLiftOver: contains static methods to read a vcf file and to deal with different file encodings
 * 
 * Workflow:
 *  1. encoding_and_tohg19() is called by ParseVCF main method
 *  2. read_and_encoding() is called (returns lines of read file)
 *      2.1 encoding of file is found in get_encoding() using the file command stored in file-5.03-bin\\bin\\file
 *      2.2 vcf file is read line by line using the found encoding
 *  3. if the vcf file is hg38 it is transformed to hg19 in liftOver()
 *     (line by line using hg38ToHg19.over.chain)
*/
public class hgLiftOver{

    /*
     * is called by main method in ParseVCF
    */
    static String chainFileHg19 = "C:\\Users\\Rudel\sFankep\\Documents\\GitHub\\HerediCaRe_VCF_Upload\\hg19ToHg38.over.chain";
    static String chainFileHg38 = "C:\\Users\\Rudel\sFankep\\Documents\\GitHub\\HerediCaRe_VCF_Upload\\hg38ToHg19.over.chain";

    // To delete
    public static ArrayList<String> encoding_and_liftover(String VCF_IN_DIR, String file){
		String infname = VCF_IN_DIR+"\\"+file;
        //encoding
		ArrayList<String> lines = read_and_encoding(infname);
		if(lines == null)
			return null;
		
        //hg19
		if(file.startsWith("hg38-")){
			System.out.println("change from hg38 to hg19");
			lines = liftOver(lines, chainFileHg38);
		} else{
            System.out.println("change from hg19 to hg38");
			lines = liftOver(lines, chainFileHg19);
        }
        return lines;
    }
    
    // To keep
    public static ArrayList<String> encoding_and_liftover(ArrayList<String> lines, String genome){
        //hg19
        switch(genome.toUpperCase()){
            case "HG38":
			    System.out.println("change from hg38 to hg19");
			    lines = liftOver(lines, chainFileHg38);
                break;
            case "HG19":
                System.out.println("change from hg19 to hg38");
			    lines = liftOver(lines, chainFileHg19);
                break;
        }
        return lines;
    }
	
    /*
     * reads file using the encoding
     * returns file as ArrayList<String> contaning the lines
    */
	public static ArrayList<String> read_and_encoding(String infname){
		ArrayList<String> lines = new ArrayList<String>();
		
        try{
            String encoding = get_encoding(infname);
            if (encoding.startsWith("Error")){
                return null;
            }
            
            File file = new File(infname);
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file), Charset.forName(encoding)));
            
            String first = String.valueOf((char) reader.read());
            Boolean firstchar = ("#".equals(first) ? true:false);
            
            String line = null;
            while ((line = reader.readLine()) != null) {
                if(line.startsWith("\"")) // " entfernen
                    line = line.substring(1, line.length()-1);
                if(line.equals("")) // leerzeilen überspringen
                    continue;
                
                if(firstchar){
                    line= "#"+line;
                    firstchar=false;
                }
				lines.add(line);
			}
			reader.close();
        }catch (Exception e) {
			return null;
		}
        return lines;
	}
	    
    /*
     * get the encoding of the vcf file using file-5.03-bin\\bin\\file.exe
    */
	public static String get_encoding(String fname){
        String errorout= "";
        try{
            //Process p = Runtime.getRuntime().exec("file-5.03-bin\\bin\\file.exe --mime-encoding "+fname);
            Process p = Runtime.getRuntime().exec("file-5.03-bin\\bin\\file.exe --mime-encoding "+fname);
            p.waitFor();
            BufferedReader pout = new BufferedReader(new InputStreamReader(p.getInputStream()));
            
            //catch error steam
            BufferedReader perr = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String line;
            errorout= "Errorout=";
            while ((line = perr.readLine()) != null) {
                if (line.isEmpty()) {
                    break;
                }
                errorout += line;
            }
            if (errorout.contains("Usage: "))
                errorout = errorout.split("Usage: ")[0];
            //System.out.println(errorout);
            //get encoding
            String encoding = pout.readLine().split("\\s+")[1];
            //System.out.println(encoding + "\t" + fname);
            return encoding.toUpperCase();
			
        } catch (Exception e){
            return "Error with file: "+errorout+" Stacktrace:"+ e.getStackTrace();
        }
    }
    
	
	/*
     * change each line to hg19 if needed
     * chrom and pos are changed
    */
	public static ArrayList<String> liftOver(ArrayList<String> lines, String chainFileName){
		ArrayList<String> out = new ArrayList<String>();
        for(String line: lines) {
			if(line.startsWith("#")){
				out.add(line);
				continue;
			}
			String[] columns= line.split("\t",3);
			String chrom = columns[0];
			//System.out.println(columns[1]);
			int pos = Integer.parseInt(columns[1]);
			ArrayList<String> chain = findChain(chainFileName, chrom, pos);
			String newchrom = chain.get(0).split("\\s+")[7];
			int newpos = searchInChain(chain,pos);
					
			out.add(newchrom + "\t" + newpos + "\t" + columns[2]);
		}
		return out;
	}


    /*
     * for a chrom and pos: finds the relevant chain entry of hg38ToHg19.over.chain
     * returns the entry as an ArrayList<String> of lines
    */
	public static ArrayList<String> findChain(String chainFileName, String chrom, int chrompos){
        ArrayList<String> outchain = new ArrayList<String>();
        try{
			BufferedReader br = new BufferedReader(new FileReader(chainFileName));
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
			br.close();
        }catch(IOException e){
			e.printStackTrace(System.out);
		}
        return outchain;
    }

    /*
     * Is given a chain block and a position
     * finds the hg19 position using the information of the chian block
    */
	public static int searchInChain(ArrayList<String> chain, int chrompos){
        String[] infos = chain.get(0).split("\\s+");
        int tpos = Integer.parseInt(infos[5]);
        int qpos = Integer.parseInt(infos[10]);
        
        //System.out.printf("\n%9s\t%9s\n","tpos","qpos");
        for(int i = 1; i< chain.size() -1; i++){
            //System.out.printf("%9d\t%9d\n",tpos,qpos);
            // {size, dt, dq}
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
        //letzte zeile
		String blockinfo = chain.get(chain.size() -1);
		int size = Integer.parseInt(blockinfo);
		
		tpos += size; // matches
        qpos += size;
            
        if(tpos >= chrompos)
            return qpos - (tpos-chrompos);
            
        return -2; //should not happen
    }
}

