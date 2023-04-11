import java.io.*;
import java.util.*;

/*
 * Utils class for several methods use in the main Parser Class
 * 
 */

public class ParseVcfUtils {

    public void createDir(String dir){
        File myDir = new File(dir);
        if (!myDir.exists()) {
            try{
                myDir.mkdir();
                System.out.println("Created directory: "+dir);
            }
            catch(SecurityException se){
                System.out.println("ERROR while creating "+dir);
            }
        }
    }

    public HashMap<String, String[]> makeTranscriptMap(String transcriptFileName){
        HashMap<String, String[]> transcripts = new HashMap<String, String[]>();
        Scanner transcriptReader = null;
        try{
            File transcriptFile = new File(transcriptFileName);
            transcriptReader = new Scanner(transcriptFile);
            while(transcriptReader.hasNextLine()){
                String transcriptLine = transcriptReader.nextLine().trim();
                if(!transcriptLine.toUpperCase().contains("NM")) continue;
                String[] keyVal = transcriptLine.split("\\s+");
                if(keyVal.length == 3){
                    transcripts.put(keyVal[1].split("\\.")[0], new String[] {keyVal[0],keyVal[2]});
                }
            }
            
        } catch(FileNotFoundException se){
            System.out.println("ERROR while reading "+transcriptFileName);
            System.exit(1);
        } finally{
            if(transcriptReader !=null){
            transcriptReader.close();}
        }
        return transcripts;
    }

    public void saveVcfContent(ArrayList<String> metadata, ArrayList<String> variants, String vcfPath){
        
        Scanner vcfReader = null;
        try{
            File vcfFile = new File(vcfPath);
            vcfReader = new Scanner(vcfFile);
            Boolean var = false;
            while(vcfReader.hasNextLine()){
                String vcfLine = myStrip(vcfReader.nextLine().trim(), "\"");
                //ignore blank lines and variants in mitochondrion
                if (vcfLine.isEmpty()|| vcfLine.startsWith("M")|| vcfLine.startsWith("chrM") ){
                    continue;            
                } else if (vcfLine.startsWith("#CHROM", 0)){
                    metadata.add(vcfLine);
                    var=true;             
                } else if(vcfLine.startsWith("#", 0)){
                    metadata.add(vcfLine);
                } else if (var) {
                    variants.add(fixVariantLine(vcfLine)); 
                } else{
                    metadata.add("##"+vcfLine);
                }
            }
            
        } catch(FileNotFoundException se){
            System.out.println("ERROR while reading "+vcfPath);
        } finally{
            if(vcfReader !=null){
                vcfReader.close();}
        }
    }

    public void snpEffFile(String tmpFile, ArrayList<String> metadata, ArrayList<String> content){
        
        BufferedWriter output=null;
        try {
            File outfile = new File(tmpFile);
            output = new BufferedWriter(new FileWriter(outfile));
            for(String line: metadata){
                output.write(line+"\n");
            }
            for(String line: content){
                output.write(line+"\n");
            }
            
        } catch(Exception e){
            System.out.println("Error while creating snpEff input file "+tmpFile);
        } finally{
            if(output !=null){
                try {
                    output.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }}
        }
    }

    public ArrayList<Variant> parseVariant(String line, HashMap<String,String[]> interest){
        ArrayList<Variant> multSave = new ArrayList<Variant>();
        String[] lineList = line.trim().split("\t");  
        String[] variantInfos = lineList[7].split(";");
        String annotation = null;
        for(String info:variantInfos){
            if(info.contains("ANN")){
                annotation = info.split("=")[1];
                break;
            } // Stop so far the smpEff annotation has been found 
        }

        // If there is no annotation just return a Variant object with the failure message.
        if (annotation==null){
            Variant toSave = new Variant();
            toSave.failure += "#No snpEff Annotation";
            multSave.add(toSave);
            return multSave;
        }
        
        //check if the annotation has a seperator and in this case go through every annotation and save one Variant per annotation
        //  else just save one Variant
        // Also check if the transcript ID is corresponding to the list of interest 
        if (annotation.contains(",")){
            ArrayList<Variant> firstStep = new ArrayList<Variant>();
            for(String ann: annotation.split(",")){
                Variant toSave = new Variant();
                toSave.chrom = fixChr(lineList[0]);
                toSave.posHg38 = lineList[1];
                toSave.refHg38 = lineList[3];
                toSave.altHg38 = lineList[4];
                toSave.patho = filterForPatho(variantInfos);
                if (toSave.patho == null){toSave.failure+="#Invalid CLASS";}
                try{
                    String[] variantFormat = lineList[9].split(":", 0);
                    toSave.zygot = filterForZygot(variantFormat);
                } catch(ArrayIndexOutOfBoundsException e){
                    toSave.zygot="###";
                }
                String[] annList = ann.split("\\|",-1);
                if (interest.containsKey(annList[6].split("\\.")[0])){
                    
                    toSave.gene = annList[3];
                    toSave.nm = annList[6];
                    toSave.hgvsC = annList[9];
                    toSave.hgvsP = annList[10];
                    toSave.annotation = annList[1];
                    int duplicate;
                    try{
                        duplicate = toSave.checkClinical(firstStep);
                    }catch(NullPointerException e){
                        duplicate = -1;
                    }
                    
                    // If the variant is not a duplicate(-1) save it
                    // else check if it is a MANE clinical transcript(1) and save, Mane select is 0.
                    if (duplicate==-1){
                        firstStep.add(toSave);
                    } else{
                        String duplicateTranscript = firstStep.get(duplicate).nm.split("\\.")[0];
                        if(Integer.parseInt(interest.get(toSave.nm.split("\\.")[0])[1]) > Integer.parseInt(interest.get(duplicateTranscript)[1])){
                            firstStep.set(duplicate, toSave);
                        }
                    }
                }else{
                    toSave.failure +="#Invalid TranscriptID";
                    firstStep.add(toSave);
                }            
            }
            //Another line filter. LATER
            for (Variant i:firstStep){
                if (i.failure.length()==0){
                    multSave.add(i);

                }
            }
            if (multSave.isEmpty()){multSave.add(firstStep.get(0));}



        } else{
            String[] ann = annotation.split("\\|",-1);
            Variant toSave = new Variant();
            toSave.chrom = fixChr(lineList[0]);
            toSave.posHg38 = lineList[1];
            toSave.refHg38 = lineList[3];
            toSave.altHg38 = lineList[4];
            toSave.patho = filterForPatho(variantInfos);
            if (toSave.patho == null){toSave.failure+="#Invalid CLASS";}
            try{
                String[] variantFormat = lineList[9].split(":", 0);
                toSave.zygot = filterForZygot(variantFormat);
            } catch(ArrayIndexOutOfBoundsException e){
                toSave.zygot="###";
            }
            if (interest.containsKey(ann[6].split("\\.")[0])){
                toSave.gene = ann[3];
                toSave.nm = ann[6];
                toSave.hgvsC = ann[9];
                toSave.hgvsP = ann[10];
                toSave.annotation = ann[1];
            }else{toSave.failure +="#Invalid TranscriptID";}
            multSave.add(toSave);
        }
        return multSave;
    }

    //Method to get the class of a variant
    private static String filterForPatho(String[] infos){
        String out="";
        for(String n: infos){
            if(n.contains("MutDB:Classification")){
                String c=n.split("=")[1];
                if(c.equals("pathogenic")) out = "5";
                if(c.equals("likelypathogenic")) out = "4";
                if(c.equals("uncertainsignificance")) out = "3";
                if(c.equals("likelybenign")) out = "2";
                if(c.equals("benign")) out = "1";
                if(c.equals("undefined")) out = "";
                return out;
            }
        }
        for(String n: infos){
            if(n.contains("CLASS")){
                return n.split("=")[1];
            }
        }
        for(String n: infos){
            if(n.contains("MT")){
                String c=n.split("=")[1];
                if(c.equals("PAT")) out = "5";
                if(c.equals("LPAT")) out = "4";
                if(c.equals("UI") || c.equals("UV")) out = "3";
                if(c.equals("LBEN")) out = "2";
                if(c.equals("BEN")) out = "1";
                if(c.equals("UD")) out = "";
                if(c.equals("A") || c.equals("FR")) out = null;
                return out;
            }
        }
        return out;
    }

    //Method to get the genotype of a variant
    private static String filterForZygot(String[] format){
        switch(format[0]){
            case "0/0":
                return "0/0";
            case "0/1":
                return "0/1";
            case "1/0":
                return "0/1";
            case "1/1":
                return "1/1";
            default:
                return "###";
        }

        
    }

    //Just to fix the notation of the chromosome column: from "chr1" to "1" or "chrM" to "M" 
    private static String fixChr(String chrom){
        
        if (chrom.contains("chr")){
            return chrom.substring(3);
        } else{ return chrom;}
    }

    // check every column of a variant line for a space character and delete + remove '"'
    private static String fixVariantLine(String line){
        String endLine = "";
        for(String entry:line.split("\t")){
            endLine+=myStrip(entry, "\"").replaceAll("\s", "")+"\t";
        }
        return endLine.trim();
   
    }

    //private strip method to delete chars at the beginning und at the end of a string
    private static String myStrip(String line, String character){
        
        if(line.startsWith(character) && line.endsWith(character)){
            return line.substring(1, line.length()-1);
        } else if(line.startsWith(character)){
            return line.substring(1);
        } else if(line.endsWith(character)){
            return line.substring(0, line.length()-1);
        } else{
            return line;
        }
    }

    /*
     * write the true variants in the output file
     */
    public void writeOutput(String member, String bogen, String mitId, String datum, ArrayList<Variant> toWrite, String filename){
        BufferedWriter output=null;
        try {
            File outfile = new File(filename);
            output = new BufferedWriter(new FileWriter(outfile));
            output.write("INSERT ALL");
            output.newLine();
            for(Variant each: toWrite){
                output.write(each.toString(member, bogen, mitId, datum));
                output.newLine();
            }
            output.write("SELECT * FROM dual;");
            output.newLine();
            output.write("COMMIT;");
            output.newLine();
        } catch(Exception e){
            System.out.println("Error while handling file "+filename);
        } finally{
            if(output !=null){
                try {
                    output.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }}
        }
    }

    /*
     * method to remove duplicate variants in a list (Now not in use)
     */
    public ArrayList<Variant> cleanVariants(ArrayList<Variant> variants, HashMap<String,String[]> interest){
        ArrayList<Variant> cleaned = new ArrayList<Variant>();
        for(Variant entry:variants){
            entry.isDuplicateIn(cleaned);
        }
        return cleaned;
    }

    /*
     * copy a file to a new destination and delete the original
     */
    public void moveFile(String source, String destination) {
        File mySource = new File(source);
        if (mySource.renameTo(new File(destination))){
            mySource.delete();
            
        } else{System.out.println("Error while moving file "+source);}
    }
}