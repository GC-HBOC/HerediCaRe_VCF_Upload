import java.io.*;
import java.util.*;
//import java.lang.*;

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
    public HashMap<String, String> makeTranscriptMap(String transcriptFileName){
        HashMap<String, String> transcripts = new HashMap<String, String>();
        try{
            File transcriptFile = new File(transcriptFileName);
            Scanner transcriptReader = new Scanner(transcriptFile);
            while(transcriptReader.hasNextLine()){
                String transcriptLine = transcriptReader.nextLine().trim();
                if(!transcriptLine.toUpperCase().contains("NM")) continue;
                String[] keyVal = transcriptLine.split("\\s+");
                if(keyVal.length == 2){
                    transcripts.put(keyVal[0], keyVal[1].split("\\.")[0]);
                }
            }
            transcriptReader.close();
        } catch(FileNotFoundException se){
            System.out.println("ERROR while reading "+transcriptFileName);
        }
        return transcripts;
    }

    public void saveVcfContent(VcfFile theFile, String vcfPath){
        ArrayList<String> vcfMetadata = new ArrayList<String>();
        ArrayList<String> vcfVariants = new ArrayList<String>();
        
        try{
            File vcfFile = new File(vcfPath);
            Scanner vcfReader = new Scanner(vcfFile);
            while(vcfReader.hasNextLine()){
                String vcfLine = vcfReader.nextLine().trim();
                if(vcfLine.startsWith("#")){
                    vcfMetadata.add(vcfLine);
                } else {
                    vcfVariants.add(vcfLine);
                }
            }
            vcfReader.close();
        } catch(FileNotFoundException se){
            System.out.println("ERROR while reading "+vcfPath);
        }
        theFile.metadata = vcfMetadata;
        theFile.variants = vcfVariants;
    }
    public void parseVariant(String line, HashMap<String, String> transcript, VcfFile file){
        Variant toSave = new Variant();
        
        if (!line.startsWith("#") && line.contains("ANN")){
            String[] lineList = line.trim().split("\t");
            toSave.chrom = lineList[0];
            toSave.posHg38 = lineList[1];
            toSave.refHg38 = lineList[3];
            toSave.altHg38 = lineList[4];
            String variantInfos = lineList[7];
            String variantFormat = lineList[9];
            
            System.out.println(lineList[0]);
            
        }
    }
}