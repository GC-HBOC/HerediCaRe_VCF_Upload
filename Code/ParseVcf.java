import java.io.*;
import java.util.*;
import java.util.ArrayList;


public class ParseVcf {
    String inputFolder;
    String processedFolder;
    String errorFolder;
    String uploadFolder;
    String referenceGenomeFasta;
    String panelTranscriptFileName;

    public ParseVcf(){
        this.inputFolder = "vcf_input";
        this.processedFolder = "vcf_processed";
        this.errorFolder = "vcf_error";
        this.uploadFolder = "vcf_upload_files";
        this.referenceGenomeFasta = "ChromFa";
        this.panelTranscriptFileName = "C:\\Users\\Rudel\sFankep\\Documents\\GitHub\\HerediCaRe_VCF_Upload\\BRCA_Panel_Transcripts_v009.tsv";

    }
    
    public static void main(String[] args){

        // ToDo parse arguments later

        ParseVcf myParser = new ParseVcf();
        ParseVcfUtils myUtils = new ParseVcfUtils();
        

        //create directories if not already existing
        String[] neededFolders = {myParser.processedFolder, myParser.errorFolder, myParser.uploadFolder}; 
        for(String i:neededFolders){ myUtils.createDir(i);};

        // create a Hashmap for the transcript
        HashMap<String, String> transcriptMap = myUtils.makeTranscriptMap(myParser.panelTranscriptFileName);
        
        //Get all available vcf files
        File vcfInputDir = new File(myParser.inputFolder);
        String[] vcfToParse = vcfInputDir.list();
        for (String file : vcfToParse) {
            if (file.endsWith(".vcf")){
                // parse the filename to get some metadata infos = [reference, patientenID, MGU-BogenNr, MitarbeiterID, Zeitstempel]
                String errorOutput="ok";
                int i = file.lastIndexOf(".");
                String[] sfile =  {file.substring(0, i), file.substring(i)};
                String[] infos=new String[5];
                if(sfile.length>0){
                    String[] temp =sfile[0].split("\\-");
                    if(temp.length == 5 ) infos = temp;
                    else{ errorOutput= "Fehlerhafter Dateiname";}
                }else{ errorOutput= "Fehlerhafter Dateiname";}

                // create a file object
                VcfFile fileObject = new VcfFile(infos);

                // save the header and the Variants from vcf file
                myUtils.saveVcfContent(fileObject, vcfInputDir + "\\" + file); 
                
                // save the lifted over part
                fileObject.liftedOver = hgLiftOver.encoding_and_liftover(fileObject.variants, fileObject.referenceGenome);
                System.out.println("origin"+ fileObject.variants);
                System.out.println("lifted"+ fileObject.liftedOver);
                hgLiftOver.l
            }
        //System.out.println(transcriptMap);
            
        }

    }
}