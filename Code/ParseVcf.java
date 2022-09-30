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
        

        //create directory if not already existing
        String[] neededFolders = {myParser.processedFolder, myParser.errorFolder, myParser.uploadFolder}; 
        for(String i:neededFolders){ myUtils.createDir(i);};

        // create a Hashmap for the transcript
        HashMap<String, String> transcriptMap = myUtils.makeTranscriptMap(myParser.panelTranscriptFileName);
        
        //Get all available vcf files
        File vcfInputDir = new File(myParser.inputFolder);
        String[] vcfToParse = vcfInputDir.list();
        for (String string : vcfToParse) {
            System.out.println(string);
            
        }

    }
}