import java.io.*;
import java.util.*;
import java.util.ArrayList;


public class ParseVcf {
    String inputFolder;
    String normalizedFolder;
    String processedFolder;
    String errorFolder;
    String uploadFolder;
    String referenceGenomeFasta;
    String panelTranscriptFileName;

    public ParseVcf(){
        this.inputFolder = "vcf_input";
        this.normalizedFolder = "vcf_normalized";
        this.processedFolder = "vcf_processed";
        this.errorFolder = "vcf_error";
        this.uploadFolder = "vcf_upload_files";
        this.referenceGenomeFasta = "ChromFa";
        this.panelTranscriptFileName = "C:\\Users\\Rudel\sFankep\\Documents\\GitHub\\HerediCaRe_VCF_Upload\\BRCA_Panel_Transcripts_v009.tsv";

    }
    
    public static void main(String[] args){
        String JAVA = "java";
        String SNPEFF = "C:\\Users\\Rudel\sFankep\\Documents\\GitHub\\HerediCaRe_VCF_Upload\\snpEff\\snpEff.jar";
        int RAMVALUE = 3;
        String GENOME = "GRCh38.mane.1.0.refseq";

        // ToDo parse arguments later

        ParseVcf myParser = new ParseVcf();
        ParseVcfUtils myUtils = new ParseVcfUtils();
        

        //create directories if not already existing
        String[] neededFolders = {myParser.processedFolder,myParser.normalizedFolder, myParser.errorFolder, myParser.uploadFolder}; 
        for(String i:neededFolders){ myUtils.createDir(i);};
            
        // create a Hashmap for the transcript
        HashMap<String, String> transcriptMap = myUtils.makeTranscriptMap(myParser.panelTranscriptFileName);
        
        //Get all available vcf files
        File vcfInputDir = new File(myParser.inputFolder);
        String[] vcfToParse = vcfInputDir.list();
        for (String file : vcfToParse) {
            if (file.endsWith(".vcf")){
                // parse the filename to get some metadata infos = [reference, patientenID, MGU-BogenNr, MitarbeiterID, Zeitstempel]
                try {
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
                    
                    // normalise the file
                    //String[] returnval_errortext = Normalize.split_lines(VCF_IN_DIR+"\\"+file,getrennt);
                    String[] returnval_errortext = Normalize.split_lines(fileObject.variants);
                    if(!returnval_errortext[0].equals("0")){
                        System.out.println("ERROR while unraveling");
                        //write_output_and_move(file,new ArrayList<Transcript>(),new StreamCatcher(returnval_errortext[1]));
                        continue;
                    }// first split_lines, then normalize und save in a file

                    System.out.println("Normalize into "+myParser.normalizedFolder);
                    fileObject.normalized = new ArrayList<>(fileObject.metadata);
                    fileObject.normalized.addAll(fileObject.variants);
                    if(Normalize.genenotation(fileObject.normalized,myParser.normalizedFolder+"\\"+file,fileObject.referenceGenome) == 1){
                        System.out.println("ERROR while normalizing");
                        //write_output_and_move(file,new ArrayList<Transcript>(),new StreamCatcher("ERROR while normalizing"));
                        continue;
                    }
                        //create a new process
                        System.out.println("Process is running...");
                        System.out.println(JAVA+" -Xmx"+RAMVALUE+"g -jar \""+SNPEFF+"\" "+ GENOME + " -noLog "+myParser.normalizedFolder+"\\"+file);
                        Process p = Runtime.getRuntime().exec(JAVA+" -Xmx"+RAMVALUE+"g -jar \""+SNPEFF+"\" "+ GENOME + " -noLog "+myParser.normalizedFolder+"\\"+file);
                        
                        
                        //List<Transcript>  trueTranscripts = new ArrayList<Transcript>();
                        
                        //catch the Error Stream
                        StreamCatcher errC = new StreamCatcher(p.getErrorStream());
                        errC.start();
                        //Interpreting the Output Stream 
                        BufferedReader pout = new BufferedReader(new InputStreamReader(p.getInputStream()));
                        
                        // for every line in snpEff output
                        for(String s=pout.readLine(); s!=null; s=pout.readLine()){
                            System.out.println("one"+s);

                        }
                        
                        //Rückgabewert des SNPEFF Programmes
                        p.waitFor();
                    

                    // save the lifted over part
                    //fileObject.liftedOver = hgLiftOver.encoding_and_liftover(fileObject.variants, fileObject.referenceGenome);
                    //System.out.println("origin"+ fileObject.variants);
                    //System.out.println("lifted"+ fileObject.liftedOver);
                    //hgLiftOver.l
                }catch (Exception ex) {
                    System.out.println("Process handeling "+file+" failed.");
                    ex.printStackTrace();
                }
                
            }
        //System.out.println(transcriptMap);
            
        }

    }
}