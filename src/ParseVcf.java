import java.io.*;
import java.util.*;

import net.sourceforge.argparse4j.*;
import net.sourceforge.argparse4j.inf.*;

public class ParseVcf {
    String inputFolder;
    String debugFolder;
    String normalizedFolder;
    String processedFolder;
    String errorFolder;
    String tmpFolder;
    String uploadFolder;
    String panelTranscriptFileName;
    String rejectedVariants;

    public ParseVcf(){
        this.inputFolder = "";
        this.debugFolder = "Debug";
        this.uploadFolder = "Output";
        this.normalizedFolder = this.debugFolder+"\\Normalized\sVCFs";
        this.processedFolder = this.debugFolder+"\\Processed";
        this.errorFolder = this.debugFolder+"\\Error";
        this.tmpFolder = this.debugFolder+"\\.temp";
        this.rejectedVariants = this.debugFolder+"\\Rejected\sVariants";
        this.panelTranscriptFileName = "resources\\BRCA_Panel_Transcripts_v009.tsv";

    }
    
    public static void main(String[] args) throws ArgumentParserException{

        /*
         * initialize some needed values and fill it later with user arguments
         */
        String JAVA = "java";
        String SNPEFF = "resources\\snpEff\\snpEff.jar";
        int RAMVALUE = 3;
        String GENOME = "GRCh38.mane.1.0.refseq";

        ParseVcf myParser = new ParseVcf();
        ParseVcfUtils myUtils = new ParseVcfUtils();

        /*
         * Create arguments with the argparse4j package
         */
        ArgumentParser parser = ArgumentParsers.newFor("ParseVCF").build();

        parser.addArgument("-i")
            .required(true)
            .metavar("<input path>")
            .help("Folder with VCF files to parse (default: "+ myParser.inputFolder +")");
        
        parser.addArgument("-o")
            .setDefault(myParser.uploadFolder)
            .metavar("<output path>")
            .help("Output folder for final .txt files (default: "+ myParser.uploadFolder +")");
        
        parser.addArgument("-sp")
            .setDefault(SNPEFF)
            .metavar("<snpEff path>")
            .help("Path to snpEff JAR file (default: "+SNPEFF+")");

        parser.addArgument("-genome")
            .setDefault(GENOME)
            .metavar("<database>")		
            .help("snpEff database name (default: "+GENOME+")");

        parser.addArgument("-jp")
            .setDefault(JAVA)
            .metavar("<java path>")
            .help("Path to the java executable. Note that Java v>=12 is required to run snpEff! (default: "+ JAVA +")");

        parser.addArgument("-d")
            .setDefault(myParser.debugFolder)
            .metavar("<debug path>")
            .help("Debug Folder. Contains processed, erroneous & normalized VCFs + Rejected Variants(.tsv) (default: "+ myParser.debugFolder +")");

        parser.addArgument("-ram")	
            .setDefault(RAMVALUE)
            .metavar("<value>")
            .help("Accessible RAM (GB) for java virtual machine(default: "+ RAMVALUE +")");

        parser.addArgument("-t")	
            .setDefault(myParser.panelTranscriptFileName)
            .metavar("<transcript>")
            .help("Transcript file (default: "+ myParser.panelTranscriptFileName +")");



        /*
         * Get the parser values, either filled from the user or from the the default values
         */
        try { 
            Namespace res = parser.parseArgs(args);
            myParser.inputFolder = res.get("i");
            myParser.uploadFolder = res.get("o");
            myParser.debugFolder = res.get("d");
            myParser.panelTranscriptFileName = res.get("t");
            SNPEFF = res.get("sp");
            JAVA = res.get("jp");
            GENOME = res.get("genome");
            RAMVALUE = res.get("ram");
        } catch (ArgumentParserException e) {
            parser.handleError(e);System.exit(1);
        }

        /*
         * Check the ressource folder
        */
        if(!(new File("resources").exists()) || ((new File("resources")).list().length ==0) ){System.err.println("Resources folder missing/empty. Check the README and follow the instructions"); System.exit(1);}
        
        File vcfInputDir = new File(myParser.inputFolder);
        String[] vcfToParse = vcfInputDir.list();
        if(!(vcfInputDir.exists()) || (vcfToParse.length ==0) ){System.err.println(myParser.inputFolder + " folder missing/empty"); System.exit(1);}
        /*
         * create a Hashmap for the transcript
         */
        System.out.println("Parsing transcript file " + myParser.panelTranscriptFileName);
        HashMap<String, String[]> transcriptMap = myUtils.makeTranscriptMap(myParser.panelTranscriptFileName);
        
        /*
         * create directories if not already existing
         */
        System.out.println("\nCreating mandatory directories, if not existing");

        String[] neededFolders = {myParser.debugFolder, myParser.processedFolder,myParser.normalizedFolder, myParser.errorFolder, myParser.uploadFolder, myParser.tmpFolder,myParser.rejectedVariants}; 
        for(String i:neededFolders){ myUtils.createDir(i);};
        
        /*
         * Get all available files and work only with files ending by .vcf
         */
        
        
        System.out.println("Found a total of "+ vcfToParse.length +" files to process in the folder: "+myParser.inputFolder);
        for (String file : vcfToParse) {
            if (file.endsWith(".vcf")){
                System.out.println("\n#################### " + file + " processing... ###########################");

                // parse the filename to get some metadata infos = [reference, patientenID, MGU-BogenNr, MitarbeiterID, Zeitstempel]
                try {
                    
                    int i = file.lastIndexOf(".");
                    String[] sfile =  {file.substring(0, i).toLowerCase(), file.substring(i)};
                    String[] infos=new String[5];
                    if(sfile.length>0){
                        String[] temp =sfile[0].split("\\-");
                        if(temp.length == 5 ) infos = temp;
                        else{ throw new Exception("Filename not corresponding to the requirements: 'reference-patId-mguNr-emplID-timeStamp.vcf'");}
                    }else{throw new Exception("Filename not corresponding to the requirements: 'reference-patId-mguNr-emplID-timeStamp.vcf'");}

                    // create a file object passing the filename informations
                    VcfFile fileObject = new VcfFile(infos);

                    // save the header and the Variants from vcf file
                    myUtils.saveVcfContent(fileObject.metadata, fileObject.variants, vcfInputDir + "\\" + file); 
                    
                    /*
                     * normalize the file
                     *   first split_lines 
                     *   then normalize und save in a file
                     */
                    System.out.println("Normalizing...");
                    String[] returnval_errortext = Normalize.split_lines(fileObject.variants, fileObject.normalized);
                    if(!returnval_errortext[0].equals("0")){
                        System.out.println("ERROR while unraveling");
                        continue;
                    }
                    if(Normalize.genenotation(fileObject.metadata, fileObject.normalized, myParser.normalizedFolder+"\\"+file, fileObject.referenceGenome) == 1){
                        System.out.println("ERROR while normalizing");
                        continue;
                    }

                    /*
                     * LiftOver Variants depending of the reference genome of the original input file
                     *  Create a temporary file, which will be deleted after snpEff running
                     */
                    fileObject.liftedOver = hgLiftOver.encoding_and_liftover(fileObject.variants, fileObject.referenceGenome);
                    
                    if (fileObject.referenceGenome.equals("hg19")){
                        myUtils.snpEffFile(myParser.tmpFolder+"\\"+file, fileObject.metadata, fileObject.liftedOver);

                    }else{
                        myUtils.snpEffFile(myParser.tmpFolder+"\\"+file, fileObject.metadata, fileObject.normalized);
                    }

                    //create a new process to run snpEff
                    System.out.println("snpEff is running...");
                    //System.out.println(JAVA+" -Xmx"+RAMVALUE+"g -jar \""+SNPEFF+"\" "+ GENOME +" -noStats -noLog "+myParser.tmpFolder+"\\"+file);
                    Process p = Runtime.getRuntime().exec(JAVA+" -Xmx"+RAMVALUE+"g -jar \""+SNPEFF+"\" "+ GENOME +" -noStats -noLog "+myParser.tmpFolder+"\\"+file);
                    
                    
                    
                    
                    //catch the Error Stream
                    StreamCatcher errC = new StreamCatcher(p.getErrorStream());
                    errC.start();
                    //Interpreting the Output Stream 
                    BufferedReader pout = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    
                    //important to retrieve informations on every line.
                    int pos = 0;

                    /*
                     * for every line in snpEff output
                     *  create a file for rejected variants
                     *  Save variants with no "failure" property
                     */

                    BufferedWriter output=null;
                    File outfile = new File(myParser.rejectedVariants+"\\"+file.replaceFirst(".vcf", ".tsv"));
                    output = new BufferedWriter(new FileWriter(outfile));
                    output.write("CHROM("+fileObject.referenceGenome+")\tPOS\tREF\tALT\tFAILURE\n");
                    for(String s=pout.readLine(); s!=null; s=pout.readLine()){
                        if (!s.startsWith("#")){
                            
                            ArrayList<Variant> allVar = myUtils.parseVariant(s,transcriptMap);
                            for(Variant tmpVar:allVar){                           
                                if (tmpVar.failure.length()==0){
                                    
                                    if (fileObject.referenceGenome.equals("hg19")){
                                        String[] moreInfo = fileObject.normalized.get(pos).trim().split("\t"); 
                                        tmpVar.posHg19 = moreInfo[1];
                                        tmpVar.refHg19 = moreInfo[3];
                                        tmpVar.altHg19 = moreInfo[4];
                                    } else{
                                        String[] moreInfo = fileObject.liftedOver.get(pos).trim().split("\t"); 
                                        tmpVar.posHg19 = moreInfo[1];
                                        tmpVar.refHg19 = moreInfo[3];
                                        tmpVar.altHg19 = moreInfo[4];
                                    }
                                    fileObject.passedVariants.add(tmpVar);

                                } else {
                                    String[] variant = fileObject.normalized.get(pos).trim().split("\t");
                                    output.write(variant[0]+"\t"+variant[1]+"\t"+variant[3]+"\t"+variant[4]+"\t"+tmpVar.failure+"\n");
                                } 
                                
                            }
                            pos+=1;
            
                        }
                    }
                    output.close();
                    
                    //output value from snpEff program
                    p.waitFor();

                    //remove tempfile from snpEff
                    File tmp = new File(myParser.tmpFolder+"\\"+file);
                    tmp.delete();
                    
                    // write the file to upload in the database                  
                    myUtils.writeOutput(fileObject.memberId, fileObject.bogenNr, fileObject.erfmit, fileObject.datum, fileObject.passedVariants, myParser.uploadFolder+"\\"+file.replaceFirst(".vcf", ".txt"));
                    
                    // move the file to the processed folder 
                    myUtils.moveFile(myParser.inputFolder+"\\"+file, myParser.processedFolder+"\\"+file);
                    System.out.println("Done ");
                } catch (Exception ex) {

                    // move the file to the error folder
                    myUtils.moveFile(myParser.inputFolder+"\\"+file, myParser.errorFolder+"\\"+file);
                    System.out.println("Process handeling "+file+" failed.");
                    ex.printStackTrace();
                }

                
            } else{System.out.println("\n"+file + " must end with '.vcf' to be processed. SKIPPED");}
        }

        //delete the tmp folder
        File tmp = new File(myParser.tmpFolder);
        tmp.delete();
        System.out.println("\nProcessed original VCFs in "+myParser.debugFolder+"\\"+myParser.processedFolder);
        System.out.println("Failed original VCFs in "+myParser.debugFolder+"\\"+myParser.errorFolder);

    }//end main
    
}//end class