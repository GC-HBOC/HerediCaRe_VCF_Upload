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
    public HashMap makeTranscriptMap(String transcriptFileName){
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
}