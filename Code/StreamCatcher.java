import java.io.*;

/*
 * StreamCatcher is used to handle the snpEff process error output stream. 
 * It can also get a String which is then just returned (used to write different error messages in the output file)
 *
*/
class StreamCatcher extends Thread{
    private InputStream is;
    private String output ="";
    
    StreamCatcher(InputStream is){
        this.is = is;
    }
    
    StreamCatcher(String output){
        this.output = output;
    }
    
    //reads the stream and dumps it in output
    public void run(){
        try{
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line=null;
            while ( (line = br.readLine()) != null){
                output+=line;
            }
        } catch (IOException ioe){
            System.out.println("ERROR while reading Stream");
        }
    }
    
    public String getOutput(){
        return output;
    }
}