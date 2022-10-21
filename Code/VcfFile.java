import java.util.*;

public class VcfFile {
    String referenceGenome;
    String memberId;
    String bogenNr;
    String erfmit;
    String datum;
    ArrayList<String> metadata;
    ArrayList<String> variants;
    ArrayList<String> liftedOver;
    ArrayList<String> normalized;
    ArrayList<Variant> passedVariants;

    public VcfFile(String[] infos){
        this.referenceGenome = infos[0];
        this.memberId = infos[1];
        this.bogenNr = infos[2];
        this.erfmit = infos[3];
        this.datum = infos[4];

    }

}
