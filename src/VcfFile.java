import java.util.*;

/*
 * Class to save several informations of a vcf file
 * initialize some values with the filename and get the rest from the content while parsing
 */
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
    ArrayList<String> snpEffVariants;

    public VcfFile(String[] infos){
        this.referenceGenome = infos[0];
        this.memberId = infos[1];
        this.bogenNr = infos[2];
        this.erfmit = infos[3];
        this.datum = infos[4];
        this.metadata = new ArrayList<String>();
        this.variants = new ArrayList<String>();
        this.liftedOver = new ArrayList<String>();
        this.normalized = new ArrayList<String>();
        this.passedVariants = new ArrayList<Variant>();
        this.snpEffVariants = new ArrayList<String>();
    }

}
