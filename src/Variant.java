import java.util.*;

/*
 * A Variant object can contain all information of a given Variant
 * equals and toString methods are implemented
 * one object will result in one line in the output file
*/
public class Variant{
    String gene;
    String nm;
    String hgvsC;
    String hgvsP;
    String annotation;
    String patho;
    String chrom;
    String posHg19;
    String refHg19;
    String altHg19;
    String posHg38;
    String refHg38;
    String altHg38;
    String zygot;
    String failure;
    
    public Variant(){
        this.failure ="";
    }
    
    //To see major informations in case of debug
    public String toString(){
        return gene+"\t"+nm+"\t"+hgvsC+"\t"+hgvsP+"\t"+annotation+"\t"+patho+"\t"+chrom+"\t"+posHg19+"\t"+refHg19+"\t"+altHg19+"\t"+posHg38+"\t"+refHg38+"\t"+altHg38+"\t"+zygot+"\t"+failure;
    }
    
    // Displays as string. will be one line in the output file
    public String toString(String member, String bogen, String mitId, String datum){
        return ("\tINTO VCF_UPLOAD (MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT,GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,PATH,CHROM,POS_HG19,REF_HG19,ALT_HG19,POS_HG38,REF_HG38,ALT_HG38,ZYGOT) VALUES ('" +
        member+"','"+
        bogen+"','"+
        mitId+"','"+
        datum+"','"+
        gene+"','"+
        nm+"','"+
        hgvsC+"','"+
        hgvsP+"','"+
        annotation+"','"+
        patho+"','"+
        chrom+"','"+
        posHg19+"','"+
        refHg19+"','"+
        altHg19+"','"+
        posHg38+"','"+
        refHg38+"','"+
        altHg38+"','"+
        zygot+"')").replaceAll("&", "|");
    }

	
	
	//checks if this variant is a duplicate
	public boolean isDuplicateIn(ArrayList<Variant> trueTranscripts){
		for(Variant t: trueTranscripts){
			if(t.equals(this)){//checks if not a duplicate
				return true;
			}
		}
		return false;
	}

    /*
     * A check for Mane select and Mane select clinical variants
     * return the position of the duplicate in the original list to compare the two variants
     */
    public int checkClinical(ArrayList<Variant> trueTranscripts){
        int pos=0;
		for(Variant t: trueTranscripts){
			if(t.softEquals(this)){//checks if gene and hgvsC are same
				return pos;
			}
            pos+=1;
		}
		return -1;
	}
	
	
	//compares gene nm hgvs_c hgvs_p chrom pos ref and alt 
	//does not compare annotation
	private boolean equals(Variant other){
		return (this.gene.equals(other.gene)) &&
			   (this.nm.equals(other.nm)) &&
			   (this.hgvsC.equals(other.hgvsC)) &&
			   (this.hgvsP.equals(other.hgvsP)) &&
			   (this.chrom.equals(other.chrom)) &&
			   (this.posHg19.equals(other.posHg19)) &&
			   (this.refHg19.equals(other.refHg19)) &&
			   (this.altHg19.equals(other.altHg19));
	}
    private boolean softEquals(Variant other){
		return (this.gene.equals(other.gene)) &&
			   (this.hgvsC.equals(other.hgvsC));
	}
}