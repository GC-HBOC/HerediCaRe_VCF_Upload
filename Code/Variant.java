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
    
    //add pclass 
    public void addFilteredClass(String s){
        if(!s.equals("")) patho=s;
    }
    

    
    // Displays as string. will be one line in the output file
    public String toString(){
        return ("into VCF_UPLOAD (MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT,GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,CHROM,GPOS,REF,ALT,PATH) values ('" +
                gene+"','"+
                nm+"','"+
                hgvsC+"','"+
                hgvsP+"','"+
                annotation+"','"+
                chrom+"','"+
                posHg19+"','"+
                refHg19+"','"+
                altHg19+"','"+
                patho+"')").replaceAll("&", "|");
    }
	
	
	//checks if this transk is a duplicate
	public boolean isDuplicateIn(List<Variant> trueTranscripts){
		for(Variant t: trueTranscripts){
			if(t.equals(this)){//checks if not a duplicate
				return true;
			}
		}
		return false;
	}
	
	
	//compares gene nm hgvs_c hgvs_p chrom pos ref and alt 
	//does not compare annotation
	public boolean equals(Variant other){
		return (this.gene.equals(other.gene)) &&
			   (this.nm.equals(other.nm)) &&
			   (this.hgvsC.equals(other.hgvsC)) &&
			   (this.hgvsP.equals(other.hgvsP)) &&
			   (this.chrom.equals(other.chrom)) &&
			   (this.posHg19.equals(other.posHg19)) &&
			   (this.refHg19.equals(other.refHg19)) &&
			   (this.altHg19.equals(other.altHg19));
	}
}