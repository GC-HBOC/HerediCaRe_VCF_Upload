import java.util.*;

/*
 * A Variant object can contain all information of a given Variant
 * equals and toString methods are implemented
 * one object will result in one line in the output file
*/
public class Variant{
    
    private String memberId;
    private String bogenNr;
    private String erfmit;
    private String datum;
    private String gene;
    private String nm;
    private String hgvsC;
    private String hgvsP;
    private String annotation;
    private String patho;
    private String chrom;
    private String posHg19;
    private String refHg19;
    private String altHg19;
    private String posHg38;
    private String refHg38;
    private String altHg38;
    private String zygot;
    public Variant(){

    }
    public Variant(String[] infos,String[] cgra){
        nm = infos[6];//NM_...
        addcgra(cgra);
        addInfos(infos);
    }
    public Variant(String[] infos,String filteredTranscript,String[] cgra){
        nm=filteredTranscript; //NM
        addcgra(cgra);
        addInfos(infos);
    }
    
    
    public void addcgra(String[] cgra){
        chrom = cgra[0];
        posHg19 = cgra[1];
        refHg19 = cgra[2];
        altHg19 = cgra[3];
    }
    
    //add gene, hgvs_c, hgvs_p, annotation 
    public void addInfos(String[] infos){
        gene = infos[3];//Gene Name
        for(String s: infos){
            if(s.contains("c.")) hgvsC=s; //HGVS.c
            if(s.contains("p.")) hgvsP=s; //HGVS.p
        }
        annotation = infos[1];//Annotation
    }
    
    //add pclass 
    public void addFilteredClass(String s){
        if(!s.equals("")) patho=s;
    }
    
    //add member_id, bogen_nr, erfmit 
    public void addFilename(String file){
        String[] sfile= file.split("\\.");
        if(sfile.length>0){
            String[] infos=sfile[0].split("\\-");
            if(infos.length == 5){
                memberId=infos[1];
                bogenNr=infos[2];
                erfmit=infos[3];
                datum=infos[4];
            }
        }
    }
    
    // Displays as string. will be one line in the output file
    public String toString(){
        return ("into VCF_UPLOAD (MEMBER_ID,BOGEN_NR,ERFMIT,ERFDAT,GEN2,REFSEQ,HGVS_DNA,HGVS_PROT,ART,CHROM,GPOS,REF,ALT,PATH) values ('" +
                memberId+"','"+
                bogenNr+"','"+
                erfmit+"','"+
                datum+"','"+
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