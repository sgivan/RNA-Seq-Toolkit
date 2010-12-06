#!/bin/bash
# auxilary script of the RNAseq package
# primarily used by SAG in FGMG project
#

echo "primer_trim.pl --infile set1.fq --fastq --idlist --outfile set1_short.fq --primerseq AAGCAGTGGTATCAACGCAGAGTACATGGG --overwrite --idlistfile set1_short_idlist.txt"
primer_trim.pl --infile set1.fq --fastq --idlist --outfile set1_short.fq --primerseq AAGCAGTGGTATCAACGCAGAGTACATGGG --overwrite --idlistfile set1_short_idlist.txt 
    
echo "primer_trim.pl --infile set2.fq --fastq --idlist --outfile set2_short.fq --primerseq AAGCAGTGGTATCAACGCAGAGTACATGGG --overwrite --idlistfile set2_short_idlist.txt"
primer_trim.pl --infile set2.fq --fastq --idlist --outfile set2_short.fq --primerseq AAGCAGTGGTATCAACGCAGAGTACATGGG --overwrite --idlistfile set2_short_idlist.txt
    
echo "primer_trim.pl --infile set1.fq --fastq --idlist --outfile set1_medium.fq --primerseq AAGCAGTGGTATCAACGCAGAGTAC --overwrite --idlistfile set1_medium_idlist.txt"
primer_trim.pl --infile set1.fq --fastq --idlist --outfile set1_medium.fq --primerseq AAGCAGTGGTATCAACGCAGAGTAC --overwrite --idlistfile set1_medium_idlist.txt 
    
echo "primer_trim.pl --infile set2.fq --fastq --idlist --outfile set2_medium.fq --primerseq AAGCAGTGGTATCAACGCAGAGTAC --overwrite --idlistfile set2_medium_idlist.txt"
primer_trim.pl --infile set2.fq --fastq --idlist --outfile set2_medium.fq --primerseq AAGCAGTGGTATCAACGCAGAGTAC --overwrite --idlistfile set2_medium_idlist.txt 
     
echo "notseq.pl --idfile set1_short_idlist.txt --idcol 1 --seqfile set1_medium.fq --outfile set1_medium_noprimer.fq --fastq"
notseq.pl --idfile set1_short_idlist.txt --idcol 1 --seqfile set1_medium.fq --outfile set1_medium_noprimer.fq --fastq 
     
echo "notseq.pl --idfile set2_short_idlist.txt --idcol 1 --seqfile set2_medium.fq --outfile set2_medium_noprimer.fq --fastq"
notseq.pl --idfile set2_short_idlist.txt --idcol 1 --seqfile set2_medium.fq --outfile set2_medium_noprimer.fq --fastq
    
mv -f set1_medium_noprimer.fq set1_medium.fq
mv -f set2_medium_noprimer.fq set2_medium.fq
     
echo "cat set1_medium_idlist.txt set1_short_idlist.txt | sort | uniq > set1_idlist_uniq.txt"
cat set1_medium_idlist.txt set1_short_idlist.txt | sort | uniq > set1_idlist_uniq.txt 
  
echo "cat set2_medium_idlist.txt set2_short_idlist.txt | sort | uniq > set2_idlist_uniq.txt"
cat set2_medium_idlist.txt set2_short_idlist.txt | sort | uniq > set2_idlist_uniq.txt
    
echo "notseq.pl --idfile set1_idlist_uniq.txt --seqfile set1.fq --fastq --outfile set1_80mers.fq --idcol 1"
notseq.pl --idfile set1_idlist_uniq.txt --seqfile set1.fq --fastq --outfile set1_80mers.fq --idcol 1 
    
echo "notseq.pl --idfile set2_idlist_uniq.txt --seqfile set2.fq --fastq --outfile set2_80mers.fq --idcol 1"
notseq.pl --idfile set2_idlist_uniq.txt --seqfile set2.fq --fastq --outfile set2_80mers.fq --idcol 1

touch set1_trimmed.fq set2_trimmed.fq
cat set1_short.fq set1_medium.fq set1_80mers.fq >> set1_trimmed.fq
cat set2_short.fq set2_medium.fq set2_80mers.fq >> set2_trimmed.fq
rm set1_short.fq set1_medium.fq set1_80mers.fq set2_short.fq set2_medium.fq set2_80mers.fq
rm set1_idlist_uniq.txt set2_idlist_uniq.txt
ln -sf set1_trimmed.fq set1.fq
ln -sf set2_trimmed.fq set2.fq

