## Fast-tracking bespoke DNA reference database generation from museum collections for biomonitoring and conservation

Andrew Dopheide, Talia Brav-Cubitt, Anastasija Podolyan, Richard A. B. Leschen, Darren Ward, Thomas R. Buckley, and Manpreet K. Dhami,

Data processing used to obtain invertebrate DNA barcodes from Illumina MiSeq output. These processes require [Illumina bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html), [Claident](https://github.com/astanabe/Claident/), [VSEARCH](https://github.com/torognes/vsearch), and [EMBOSS](http://emboss.sourceforge.net/). 

### Generation of sequence data

To generate the MiSeq data, two overlapping fragments (FC and BR) of the mitochodrial COI gene were amplified by PCR from each of 450 invertebrate specimens, using the primer pairs Ill_LCO1490 with Ill_FC-R, and Ill_BR-F with Ill_HCO2198, based on [Shokralla et al. 2015, Massively parallel multiplex DNA sequencing for specimen identification using an Illumina MiSeq platform. Scientific Reports 5:9687, DOI: 10.1038/srep09687](http://www.nature.com/srep/2015/150408/srep09687/full/srep09687.html). The amplicons were tagged with Illumina Nextera XT forward or reverse adapters and padded with 0-4 mer nucleotide spacers to increase sequence heterogeneity. Unique 8-mer multiplex tags were added to the amplicons from each specimen in a second PCR step. The amplicons were processed on an Illumina MiSeq system in 2 x 300 b.p. paired-end sequencing mode.  

### Data processing steps

1. Raw MiSeq data in bcl format was converted into fastq format, without demultiplexing, using `bcl2fastq`. 

2. The forward and reverse fastq sequences were each demultiplexed into separate FC and BR fastq files for each of 450 invertebrate specimens, using `clsplitseq` from Claident. 
For example, to demultiplex/extract FC amplicon sequences from the raw fastq data:
``` 
clsplitseq \
	--runname=COI \
	--truncateN=enable \
    --index1file=barcodes_i7_RC.txt # Multiplex indices
	--index2file=barcodes_i5.txt \
    --primerfile=LCO1490 \
	--reverseprimerfile=FC-R \ 
	--minqualtag=20 --numthreads=8 \
    Undetermined_S0_L001_R1_001.fastq.gz \ # Un-demultiplexed sequence files
	Undetermined_S0_L001_I1_001.fastq.gz \
    Undetermined_S0_L001_I2_001.fastq.gz \
	Undetermined_S0_L001_R2_001.fastq.gz \
    Demultiplexed_LCO1490-FC # Output folder for demultiplexed sequences
```
This resulted in a pair of R1 and R2 fastq files for the FC amplicon, and a pair of R1 and R2 fastq files for the BR amplicon, for each of the 450 specimens. 

Steps 3-7 were carried out using code in `merge_filter_seqs.sh`

3. For each amplicon (FC and BR) and specimen, the R1 and R2 sequences contained in each pair of fastq files were merged and error-filtered using `-fastq_mergepairs` in VSEARCH.

4. The filtered sequences from each amplicon and specimen were then dereplicated using `-derep_fulllength`, denoised using `-cluster_unoise` with `-minsize 1`, and filtered for chimeras using `-uchime3_denovo`, all in VSEARCH.

5. A taxonomic identity for each denoised FC and BR sequence was obtained by BLAST against GenBank nr, accepting the top match in each case.

6. Up to n = 20 most abundant denoised FC and BR sequences (if any) were output as fasta files for each specimen, using the VSEARCH command `-derep_fulllength` with option `-topn n`. The FC and BR sequence files were concatenated together into a single fasta file per specimen.  

7. All possible pairwise sequence alignments with the expected overlap between FC and BR amplicons (overlap of 85 bp with 100 % pairwise sequence identity) were identified between the sequences in the concatenated FC and BR fasta file for each specimen using the VSEARCH command `-allpairs_global` with options -id 1 and -mincols 85. 

8. If no alignments were detected among the top 20 sequences from a specimen, steps 6 and 7 were repeated with all denoised FC and BR sequences from that specimen.  

Steps 9-10 were carried out in python using `sequences_to_barcodes.py`. This requires a table of specimen taxonomic identification details, as well as the NBCI taxonomy database taxdmp file `rankedlineage.dmp` (available from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) to parse BLAST identification results.  

9. Putatively correct full length or partial barcodes were identified by comparing BLAST taxonomic identities for FC and BR sequences, and aligned FC-BR sequence pairs, to _a_ _priori_ specimen identification data, prioritising the lowest taxonomic rank level of identifications. The correct FC sequence, BR sequence, or aligned FC-BR sequence pair for each specimen was deemed to be whichever of these had the lowest taxonomic rank identification matching the _a_ _priori_ taxonomic identity of the specimen. The sequences were also required to have lengths between 324 and 326 b.p. for FC and between 417 and 419 b.p. for BR, and BLAST identification bitscores ≥ 200 for FC and ≥ 250 for BR, to avoid spurious matches. If more than one FC sequence, BR sequence, or aligned FC-BR sequence pair per specimen met these criteria, the most abundant sequence or sequence pair was selected.

10. If the selected barcode for a specimen was an aligned FC-BR sequence pair, these sequences were  combined by EMBOSS merger. This carries out a pairwise alignment according to the Needleman-Wunsch algorithm, outputting a merged barcode sequence, with an alignment summary that was checked to ensure the alignment attributes were as expected (overlap of 85 bp and score of 425). Otherwise, if an FC or BR sequence was selected, the FC or BR sequence was output as a partial barcode in fasta format.
 
[![manaakiwhenua-standards](https://github.com/manaakiwhenua/<Fast_DNA_barcoding>/workflows/manaakiwhenua-standards/badge.svg)](https://github.com/manaakiwhenua/manaakiwhenua-standards)
