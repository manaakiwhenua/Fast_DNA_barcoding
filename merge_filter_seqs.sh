#!/bin/bash

module load VSEARCH/2.14.2-GCC-9.2.0

# Merge forward/reverse sequences with VSEARCH, filter for errors and length, and denoise.

nSeqs=20 # Set number of sequences to examine for alignments
amplicons=("LCO1490-FC" "BR-HCO2198")

for amp in ${amplicons[@]}; do
  echo "Processing amplicon '$amp'"
  files=$( ls Demultiplexed_$amp/*.forward.fastq.gz )
  mkdir Processing_$amp # Make folder for results
  for f in $files; do
    echo "Starting case '$f'"
    #continue
    r=$(sed -e "s/forward/reverse/" <<< "$f") # Identifies the corresponding R2 file
    s=$(cut -d_ -f4 <<< "$f") # Identifies sample ID from filename
    echo $s
    if [[ $amp == "LCO1490-FC" ]]
    then
      vsearch -fastq_mergepairs $f -reverse $r -fastq_minovlen 200 -fastq_minmergelen 250 -fastq_maxmergelen 400 -fastq_maxee 1.0 -fastq_maxns 0 -fastaout Processing_$amp/$s.filtered.fasta # Merge and filter LCO1490-FC reads. Expected length is 325.
    elif [[ $amp == "BR-HCO2198" ]]
    then
      vsearch -fastq_mergepairs $f -reverse $r -fastq_minovlen 100 -fastq_minmergelen 350 -fastq_maxmergelen 500 -fastq_maxee 1.0 -fastq_maxns 0 -fastaout Processing_$amp/$s.filtered.fasta # Merge and filter BR-HCO2198 reads. Expected length is 418.
    fi
    f_count=$(grep -c '>' Processing_$amp/$s.filtered.fasta) # Check how many sequences remain after assembly and filtering
    echo $s: $f_count filtered sequences
    if [[ $f_count -gt 0 ]] 
      then
      vsearch -derep_fulllength Processing_$amp/$s.filtered.fasta -relabel $s.Uniq -sizeout -output Processing_$amp/$s.uniques.fasta
      u_count=$(grep -c '>' Processing_$amp/$s.uniques.fasta) # Check how many sequences remain after dereplication
      echo $s: $u_count unique sequences
      echo $s: denoising...
      vsearch -cluster_unoise Processing_$amp/$s.uniques.fasta -relabel $s.denoise -minsize 1 -sizein -sizeout -centroids Processing_$amp/$s.denoise.fasta
      vsearch -uchime3_denovo Processing_$amp/$s.denoise.fasta -nonchimeras Processing_$amp/$s.denoise.good.fasta # Check for chimeras
      d_count=$(grep -c '>' Processing_$amp/$s.denoise.good.fasta) # Check how many sequences remain after denoising
      echo $s: $d_count denoised sequences
      # vsearch -derep_fulllength Processing_$amp/$s.filtered.fasta -relabel $s.Uniq -sizeout -topn 1 -output Processing_$amp/$s.top1.fasta # Output the single most abundant sequence
      if [[ $d_count -gt $nSeqs ]]
      then 
        d_count=$nSeqs; 
      fi # Limit the number of output sequences
      vsearch -derep_fulllength Processing_$amp/$s.denoise.good.fasta -sizein -sizeout -topn $d_count -output Processing_$amp/$s.denoise.topn.fasta # Output n most abundant sequences
     else
      echo $s: No sequences to denoise; starting next...
    fi
    #echo $f: $? >> exitcodes.txt # Not sure what this does?
  done
done

# Clean up unneeded intermediate files?
#rm Processing_$amp/*discarded* 
#rm Processing_$amp/*unassembled*
#rm Processing_$amp/*filtered*
#rm Processing_$amp/*uniques*
#rm Processing_$amp/*denoise.fasta*

mkdir Combined_amplicons # Make folder for merged sequence results
# Combined LCO/BR denoised sequences into one file for later use
cat  Processing_LCO1490-FC/*.denoise.good.fasta > Processing_LCO1490-FC/ALL_LCO.denoise.good.fasta
cat  Processing_BR-HCO2198/*.denoise.good.fasta > Processing_BR-HCO2198/ALL_BR.denoise.good.fasta
cat  Processing_LCO1490-FC/*.denoise.topn.fasta > Processing_LCO1490-FC/ALL_LCO.denoise.topn.fasta
cat  Processing_BR-HCO2198/*.denoise.topn.fasta > Processing_BR-HCO2198/ALL_BR.denoise.topn.fasta
# Add amplicon labels to sequence headers, to allow easy matching with alignment results
sed 's/>/>LCO./g' Processing_LCO1490-FC/ALL_LCO.denoise.good.fasta > Combined_amplicons/ALL_LCO.denoise.good.fasta
sed 's/>/>BR./g' Processing_BR-HCO2198/ALL_BR.denoise.good.fasta > Combined_amplicons/ALL_BR.denoise.good.fasta
sed 's/>/>LCO./g' Processing_LCO1490-FC/ALL_LCO.denoise.topn.fasta > Combined_amplicons/ALL_LCO.denoise.topn.fasta
sed 's/>/>BR./g' Processing_BR-HCO2198/ALL_BR.denoise.topn.fasta > Combined_amplicons/ALL_BR.denoise.topn.fasta

###############################################################################
# Identify optimal pairwise alignments between denoised sequences
for i in $(seq -f "%03g" 1 450); do
  sp_id=i$i
  if [[ -f "Processing_LCO1490-FC/"$sp_id".denoise.topn.fasta" && -f "Processing_BR-HCO2198/"$sp_id".denoise.topn.fasta" ]]; then # LCO and BR both found
    echo $sp_id: Both files exist. Finding alignments...
    # Add amplicon labels to sequence labels (so they can be distinguished later)
    sed 's/>/>LCO./g' Processing_LCO1490-FC/$sp_id.denoise.topn.fasta > Combined_amplicons/LCO_$sp_id.denoise.topn.fasta
    sed 's/>/>BR./g' Processing_BR-HCO2198/$sp_id.denoise.topn.fasta > Combined_amplicons/BR_$sp_id.denoise.topn.fasta
    # Combine LCO and BR sequences into one file
    cat Combined_amplicons/LCO_$sp_id.denoise.topn.fasta Combined_amplicons/BR_$sp_id.denoise.topn.fasta > Combined_amplicons/both_$sp_id.denoise.topn.fasta
    # Use allpairs_global to identify optimal alignments (overlap of 85 bp and 100 % identity). Doesn't output merged sequences though...
    vsearch -allpairs_global Combined_amplicons/both_$sp_id.denoise.topn.fasta -id 1 -mincols 85 -blast6out Combined_amplicons/both_$sp_id.denoise.topn.matches.txt
    # If no alignments found among topn sequences, check for alignments among ALL sequences
    gr1="$( grep -c LCO Combined_amplicons/both_$sp_id.denoise.topn.matches.txt )"
    if [[ $gr1 -eq 0 ]]; then
      gr2="$( grep -c \> Processing_LCO1490-FC/$sp_id.denoise.good.fasta )"
      gr3="$( grep -c \> Processing_BR-HCO2198/$sp_id.denoise.good.fasta )"
      if [[ $gr2 -gt $d_count || $gr3 -gt $d_count ]]; then
        echo $sp_id: No alignments found among topn sequences. Checking for alignments among ALL sequences...
        sed 's/>/>LCO./g' Processing_LCO1490-FC/$sp_id.denoise.good.fasta > Combined_amplicons/LCO_$sp_id.denoise.good.fasta
        sed 's/>/>BR./g' Processing_BR-HCO2198/$sp_id.denoise.good.fasta > Combined_amplicons/BR_$sp_id.denoise.good.fasta
        cat Combined_amplicons/LCO_$sp_id.denoise.good.fasta Combined_amplicons/BR_$sp_id.denoise.good.fasta > Combined_amplicons/both_$sp_id.denoise.good.fasta
        vsearch -allpairs_global Combined_amplicons/both_$sp_id.denoise.good.fasta -id 1 -mincols 85 -blast6out Combined_amplicons/both_$sp_id.denoise.ALL.matches.txt
        # If alignments found among ALL sequences, output combined sequences for BLAST
        gr4="$( grep -c LCO Combined_amplicons/both_$sp_id.denoise.ALL.matches.txt )"
        if [[ $gr4 -gt 0 ]]; then # 
          cat Combined_amplicons/LCO_$sp_id.denoise.good.fasta >> Combined_amplicons/TopnNoAlign_LCO.denoise.good.fasta
          cat Combined_amplicons/BR_$sp_id.denoise.good.fasta >> Combined_amplicons/TopnNoAlign_BR.denoise.good.fasta
        fi
      fi
    fi
  elif [[ -f "Processing_LCO1490-FC/"$sp_id".denoise.topn.fasta" && ! -f "Processing_BR-HCO2198/"$sp_id".denoise.topn.fasta" ]]; then # Only LCO found
    echo $sp_id: Found LCO file only...
    #sed 's/>/>LCO./g' Processing_LCO1490-FC/$sp_id.denoise.topn.fasta > Combined_amplicons/LCO_only_$sp_id.denoise.topn.fasta
  elif [[ ! -f "Processing_LCO1490-FC/"$sp_id".denoise.topn.fasta" && -f "Processing_BR-HCO2198/"$sp_id".denoise.topn.fasta" ]]; then # Only BR found
    echo $sp_id: Found BR file only...
    #sed 's/>/>BR./g' Processing_BR-HCO2198/$sp_id.denoise.topn.fasta > Combined_amplicons/BR_only_$sp_id.denoise.topn.fasta
  elif [[ ! -f "Processing_LCO1490-FC/"$sp_id".denoise.topn.fasta" && ! -f "Processing_BR-HCO2198/"$sp_id".denoise.topn.fasta" ]]; then # Neither LCO or BR found
    echo $sp_id: Neither file found...
  fi
done

###############################################################################
# Later, after sequences assembled into barcodes
cat Barcode_picks/*fasta > Barcode_picks/All_merged_plus_LCO_BR_only.fasta
# To make a phylogeny of barcode sequences
module load MAFFT/7.429-gimkl-2020a
module load FastTree/2.1.11-GCCcore-9.2.0
#mafft --auto All_merged_plus_LCO_BR_only_barcodes.fasta > All_merged_plus_LCO_BR_only_barcodes.mafft
# L-INS-1 method (slow but accurate; acceptable speed for ~450 sequennces)
mafft --maxiterate 1000 --localpair Barcode_picks/All_merged_plus_LCO_BR_only.fasta > Barcode_picks/All_merged_plus_LCO_BR_only.mafft
FastTree -nt Barcode_picks/All_merged_plus_LCO_BR_only.mafft > Barcode_picks/All_merged_plus_LCO_BR_only.fasttree