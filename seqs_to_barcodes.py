# -*- coding: utf-8 -*-

# Python code to find target LCO1490-FC and BR-HCO2198 sequences and combine into full length barcodes.
# Requires Python 3.8 with Biopython, and EMBOSS to be available.
# Input data requirements are:
# - a tab-delimited table of specimen taxonomic identification details.
# - a tab-delimited table of top BLAST results for each denoised FC and BR sequence.
# - rankedlineage.dmp from NCBI taxdump.
# - denoised sequence files
module load Python/3.8.2-gimkl-2020a
module load EMBOSS/6.6.0-gimkl-2020a
# vsearch = "/nesi/project/landcare00060/bioinformatics/vsearch-2.21.1-linux-x86_64/bin/vsearch" # not used

import os, subprocess, re
from Bio import SeqIO, SeqRecord

#os.chdir("Processed_data/Demultiplexed/")

# Load new NCBI taxdump ranked lineages for filtering of BLAST results
rankedlineages = {}
with open("rankedlineage.dmp", "r", encoding='utf-8') as txranks:
    print("Loading tax dump ranked lineages...")
    # taxid|name|species|genus|family|order|class|phylum|kingdom|superkingdom
    for row in txranks:
        bits = row.split("\t|\t")
        bits[len(bits)-1] = bits[len(bits)-1].split("\t")[0]
        rankedlineages[bits[0]] = row
    print("...done.")

# Load denoised sequences and BLAST results
seqs_LCO = {}
with open("Combined_amplicons/ALL_LCO.denoise.good.fasta","r") as seqfile:
    print("Loading LCO1490-FC denoised sequences...")
    for rec in SeqIO.parse(seqfile, "fasta"):
        seqs_LCO[rec.id] = rec.seq
    print("...done.")

seqs_BR = {}
with open("Combined_amplicons/ALL_BR.denoise.good.fasta","r") as seqfile:
    print("Loading BR-HCO2198 denoised sequences...")
    for rec in SeqIO.parse(seqfile, "fasta"):
        seqs_BR[rec.id] = rec.seq
    print("...done.")

# Load BLAST results and get taxonomy details for each sequence 
def load_BLAST_results(BLASTfile):
    BLAST_res = {}
    for row in BLASTfile:
        bits = row.split("\t")
        seq_id, taxid = bits[0], bits[10].split(";")[0]
        lineage = rankedlineages.get(taxid, "None").split("\t|\t")
        BLAST_details = "{}\t{}\t{}\t{}\t{}".format(bits[3],bits[6],bits[7],bits[9],bits[11])
        if lineage != ["None"]:
            lineage_details = "{}\t{}\t{}\t{}".format(lineage[4],lineage[5],lineage[6],lineage[7])
        else:
            lineage_details = "\t\t\t"
        #print("{}: {} {}".format(seq_id, BLAST_details, lineage_details))
        BLAST_res[seq_id] = "{}\t{}\t{}".format(seq_id, BLAST_details, lineage_details)
    return(BLAST_res)

# Load BLAST results for topn sequences
with open("Combined_amplicons/ALL_LCO.denoise.topn.fasta.BLASTx1.txt","r") as BLASTfile:
    print("Loading LCO1490-FC BLAST results...")
    BLAST_LCO = load_BLAST_results(BLASTfile)
    print("...done.")

with open("Combined_amplicons/ALL_BR.denoise.topn.fasta.BLASTx1.txt","r") as BLASTfile:
    print("Loading BR-HCO2198 BLAST results...")
    BLAST_BR = load_BLAST_results(BLASTfile)
    print("...done.")

# Load BLAST results from specimens with alignments among all sequences but not among topn sequences
with open("Combined_amplicons/TopnNoAlign_LCO.denoise.good.fasta.BLASTx1.txt","r") as BLASTfile:
    print("Loading LCO1490-FC BLAST results for ALL sequences, where no alignments found among topn...")
    BLAST_LCO2 = load_BLAST_results(BLASTfile)
    print("...done.")

with open("Combined_amplicons/TopnNoAlign_LCO.denoise.good.fasta.BLASTx1.txt","r") as BLASTfile:
    print("Loading BR-HCO2198 BLAST results for ALL sequences, where no alignments found among topn...")
    BLAST_BR2 = load_BLAST_results(BLASTfile)
    print("...done.")

# Merge the two sets of BLAST results
BLAST_LCO_both = {**BLAST_LCO, **BLAST_LCO2}
BLAST_BR_both = {**BLAST_BR, **BLAST_BR2}

        
# Load taxonomy reference details for specimens
tax_details = {}
i = 1
with open("Taxonomy_details.txt", "r") as tx_ref:
    print("Loading specimen taxonomy details...")
    for row in tx_ref:
        if i > 1:
            bits = row.split("\t")
            bits_tx = [bit.lower().capitalize().strip() for bit in bits if bits.index(bit) > 1]
            tax_details[bits[0]] = bits_tx
        i += 1
    print("...done.")

# Load pairwise alignment details
alignments = {}
print("Loading alignments...")
for i in range(1, 450+1):
    fALL = "Combined_amplicons/both_i{:0>3d}.denoise.ALL.matches.txt".format(i)
    ftopn = "Combined_amplicons/both_i{:0>3d}.denoise.topn.matches.txt".format(i)
    if os.path.isfile(fALL) or os.path.isfile(ftopn):
        if os.path.isfile(fALL):
            alignf = open(fALL, "r")
        elif os.path.isfile(ftopn) and not os.path.isfile(fALL):
            alignf = open(ftopn, "r")
        matches = []
        for row in alignf:
            LCO_id = row.split("\t")[0]#.split(".", 1)[1]
            BR_id = row.split("\t")[1]#.split(".", 1)[1]
            LCO_tx = BLAST_LCO_both.get(LCO_id, "None")
            BR_tx = BLAST_BR_both.get(BR_id, "None")
            matches.append([LCO_id, LCO_tx, BR_id, BR_tx])
        alignments["i{:0>3d}".format(i)] = matches

# Select target sequences based on BLAST results and combine
# i001-i222 & i355-i450 = insects; i223-i354 = worms
def pick_seq_by_taxonomy(BLAST_sp, sp_tx, sp_id, amp, breakties = "abundances"):
    picks = ""
    if amp == "LCO":
        thr = 200
    elif amp == "BR":
        thr = 250
    print("Finding {} sequence for {}...".format(amp, sp_id))
    if "" in sp_tx:
        sp_tx.remove("") # Remove any empty taxonomy values
    for target in reversed(sp_tx): # Search for taxonomy identification details from most to least specific rank
        picks = [v for v in BLAST_sp if re.search(target, v) and float(v.split("\t")[4]) > thr]
        if len(picks) > 0:
            print("Target found: {}".format(target))
            if len(picks) == 1:
                pick = picks[0]
            elif len(picks) > 1: # If > 1 match to target, take the one with the max abundance, or bitscore?
                if breakties == "bitscores":
                    bitscores = [float(pick.split("\t")[4]) for pick in picks]
                    pick = [pick for pick in picks if float(pick.split("\t")[4]) == max(bitscores)]
                elif breakties == "abundances":
                    abundances = [float(pick.split("\t")[0].split("size=")[1]) for pick in picks]
                    pick = [pick for pick in picks if float(pick.split("\t")[0].split("size=")[1]) == max(abundances)]
                if len(pick) > 0: # If > 1 match with the same target AND bitscores/abundances (unlikely), take the first one
                    pick = pick[0]
            break
    if len(picks) == 0:
        print("No target sequence found")
        pick = "No target sequence"
    return(pick)

def find_alignments(sp_id, use_tx = True, breakties = "abundances"):
    print("Finding alignments for {}...".format(sp_id))
    pick1 = "None"
    pick2 = "None"
    sp_tx = tax_details[sp_id] # Get specimen taxonomy details
    if "" in sp_tx:
        sp_tx.remove("") # Remove any empty taxonomy values
    sp_align = alignments.get(sp_id, "None")
    if sp_align != "None":
        # Select alignment with expected class identifications and highest abundance
        for target in reversed(sp_tx[1:5]):
            keep_al = [al for al in sp_align if \
                    (re.search(target, al[1]) and int(al[1].split("\t")[4]) >= 200) or \
                    (re.search(target, al[3]) and int(al[3].split("\t")[4]) >= 250) and \
                    324 <= len(seqs_LCO[al[0]]) <= 326 and 417 <= len(seqs_BR[al[2]]) <= 419 ]
            if len(keep_al) > 0:
                # If more than one alignment found at given taxonomic rank level,
                # take the one with highest mean abundance, or highest mean bitscore?
                if len(keep_al) > 1:
                    if breakties == "abundances":
                        sizes = []
                        for al in keep_al:
                            LCO_size = float(al[0].split("=")[1])
                            BR_size = float(al[2].split("=")[1])
                            sizes.append((LCO_size + BR_size)/2)
                        pick = keep_al[sizes.index(max(sizes))]
                    elif breakties == "bitscores":
                        scores = []
                        for al in keep_al:
                            # Because a few sequences have missing taxonomy/no bitscore:
                            if al[1] == "None":
                                LCO_score = 0
                            else:
                                LCO_score = int(al[1].split("\t")[4])
                            if al[3] == "None":
                                BR_score = 0
                            else:
                                BR_score = int(al[3].split("\t")[4])
                            scores.append((LCO_score + BR_score)/2)
                        pick = keep_al[scores.index(max(scores))]
                elif len(keep_al) == 1:
                    pick = keep_al[0]
                break
        if len(keep_al) == 0:
            pick = "None"
        if pick != "None":
            print("{} alignment pick: {} + {}.".format(sp_id, pick[0], pick[2]))
        elif pick == "None":
            print("{}: no alignment picks found".format(sp_id))
        return(pick)
    elif sp_align == "None":
        print("{}: no alignments found".format(sp_id))
        return("None")

def merge_seqs(LCO_pick_id, BR_pick_id):
    LCO_seq = SeqRecord.SeqRecord(seqs_LCO[LCO_pick_id], id = sp_id, description = "{} + {}".format(LCO_pick_id, BR_pick_id)) 
    BR_seq = SeqRecord.SeqRecord(seqs_BR[BR_pick_id], id = sp_id, description = "{} + {}".format(LCO_pick_id, BR_pick_id))
    SeqIO.write(LCO_seq, "Combined_amplicons/{}_LCO_pick.fasta".format(sp_id), "fasta")
    SeqIO.write(BR_seq, "Combined_amplicons/{}_BR_pick.fasta".format(sp_id), "fasta")
    # Merge sequence picks using EMBOSS merger 
    merge = ["merger -asequence=Combined_amplicons/{0}_LCO_pick.fasta " \
                    "-bsequence=Combined_amplicons/{0}_BR_pick.fasta " \
                    "-outfile=Barcode_picks/{0}_align.txt " \
                    "-outseq=Barcode_picks/{0}_merged.fasta".format(sp_id)]
    subprocess.run(merge, shell=True)
    # Obtain pairwise alignment details
    with open("Barcode_picks/{0}_align.txt".format(sp_id), "r") as align:
        stats = [line for line in align if re.search("Length|Identity|Gaps|Score", line)]
        mlength = stats[0].split(":")[1].strip()
        mid = stats[1].split(":")[1].strip()
        mgaps = stats[2].split(":")[1].strip()
        mscore = stats[3].split(":")[1].strip()
        merge_stats = [mlength, mid, mgaps, mscore]
    return(merge_stats)
 
os.mkdir("Barcode_picks/") # Make a new folder for final barcode output
bt = "abundances"
with open("Barcode_picks/Seqs_to_barcodes_summary.txt", "w") as summary:
    summary.write("Specimen\t"\
                  "LCO_seq\tLCO_GI\tLCO_length\tLCO_id\tLCO_bitscore\tLCO_species\t"\
                  "LCO_Family\tLCO_Order\tLCO_Class\tLCO_Phylum\t"\
                  "BR_seq\tBR_GI\tBR_Length\tBR_ID\tBR_bitscore\tBR_species\t"\
                  "BR_Family\tBR_Order\tBR_Class\tBR_Phylum\t"\
                  "Outcome\tLength\tIdentity\tGaps\tScore\n")
    merge_count = 0
    LCO_only_count = 0
    BR_only_count = 0
    both_nomerge_count = 0
    neither_count = 0
    merge_stats = [0,0,0,0]
    for i in range(1, 450+1):
        sp_id = "i{:0>3d}".format(i)
        sp_tx = tax_details[sp_id]
        LCO_BLAST_sp = [BLAST_LCO_both[k] for k in BLAST_LCO_both.keys() if k.split(".")[1] == sp_id]
        LCO_pick_detail = pick_seq_by_taxonomy(LCO_BLAST_sp, sp_tx, sp_id, amp = "LCO", breakties = bt)
        LCO_pick_id = LCO_pick_detail.split("\t")[0]
        BR_BLAST_sp = [BLAST_BR_both[k] for k in BLAST_BR_both.keys() if k.split(".")[1] == sp_id]
        BR_pick_detail = pick_seq_by_taxonomy(BR_BLAST_sp, sp_tx, sp_id, amp = "BR", breakties = bt)
        BR_pick_id = BR_pick_detail.split("\t")[0]
        sp_al = find_alignments(sp_id, breakties = bt)
        # Check if non-aligned picks are better (i.e. have more specific taxonomy) than alignment
        if LCO_pick_detail != "No target sequence":
            LCO_pick_ranks = [sp_tx.index(tx) for tx in reversed(sp_tx) if tx in LCO_pick_detail]
        else: # No target LCO sequence found
            LCO_pick_ranks = [0]
        if BR_pick_detail != "No target sequence":
            BR_pick_ranks = [sp_tx.index(tx) for tx in reversed(sp_tx) if tx in BR_pick_detail]
        else: # No target BR sequence found
            BR_pick_ranks = [0]
        if sp_al != "None":
            al_LCO_ranks = [sp_tx.index(tx) for tx in reversed(sp_tx) if tx in sp_al[1]]
            al_BR_ranks = [sp_tx.index(tx) for tx in reversed(sp_tx) if tx in sp_al[3]]
            al_ranks = al_LCO_ranks + al_BR_ranks
        else: # No alignment found
            al_ranks = [0]
        if max(al_ranks) >= max(LCO_pick_ranks + BR_pick_ranks) and max(al_ranks) > 0: # Alignment better
            print("Combining seqs for {}: {}".format(sp_id, sp_tx))
            LCO_pick_id_al = sp_al[0]
            LCO_pick_detail_al = sp_al[1]
            BR_pick_id_al = sp_al[2]
            BR_pick_detail_al = sp_al[3]
            merge_stats = merge_seqs(LCO_pick_id_al, BR_pick_id_al)
            print("{}: merged length = {}, score = {}".format(sp_id, merge_stats[0], merge_stats[3]))
            merge_count += 1
            outcome = "Aligned"
        elif max(LCO_pick_ranks + BR_pick_ranks) > max(al_ranks): # Un-aligned sequence pick(s) better
            merge_stats = [0,0,0,0]
            if max(LCO_pick_ranks) > 0 and max(BR_pick_ranks) == 0: # LCO found but not BR
                LCO_only_count += 1
                LCO_seq = SeqRecord.SeqRecord(seqs_LCO[LCO_pick_id], id = LCO_pick_id, description = LCO_pick_id)
                SeqIO.write(LCO_seq, "Barcode_picks/{}_LCO_only.fasta".format(sp_id), "fasta")
                #BR_pick_detail = "No target sequence\t\t\t\t\t\t\t\t\t"
                outcome =  "LCO1490-FC sequence only"
            elif max(BR_pick_ranks) > 0 and max(LCO_pick_ranks) == 0: # BR found but not LCO
                BR_only_count += 1
                BR_seq = SeqRecord.SeqRecord(seqs_BR[BR_pick_id], id = BR_pick_id, description = BR_pick_id)
                SeqIO.write(BR_seq, "Barcode_picks/{}_BR_only.fasta".format(sp_id), "fasta")
                #LCO_pick_detail = "No target sequence\t\t\t\t\t\t\t\t\t"
                outcome = "BR-HCO2198 sequence only"
            elif max(LCO_pick_ranks) > 0 and max(BR_pick_ranks) > 0: # Both sequence picks found
                both_nomerge_count += 1
                if max(LCO_pick_ranks) > max(BR_pick_ranks): # LCO pick has more specific taxonomy
                    outcome = "No alignments, take LCO1490-FC sequence"
                    LCO_seq = SeqRecord.SeqRecord(seqs_LCO[LCO_pick_id], id = LCO_pick_id, description = LCO_pick_id)
                    SeqIO.write(LCO_seq, "Barcode_picks/{}_LCO_only.fasta".format(sp_id), "fasta")
                elif max(BR_pick_ranks) > max(LCO_pick_ranks): # BR pick has more specific taxonomy
                    outcome = "No alignments, take BR-HCO2198 sequence"
                    BR_seq = SeqRecord.SeqRecord(seqs_BR[BR_pick_id], id = BR_pick_id, description = BR_pick_id)
                    SeqIO.write(BR_seq, "Barcode_picks/{}_BR_only.fasta".format(sp_id), "fasta")
                elif max(BR_pick_ranks) == max(LCO_pick_ranks): # LCO and BR have equivalent taxonomy
                    picks_both = [LCO_pick_detail, BR_pick_detail]
                    if bt == "bitscores": # Pick sequence with highest bitscore
                        bitscores = [int(pick.split("\t")[4]) for pick in picks_both if pick != "None"]
                        pick = [pick for pick in picks_both if int(pick.split("\t")[4]) == max(bitscores)]
                    elif bt == "abundances": # Pick sequence with highest abundance
                        abundances = [float(pick.split("\t")[0].split("size=")[1]) for pick in picks_both]
                        pick = [pick for pick in picks_both if float(pick.split("\t")[0].split("size=")[1]) == max(abundances)][0]
                    if "LCO" in pick: # LCO picked
                        outcome = "No alignments, take LCO1490-FC sequence"
                        LCO_seq = SeqRecord.SeqRecord(seqs_LCO[LCO_pick_id], id = LCO_pick_id, description = LCO_pick_id)
                        SeqIO.write(LCO_seq, "Barcode_picks/{}_LCO_only.fasta".format(sp_id), "fasta")
                    elif "BR" in pick: # BR picked
                        outcome = "No alignments, take BR-HCO2198 sequence"
                        BR_seq = SeqRecord.SeqRecord(seqs_BR[BR_pick_id], id = BR_pick_id, description = BR_pick_id)
                        SeqIO.write(BR_seq, "Barcode_picks/{}_BR_only.fasta".format(sp_id), "fasta")
        elif max(LCO_pick_ranks) == 0 and max(BR_pick_ranks) == 0: # Neither LCO or BR found
            neither_count += 1
            #LCO_pick_detail = "No target sequence\t\t\t\t\t\t\t\t\t"
            #BR_pick_detail = "No target sequence\t\t\t\t\t\t\t\t\t"
            outcome = "Neither target sequence"
        if LCO_pick_detail == "No target sequence":
            LCO_pick_detail = "No target sequence\t\t\t\t\t\t\t\t\t"
        if BR_pick_detail == "No target sequence":
            BR_pick_detail = "No target sequence\t\t\t\t\t\t\t\t\t"
        merge_bits = "{}\t{}\t{}\t{}".format(merge_stats[0],merge_stats[1],merge_stats[2],merge_stats[3])
        summary.write("{}\t{}\t{}\t{}\t{}\n".format(sp_id, LCO_pick_detail, BR_pick_detail, outcome, merge_bits))

print("Finished: {} merged barcodes, {} LCO only, {} BR only, {} both but no alignment, {} neither.".format(merge_count, LCO_only_count, BR_only_count, both_nomerge_count, neither_count))