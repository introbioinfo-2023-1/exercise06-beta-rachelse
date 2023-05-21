# command02.sh
# 1. Download and setup Swiss-Prot database for MMSeqs in the ./db folder.

# 2. Run easy-search workflow for unknown_bacteria.proteins.faa against Swiss-Prot
#    and save the result as ./result/unknown_bacteria.proteins.aln.swissprot.tsv

# 3. Run Mmseqs2 easy-taxonomy workflow for unknown_bacteria.proteins.faa and
#    save the results with the prefix ./result/unknown_bacteria.taxonomy
#    Identify the strain of the unknown genome with ./result/unknown_bacteria.taxonomy_report and
#    save the scietific name of the strain at ./result/unknown_bacteria.strain.txt

# 4. Repeat easy-search workflow analysis with ./data/unknown_transcripts.fasta as in Step 2.
#    Save the result as ./result/unknown_transcripts.aln.swissprot.tsv and copy
#    the lines of best hits of two transcripts to ./result/unknown_transcripts.proteins.txt