#!/bin/bash
# ==========================================================
# 🇪🇹 Ethiopian Marburg virus - GenBank Annotation Script (v13)
# 🎯 TARGET: Nature Medicine submission
# 🧬 WORKFLOW: iVar v1.3.1 | Illumina MiSeq | 12,450x
# ==========================================================

set -euo pipefail

CONSENSUS_DIR="../results/07_consensus"
REFERENCE_GB="../reference_genomes/Marburg_reference.gb"

# Define the lead submitters for Reference 2
SUBMISSION_AUTHORS="Feyssa,M.D., Zerihun.B and Wolday,D."

for fasta in "$CONSENSUS_DIR"/*.fa; do
    base=$(basename "$fasta" .fa)
    dest_gb="$CONSENSUS_DIR/${base}.gb"

    # Extract sequence and length
    seq=$(grep -v '^>' "$fasta" | tr -d '\n' | tr -d '\r')
    seq_len=${#seq}
    today=$(date +%d-%b-%Y | tr 'a-z' 'A-Z')

    # 1. Feature Extraction (Strictly skips Reference Headers)
    awk -v seqlen="$seq_len" -v sample="$base" '
    BEGIN {source_done=0; in_translation=0; in_features=0}
    
    /^FEATURES/ {in_features=1; next}
    /^ORIGIN/ {in_features=0; next}
    
    in_features == 1 {
        # Update Source modifier with Ethiopian Metadata
        if (/^ {5}source/ && source_done==0) {
            print "     source          1.." seqlen
            print "                     /organism=\"Orthomarburgvirus marburgense\""
            print "                     /mol_type=\"viral cRNA\""
            printf "                     /isolate=\"MARV|H.Sapiens|ETH_Jinka|%s|2025\"\n", sample
            print "                     /isolation_source=\"blood\""
            print "                     /host=\"Homo sapiens\""
            print "                     /db_xref=\"taxon:3052505\""
            print "                     /geo_loc_name=\"Ethiopia: Jinka\""
            print "                     /collection_date=\"Nov-2025\""
            source_done=1
            while(getline && !/^ {5}[a-z]/) {} 
        }

        # Update Feature coordinates (Gene/CDS)
        if (/^ {5}[a-z]/) {
            in_translation=0
            feat_type = $1; split($2, pos, "\\.\\.")
            start = pos[1]; end = pos[2]
            gsub(/[^0-9]/, "", start); gsub(/[^0-9]/, "", end)
            
            # Skip features that are entirely outside the consensus range
            if (start > seqlen) { next }
            new_end = (end > seqlen ? seqlen : end)
            
            printf "     %-15s %d..%d\n", feat_type, start, new_end
            next
        }

        # Handle Locus Tags and remove old protein/GI tags
        if (/^ {21}\//) {
            if ($0 ~ /\/locus_tag=/) {
                match($0, /[0-9]+/, arr); gene_num = arr[0]
                printf "                     /locus_tag=\"MARVETH_%02d\"\n", gene_num
                next
            }
            if ($0 ~ /\/translation=/) { in_translation=1; next }
            if (in_translation == 1 && $0 ~ /^ {21}[^/]/) { next }
            if (in_translation == 1 && $0 ~ /^ {21}\// && !($0 ~ /\/translation=/)) { in_translation=0 }
            if ($0 ~ /\/protein_id/ || $0 ~ /\/db_xref="GeneID/ || $0 ~ /\/db_xref="GI/) { next }
            
            if (in_translation == 0) { print $0 }
            next
        }
        if (in_translation == 0) print $0
    }' "$REFERENCE_GB" > "$dest_gb.tmp_features"

    # 2. Construct the GenBank Flatfile
    {
        printf "LOCUS       %-15s %d bp    cRNA    linear   VRL %s\n" "$base" "$seq_len" "$today"
        echo "DEFINITION  Orthomarburgvirus marburgense isolate"
        echo "            MARV|H.Sapiens|ETH_Jinka|${base}|2025, complete genome."
        echo "ACCESSION   $base"
        echo "VERSION     ${base}.1"
        echo "KEYWORDS    ."
        echo "SOURCE      Orthomarburgvirus marburgense"
        echo "  ORGANISM  Orthomarburgvirus marburgense"
        echo "            Viruses; Riboviria; Orthornavirae; Negarnaviricota;"
        echo "            Haploviricotina; Monjiviricetes; Mononegavirales; Filoviridae;"
        echo "            Orthomarburgvirus."
        echo "REFERENCE   1  (bases 1 to ${seq_len})"
        echo "  AUTHORS   Feyssa,M.D., Tollera,G., Abdella,S., Bashea,C., Zerihun,B., Getu,M., Tasew,G., Tessema,M. and Wolday,D."
        echo "  TITLE     Genomic characterization and transmission dynamics of the first"
        echo "            cases of the 2025 Marburg virus outbreak in Ethiopia"
        echo "  JOURNAL   Nature Medicine"
        echo "REFERENCE   2  (bases 1 to ${seq_len})"
        echo "  AUTHORS   $SUBMISSION_AUTHORS"
        echo "  TITLE     Direct Submission"
        echo "  JOURNAL   Submitted ($(date +%d-%b-%Y)) Ethiopian Public Health Institute (EPHI),"
        echo "            Swaziland Street, Addis Ababa, Ethiopia"
        echo "COMMENT     ##Genome-Assembly-Data-START##"
        echo "            Assembly Date         :: 15-Dec-2025"
        echo "            Assembly Method       :: iVar v1.3.1"
        echo "            Genome Coverage       :: 12,450x"
        echo "            Sequencing Technology :: Illumina MiSeq"
        echo "            ##Genome-Assembly-Data-END##"
        echo "FEATURES             Location/Qualifiers"
        cat "$dest_gb.tmp_features"
        echo "ORIGIN"
        echo "$seq" | fold -w60 | awk '{printf "%9d %s\n", (NR-1)*60+1, $0}'
        echo "//"
    } > "$dest_gb"

    rm "$dest_gb.tmp_features"
done

echo "✅ Annotation complete. Files are in $CONSENSUS_DIR"
