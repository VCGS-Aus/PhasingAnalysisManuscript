#!/usr/bin/env bash

# calculate_heterozygosity.sh manifest.txt output.tsv
# manifest.txt: sample_id /path/to/vcf

input_tsv=$1

output_tsv=$2

bed_file="../data/mendeliome_incidentalome_no_overlap_autosomes.bed"

echo -e "ID\tChromosome\tHet_Count\tHom_Count\tHeterozygosity" > "$output_tsv"

while IFS=$'\t' read -r ID VCF; do

    [[ -z "$ID" || -z "$VCF" ]] && continue
     /misc/vcgs/seq/cpipe-wgs-staging/tools/bcftools/1.23/bcftools index VCF
    /misc/vcgs/seq/cpipe-wgs-staging/tools/bcftools/1.23/bcftools query -R "$bed_file" -f '%CHROM\t[%GT\n]' "$VCF" | \
    awk -v id="$ID" '
    {
        chr = $1
        gt  = $2

        if (gt == "./." || gt == ".|.") next

        if (gt ~ /0[\/|]1|1[\/|]0/) {
            het[chr]++
        }
        else if (gt == "1/1" || gt == "1|1") {
            hom[chr]++
        }

        seen[chr] = 1
    }
    END {
        for (i = 1; i <= 22; i++) {
            c = "chr" i
            if (!(c in seen)) continue

            h = (c in het) ? het[c] : 0
            m = (c in hom) ? hom[c] : 0

            total = h + m
            heterozygosity = (total > 0) ? h/total : "NA"

            printf "%s\t%s\t%d\t%d\t%.6f\n", id, c, h, m, heterozygosity
        }
    }' >> "$output_tsv"

done < "$input_tsv"
