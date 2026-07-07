#!/bin/bash
set -euo pipefail
# Redirect both stdout (the echo progress lines) and stderr into the
# Snakemake log file, instead of only stderr.
exec > "${snakemake_log[0]}" 2>&1

VCF_FILE="${snakemake_input[vcf]}"
OUTPUT_FILE="${snakemake_output[tsv]}"

echo "VCF_FILE: $VCF_FILE"

(
echo -e "CHROM\tPOS\tREF\tALT\tSample\tGT\tGene\tAnnotation\tImpact\tHGVS.c\tHGVS.p\tAA_Change\tVariant_Class\tLOF\tNMD"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/LOF\t%INFO/NMD\t[%SAMPLE=%GT\t]\n' "$VCF_FILE" | \
awk -F'\t' -v OFS='\t' '

function abbr(p,    s, aa, map) {
    if (p == "" || p !~ /^p\./) return p;

    map["Ala"]="A"; map["Arg"]="R"; map["Asn"]="N"; map["Asp"]="D"; map["Cys"]="C";
    map["Gln"]="Q"; map["Glu"]="E"; map["Gly"]="G"; map["His"]="H"; map["Ile"]="I";
    map["Leu"]="L"; map["Lys"]="K"; map["Met"]="M"; map["Phe"]="F"; map["Pro"]="P";
    map["Ser"]="S"; map["Thr"]="T"; map["Trp"]="W"; map["Tyr"]="Y"; map["Val"]="V";
    map["Ter"]="*";

    s = p;
    sub(/^p\./, "", s);

    for (aa in map) gsub(aa, map[aa], s);

    return s;
}

function parse_effects(field, arr,    n,i,tmp,gene) {
    delete arr;
    if (field == "." || field == "") return;

    n = split(field, tmp, ",");

    for (i = 1; i <= n; i++) {
        gsub(/[()]/, "", tmp[i]);
        split(tmp[i], f, "|");
        gene = f[1];
        arr[gene] = 1;
    }
}

{
    chrom = $1;
    pos = $2;
    ref = $3;
    alt_string = $4;
    ann_string = $5;
    lof_string = $6;
    nmd_string = $7;

    split(alt_string, alt_alleles, ",");
    split(ann_string, anns, ",");

    parse_effects(lof_string, lof_genes);
    parse_effects(nmd_string, nmd_genes);

    for (j = 8; j <= NF; j++) {
        split($j, s, "=");
        sample = s[1];
        gt = s[2];

        split(gt, gts, /[\/|]/);

        n_nonref = 0;
        for (k in gts)
            if (gts[k] != "0" && gts[k] != ".") n_nonref++;

        if (n_nonref == 0) continue;

        delete seen;

        for (k in gts) {
            a = gts[k];
            if (a == "0" || a == "." || (a in seen)) continue;
            seen[a] = 1;

            for (i = 1; i <= length(anns); i++) {
                split(anns[i], f, "|");

                ann_alt = f[1];
                annotation = f[2];
                impact = f[3];
                gene = f[4];
                hgvsc = f[10];
                hgvsp = f[11];

                if (ann_alt != alt_alleles[a]) continue;

                # Variant class
                if (length(ref) < length(alt_alleles[a])) variant_class="Insertion";
                else if (length(ref) > length(alt_alleles[a])) variant_class="Deletion";
                else variant_class="SNP";

                short = abbr(hgvsp);

                # LOF / NMD detection (no percentages)
                lof_flag = (gene in lof_genes) ? 1 : 0;
                nmd_flag = (gene in nmd_genes) ? 1 : 0;

                print chrom, pos, ref, alt_alleles[a], sample, gt,
                      gene, annotation, impact,
                      hgvsc, hgvsp, short,
                      variant_class,
                      lof_flag,
                      nmd_flag;
            }
        }
    }
}
' | sort -k5,5 ) > "$OUTPUT_FILE"
