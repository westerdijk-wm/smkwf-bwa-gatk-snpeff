#!/usr/bin/env bash
#
# Flatten a SnpEff-annotated VCF into a per-sample TSV table and split into
# batches of BATCH_SIZE samples. Designed for use with Snakemake's `script:`
# directive using a directory() output.
#
# NOTE: handles both diploid GT ("0/1", "0|1") and haploid GT ("1", "0", ".")
# by explicitly detecting ploidy per-genotype and emitting it as its own
# column, instead of relying on downstream code to guess it from the shape
# of the raw GT string.

set -euo pipefail

# Redirect both stdout (the echo progress lines) and stderr into the
# Snakemake log file, instead of only stderr.
exec > "${snakemake_log[0]}" 2>&1

VCF_FILE="${snakemake_input[vcf]}"
OUTDIR="${snakemake_output[batches]}"

BATCH_SIZE=25
OUTPUT_PREFIX="${OUTDIR}/batch"

HEADER="CHROM\tPOS\tREF\tALT\tSample\tGT\tPloidy\tGene\tAnnotation\tImpact\tHGVS.c\tHGVS.p\tAA_Change\tVariant_Class\tLOF\tNMD"

mkdir -p "$OUTDIR"

TMP_FILE=$(mktemp)
trap 'rm -f "$TMP_FILE" "${TMP_FILE}.samples"' EXIT

echo "VCF: $VCF_FILE"
echo "Output directory: $OUTDIR"
echo "Batch size: $BATCH_SIZE"

##############################################################################
# Step 1. Flatten the VCF into one row per (sample, allele, annotation)
##############################################################################

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

    if (field == "." || field == "")
        return;

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

        # Detect ploidy explicitly instead of assuming a separator is present.
        is_phased_or_diploid = (gt ~ /[\/|]/);
        n_alleles = split(gt, gts, /[\/|]/);

        if (n_alleles == 1)
            ploidy = "haploid";
        else if (n_alleles == 2)
            ploidy = "diploid";
        else
            ploidy = "polyploid_" n_alleles;

        n_nonref = 0;

        for (k in gts)
            if (gts[k] != "0" && gts[k] != ".")
                n_nonref++;

        if (n_nonref == 0)
            continue;

        delete seen;

        for (k in gts) {

            a = gts[k];

            if (a == "0" || a == "." || (a in seen))
                continue;

            seen[a] = 1;

            for (i = 1; i <= length(anns); i++) {

                split(anns[i], f, "|");

                ann_alt = f[1];
                annotation = f[2];
                impact = f[3];
                gene = f[4];
                hgvsc = f[10];
                hgvsp = f[11];

                if (ann_alt != alt_alleles[a])
                    continue;

                if (length(ref) < length(alt_alleles[a]))
                    variant_class="Insertion";
                else if (length(ref) > length(alt_alleles[a]))
                    variant_class="Deletion";
                else
                    variant_class="SNP";

                short = abbr(hgvsp);

                lof_flag = (gene in lof_genes) ? 1 : 0;
                nmd_flag = (gene in nmd_genes) ? 1 : 0;

                print chrom, pos, ref, alt_alleles[a], sample, gt, ploidy,
                      gene, annotation, impact,
                      hgvsc, hgvsp, short,
                      variant_class,
                      lof_flag,
                      nmd_flag;
            }
        }
    }
}
' | sort -k5,5 > "$TMP_FILE"

##############################################################################
# Step 2. Split into batches of samples
##############################################################################

echo "Splitting into sample batches..."

mapfile -t SAMPLES < <(cut -f5 "$TMP_FILE" | uniq)

echo "Detected ${#SAMPLES[@]} samples."

batch_num=1

for ((i=0; i<${#SAMPLES[@]}; i+=BATCH_SIZE)); do

    batch_file=$(printf "%s_%03d.tsv" "$OUTPUT_PREFIX" "$batch_num")

    batch_samples=("${SAMPLES[@]:i:BATCH_SIZE}")

    echo -e "$HEADER" > "$batch_file"

    printf '%s\n' "${batch_samples[@]}" > "${TMP_FILE}.samples"

    awk -F'\t' '
        NR==FNR {
            keep[$1]=1
            next
        }
        ($5 in keep)
    ' "${TMP_FILE}.samples" "$TMP_FILE" >> "$batch_file"

    echo "Wrote $(basename "$batch_file") (${#batch_samples[@]} samples)"

    ((batch_num++))

done

echo "Finished successfully."