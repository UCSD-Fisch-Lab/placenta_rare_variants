source /home/m7huang/miniconda3/etc/profile.d/conda.sh

### https://github.com/claupomm/RNA-seq_snv_tumour_only?tab=readme-ov-file ###

########## START: placental-only upstream ###########
## RNA editing:
# cd /expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/tools
# wget -O - http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz | gunzip > TABLE1_hg38_v3.txt
# awk 'NR>1 {print $2"\t"$3}' TABLE1_hg38_v3.txt > TABLE1_hg38_v3.positions.txt

## Processing of placental RNA editing:
#cd /expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/tools
#awk '{
#    # remove leading prefix
#    sub(/^no_barcode_/, "")
#
#    # split remaining string by "_"
#    n = split($0, a, "_")
#
#    # find chromosome (everything before first numeric field)
#    chrom = ""
#    for (i=1; i<=n; i++) {
#        if (a[i] ~ /^[0-9]+$/) {   # position found
#            pos = a[i]
#            ref = a[i+1]
#            alt = a[i+2]
#            strand = a[i+3]
#            break
#        }
#        chrom = (chrom=="" ? a[i] : chrom "_" a[i])
#    }
#
#    # allele format: ref+alt+strand → e.g. "AG+"
#    allele = ref alt strand
#
#    # keep only A→G (+) or T→C (−)
#    if (allele == "AG+" || allele == "TC-") {
#        print chrom "\t" pos "\t" allele
#    }
#}' /expanse/lustre/projects/csd728/cnmason/filtered_sites_0.25.txt \
#> rna_edit_sites_0.25_ai_20251218.txt
#
#awk '{print $1"\t"$2}' rna_edit_sites_0.25_ai_20251218.txt > rna_edit_sites_0.25_ai_positions_20251218.txt

## Low complexity regions:
#wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
#unzip snpEff_latest_core.zip
#cd snpEff/
#curl -o LCR-hs38_with_chr.bed.gz -L https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
#gunzip -c LCR-hs38_with_chr.bed.gz > LCR-hs38_with_chr.bed
#sed 's/^chr//g' LCR-hs38_with_chr.bed > LCR-hs38.bed

## 1000 Genomes:
# wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase3_v4_20130502.sites.hg38.vcf
#bcftools filter -i 'INFO/AF>0.1' -o /expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/tools/1000G_phase3_v4_20130502.sites.hg38.af01.vcf /expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/tools/1000G_phase3_v4_20130502.sites.hg38.vcf

## dbSNP:
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
# zgrep "^#" GCF_000001405.40.gz > GCF_000001405.40.common.vcf
# zgrep ";COMMON" GCF_000001405.40.gz >> GCF_000001405.40.common.vcf
# refseq=$(cut -f1 gcf_refseq2chr.tsv)
# IFS=' ' read -ra refseq <<< $(echo $refseq)
# chr=$(cut -f2 gcf_refseq2chr.tsv)
# IFS=$' ' read -ra chr <<< $(echo $chr)
# touch GCF_000001405.40.common.chr.tsv
#for i in ${!refseq[@]}; do
#date
#echo ${chr[$i]}
#zgrep ${refseq[$i]} GCF_000001405.40.common.vcf.gz | sed "s/${refseq[$i]}/${chr[$i]}/g" | cut -f1,2 >> GCF_000001405.40.common.chr.tsv
#done
########### END: Placental-only upstream ############

### START: Upstream VCF quality filtering ###
conda activate bcftools
output="/expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtering_20250528/all_samples/circulating_20251208/individual_vcfs/"
input="/expanse/projects/qstore/csd728/cayla/2025-04_Moufarrej_rnavars/final/2025*"

# run bcftools quality filtering
for sample in "$input"/*.gz; do
    base="$(basename "$sample" | sed 's/-.*//')"

    bcftools filter -i 'INFO/DP>10 && GQ>20 && MQ>40 && ReadPosRankSum>-3 && MQRankSum>-10 && FS<200 && QUAL>50 && possible_rnaedit!=1' -o $output/$base.notGVCF.vcf $sample
done
### END: Upstream VCF quality filtering ###

### START: Placental-only filtering ###
path=/expanse/lustre/projects/csd728/m7huang/somatic_filtering_20250528/
cd $path
mkdir samples

# directory with samples
ls /expanse/lustre/projects/csd728/kfisch/placental_variants/new_annotations/filtered_vcfs > vcf_names.txt

# path to tool folder

rna_edit=$path/tools/TABLE1_hg38_v3.txt.gz
SNPEFF=$path/tools/snpEff
filter=$path/tools/1000G_phase3_v4_20130502.sites.hg38.af01.vcf
filter2=$path/tools/GCF_000001405.40.common.chr.tsv

while IFS= read -r sample;
do
    echo "$sample"
    conda activate vcftools
    gzip -c samples/$sample.notGVCF.vcf > samples/$sample.notGVCF.vcf.gz

    java -jar $SNPEFF/SnpSift.jar intervals \
        -noLog -x -i samples/$sample.notGVCF.vcf.gz \
        $SNPEFF/LCR-hs38.bed | gzip > samples/$sample.pass2.lcr.vcf.gz

    conda activate vcftools
    vcftools --gzvcf samples/$sample.pass2.lcr.vcf.gz \
     --exclude-positions $filter --recode --recode-INFO-all --stdout | pigz -p4 > samples/$sample.hc.pass2.lcr.1kG.vcf.gz

    vcftools --gzvcf samples/$sample.hc.pass2.lcr.1kG.vcf.gz \
        --exclude-positions $filter2 --recode --recode-INFO-all --stdout | pigz -p4 > filtered_samples/$sample.hc.pass2.lcr.1kG.dbsnp.vcf.gz

done < vcf_names.txt
### END: Placental-only filtering ###

### START: Filtering of RNA edits ###
folder="Discovery"
#folder="Validation"
#folder="Circulating"

input="/expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/$folder/RNAseq_somatic_filtered_${folder}"
output="/expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/$folder/"

cd $output
mkdir rna_edit_filtration
mkdir rna_edits
mkdir final_vcfs

for sample in "$input"/*.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.gz; do
    base=$(basename "$sample" | cut -d'.' -f1)
    echo -e "$base"
    
    ## overlap of sample and rna editing sites
    conda activate vcftools
    vcftools --gzvcf $sample \
        --recode --positions /expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/tools/TABLE1_hg38_v3.positions.txt \
        --stdout > "$output/rna_edit_filtration/${base}_redi_test_intersect.vcf"
    vcftools --gzvcf $sample \
        --recode --positions /expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/tools/rna_edit_sites_0.25_ai_positions_20251218.txt \
        --stdout > "$output/rna_edit_filtration/${base}_pr_test_intersect.vcf"
    vcf-concat "$output/rna_edit_filtration/${base}_redi_test_intersect.vcf" "$output/rna_edit_filtration/${base}_pr_test_intersect.vcf" > "$output/rna_edit_filtration/${base}_test_intersect.vcf"
    
    ## filter for only A to G or T to C sites
    conda activate bcftools
    bcftools sort -Ou "$output/rna_edit_filtration/${base}_test_intersect.vcf" \
        | bcftools norm -d all -Ov -o "$output/rna_edit_filtration/${base}_test_intersect_sorted.vcf"
    bcftools view \
        -i '(REF="A" && ALT="G") || (REF="T" && ALT="C")' \
        "$output/rna_edit_filtration/${base}_test_intersect_sorted.vcf" | bgzip -c > "$output/rna_edits/${base}_rna_edits.vcf.bgz"
    bcftools index "$output/rna_edits/${base}_rna_edits.vcf.bgz"
    
    ## remove rna editing sites
    gunzip -c "$sample" \
      | bgzip -c \
      > "$output/rna_edit_filtration/${base}.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.bgz"
    bcftools index "$output/rna_edit_filtration/${base}.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.bgz"
    
    bcftools isec -C \
      -w1 \
      -Ov \
      "$output/rna_edit_filtration/${base}.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.bgz" \
      "$output/rna_edits/${base}_rna_edits.vcf.bgz" \
      -o "$output/final_vcfs/${base}_final.vcf"

    ## undo bgzip
    bgzip -d "$output/rna_edits/${base}_rna_edits.vcf.bgz"
done

### counting for validation
##gunzip -c /expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/Discovery/RNAseq_somatic_filtered_Discovery/SampleP997.notGVCF.vcf.hc.pass2.lcr.1kG.dbsnp.vcf.gz | grep -v '^#' | wc -l
##grep -v '^#' "/expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/Discovery/rna_edits/SampleP997_rna_edits.vcf" | wc -l
##grep -v '^#' "/expanse/lustre/projects/csd728/m7huang/placental_genomics_RNAseq/somatic_filtered_20251216/Discovery/final_vcfs/SampleP997_final.vcf" | wc -l
### END: Filtering of RNA edits ###

####START: Variant Annotation #####
#### needs four cpu cores ####

mkdir annotated_samples
mkdir maf_samples

conda activate vep

for sample in "$output"/final_vcfs/*.vcf; do
    base=$(basename "${sample/_all_filters/}" .vcf)
    echo -e "$base"
    
    /home/m7huang/miniconda3/envs/vep/bin/perl /home/m7huang/software/ensembl-vep/vep \
        --input_file "$output/final_vcfs/${base}_all_filters.vcf" \
        --output_file "$output/annotated_samples/${base}.vep.vcf" \
        --offline --dir /expanse/lustre/projects/csd728/kfisch/placental_variants/ \
        --fasta /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
        --force_overwrite --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf  --pubmed --fork 4 --cache_version 107 --polyphen b --everything --af_gnomadg --af_gnomade

    perl /home/m7huang/software/vcf2maf-1.6.21/vcf2maf.pl \
        --input-vcf "$output/annotated_samples/${base}.vep.vcf" \
        --output-maf "$output/maf_samples/${base}.maf" \
        --ref-fasta /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
         --tumor-id $base --inhibit-vep --ncbi-build GRCh38 --retain-info AC,AF,AN,BaseQRankSum,ClippingRankSum,DP,ExcessHet,FS,MLEAC,MLEAF,MQ,MQ0,MQRankSum,QD,ReadPosRankSum,SOR --retain-fmt GT,AD,DP,GQ,PL --retain-ann QUAL,FILTER,MAX_AF,MAX_AF_POPS,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_OTH_AF,gnomADg_SAS_AF,LoF,LoF_filter,LoF_flags,LoF_info,MANE_SELECT,MANE_PLUS_CLINICAL,APPRIS,UNIPROT_ISOFORM,miRNAdone
exit 0
