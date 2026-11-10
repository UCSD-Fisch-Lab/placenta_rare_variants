set -euo pipefail

source /home/m7huang/miniconda3/etc/profile.d/conda.sh

##################### UPSTREAM ####################
## wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
## wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
#
##path=/expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/
#conda activate gatk
##gatk SelectVariants \
##        -V $path/tools/af-only-gnomad.hg38.vcf.gz \
##        -select-type SNP -restrict-alleles-to BIALLELIC \
##        -select "AF > 0.05" \
##        -O $path/tools/af-only-gnomad-common-biallelic.grch38.main.vcf.gz \
##        --lenient
#
##bcftools filter -e 'INFO/AF > 0.01' \
##    -o af-only-gnomad-common-0.01.grch38.main.vcf.gz \
##    -O z \
##    af-only-gnomad.hg38.vcf.gz

#conda activate gatk
#bwa index /expanse/lustre/projects/csd728/m7huang/raw_data/hg38.fa

#################### DNA-Seq Alignment Command Line Parameters ####################
### from: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
### https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
## ensure that the average read length > 70 bps
# samtools stats $input/${sample}_R1.fastq.gz.bam  | grep "average length"

#input="/expanse/lustre/projects/csd728/m7huang/raw_data/msc_wgs/merged_files"
#output="/expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/samples_20251105"

sample=$1

echo $sample
cd $output

conda activate samtools
samtools stats $input/${sample}_L001_R1_001.fastq.gz  | grep "average length"

mkdir $sample
mkdir $sample/tmp

## run with 8 threads!
conda activate gatk
bwa mem \
    -t 8 \
    -T 0 \
    /expanse/lustre/projects/csd728/m7huang/raw_data/hg38.fa \
    $input/${sample}_R1_001.fastq.gz $input/${sample}_R2_001.fastq.gz > $output/$sample/$sample.sam

conda activate samtools
samtools view $output/$sample/$sample.sam \
    -Shb \
    -o $output/$sample/$sample.bam

 rm $output/$sample/$sample.sam

###################### WGS Co-Cleaning Command Line Parameters ####################
conda activate gatk
gatk SortSam \
    CREATE_INDEX=true \
    INPUT=$input/$sample-ready.bam \
    OUTPUT=$output/$sample/$sample.sort.bam \
    SORT_ORDER=coordinate \
    TMP_DIR=$output/$sample/tmp \
    VALIDATION_STRINGENCY=STRICT

gatk --java-options "-Xmx8G" MarkDuplicates \
    CREATE_INDEX=true \
    INPUT=$output/$sample/$sample.sort.bam \
    O=$output/$sample/$sample.sort.marked_duplicates.bam \
    M=$output/$sample/$sample.marked_dup_metrics.txt \
    TMP_DIR=$output/$sample/tmp \
    VALIDATION_STRINGENCY=STRICT

#### add random read group ####
gatk --java-options "-Xmx8G" AddOrReplaceReadGroups \
   I=$output/$sample/$sample.sort.marked_duplicates.bam \
   O=$output/$sample/$sample.sort.marked_duplicates.rg.bam \
   TMP_DIR=$output/$sample/tmp \
   RGID=1 \
   RGLB=1 \
   RGPL=1 \
   RGPU=1 \
   RGSM=1

gatk --java-options "-Xmx8G" BaseRecalibrator \
   --input $output/$sample/$sample.sort.marked_duplicates.rg.bam \
   -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
   --tmp-dir $output/$sample/tmp \
   --known-sites /expanse/lustre/projects/csd728/m7huang/iPSC_cell_line_whole_genome/tools/1000G_phase3_v4_20130502.sites.hg38.vcf \
   -O $output/$sample/$sample.recal_data.table

gatk --java-options "-Xmx8G" ApplyBQSR \
   -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
   -I $output/$sample/$sample.sort.marked_duplicates.rg.bam \
   --tmp-dir $output/$sample/tmp \
   --bqsr-recal-file $output/$sample/$sample.recal_data.table \
 -O $output/$sample/$sample.sort.marked_duplicates.bqsr.bam

 ################### HaplotypeCaller Call Command-Line Parameters ####################
conda activate gatk

gatk --java-options "-Xmx8g" HaplotypeCaller  \
    -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa\
    -I $output/$sample/$sample.sort.marked_duplicates.bqsr.bam \
    -O $output/$sample/$sample.hc2.vcf.gz \
    --native-pair-hmm-threads 4

gatk VariantFiltration \
   -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
   -V $output/$sample/$sample.hc2.vcf.gz \
   -O $output/$sample/$sample.hc2.filt.vcf.gz \
    --filter-expression "FS > 30.0" --filter-expression "QD < 2.0" \
    --filter-name FS --filter-name QD \
    --cluster-size 3 --cluster-window-size 35

conda activate bcftools
bcftools filter -i 'INFO/DP>10 && GQ>20 && MQ>40 && ReadPosRankSum>-3 && MQRankSum>-10 && FS<200 && QUAL>50' -o $output/$sample/$sample.notGVCF.filt.final.vcf.gz $output/$sample/$sample.hc2.filt.vcf.gz

#################### Tumor-Only Variant Call Command-Line Parameters ####################
### updated workflow from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
### https://gatk.broadinstitute.org/hc/en-us/community/posts/360067310872-How-to-find-or-generate-common-germline-variant-sites-VCF-required-by-GetPileupSummaries
### wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz

conda activate gatk
gatk --java-options "-Xmx10G" Mutect2 \
    -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa\
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
    -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    -I $input/$sample/$sample.sort.marked_duplicates.bqsr.bam \
    -O $output/$sample/$sample.mutect2.vcf \
    --native-pair-hmm-threads 4 \
    --tmp-dir $output/$sample/tmp \
    --f1r2-tar-gz $output/$sample/$sample.f1r2.tar.gz \
    --germline-resource /expanse/lustre/projects/csd728/m7huang/iPSC_cell_line_whole_genome/tools/af-only-gnomad.hg38.vcf.gz \
    -pon /expanse/lustre/projects/csd728/m7huang/iPSC_cell_line_whole_genome/tools/1000g_pon.hg38.vcf.gz

gatk --java-options "-Xmx10G" LearnReadOrientationModel \
    -I $output/$sample/$sample.f1r2.tar.gz \
    -O $output/$sample/$sample.read.orientation.model.tar.gz \
    --tmp-dir $output/$sample/tmp

gatk --java-options "-Xmx10G" GetPileupSummaries \
    -I $output/$sample/$sample.sort.marked_duplicates.bqsr.bam \
    -O $output/$sample/$sample.targeted_sequencing.table \
    -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
    --tmp-dir $output/$sample/ \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
    -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM \
    -V /expanse/lustre/projects/csd728/m7huang/iPSC_cell_line_whole_genome/tools/af-only-gnomad-common-biallelic.grch38.main.vcf.gz

gatk --java-options "-Xmx10G" CalculateContamination \
    --tmp-dir $output/$sample/ \
    -I $output/$sample/$sample.targeted_sequencing.table \
    -O $output/$sample/$sample.targeted_sequencing.contamination.table

gatk --java-options "-Xmx10G" FilterMutectCalls \
    -O $output/$sample/$sample.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \
    -R /expanse/lustre/projects/csd728/kfisch/raw_data/hg38.fa \
    -V $output/$sample/$sample.mutect2.vcf \
    --contamination-table $output/$sample/$sample.targeted_sequencing.contamination.table \
    --tmp-dir $output/$sample/tmp \
    --ob-priors $output/$sample/$sample.read.orientation.model.tar.gz \
    -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
    -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
    -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM

echo $sample > $output/$sample/$sample.sample_name.txt

conda activate bcftools
bcftools reheader -s $output/$sample/$sample.sample_name.txt $output/$sample/$sample.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz > $output/$sample/$sample.vcf.gz

final_output="/expanse/lustre/projects/csd728/m7huang/cell_line_wgs_somatic/hc_placenta_20261005/"

conda activate bcftools
bcftools filter -e '(((REF="C" && ALT="A") || (REF="G" && ALT="T")) && ((INFO/AS_UNIQ_ALT_READ_COUNT[0] < 3) || (FORMAT/AD[0:1] < 4) || (FORMAT/F1R2[0:1] == 0 || FORMAT/F2R1[0:1] == 0) || (FORMAT/SB[0:2] == 0 || FORMAT/SB[0:3] == 0) || (INFO/TLOD[0] < 8) || (INFO/MPOS[0] < 8) || (INFO/MBQ[1] < 25) || (INFO/MMQ[1] < 40) || (INFO/ROQ && INFO/ROQ < 30)))' \
  -Oz \
  -o "${final_output}/final_filtered/${sample}.filtered.vcf.gz" \
  "${output}/${sample}/${sample}.vcf.gz"

bcftools view -f 'PASS,.' -Oz \
    -o "${final_output}/final_filtered/${sample}.filtered.pass.vcf.gz" \
    "${final_output}/final_filtered/${sample}.filtered.vcf.gz"

exit 0
