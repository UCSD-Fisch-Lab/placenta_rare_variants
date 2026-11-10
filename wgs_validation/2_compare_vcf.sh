source /home/m7huang/miniconda3/etc/profile.d/conda.sh

conda activate bcftools

## HaplotypeCaller vs HaplotypeCaller
rna_dir="placenta_hc_rnaseq"
wgs_dir="placenta_hc_wgs"
print_dir="hc-rna_hc-wgs"

while IFS= read -r line; do
    sample=$(echo "$line" | awk '{print $1 "_" $2}')
    rna_file=$(echo "$line" | awk '{print $3}')
    wgs_file=$(echo "$line" | awk '{print $4}')

    echo "Sample: $sample"
    echo "RNA: ${rna_dir}/${rna_file}"
    echo "WGS: ${wgs_dir}/${wgs_file}"

    bcftools index "samples/${rna_dir}/${rna_file}"
    bcftools index "samples/${wgs_dir}/${wgs_file}"
    bcftools isec "samples/${rna_dir}/${rna_file}" "samples/${wgs_dir}/${wgs_file}" -p results/${print_dir}/${sample}_rna_wgs

done < $1

## HaplotypeCaller vs Mutect2
rna_dir="placenta_hc_rnaseq"
wgs_dir="placenta_mutect2_wgs"
print_dir="hc-rna_mutect2-wgs"

while IFS= read -r line; do
    sample=$(echo "$line" | awk '{print $1 "_" $2}')
    rna_file=$(echo "$line" | awk '{print $3}')
    wgs_file=$(echo "$line" | awk '{print $5}')

    echo "Sample: $sample"
    echo "RNA: ${rna_dir}/${rna_file}"
    echo "WGS: ${wgs_dir}/${wgs_file}"

    bcftools index "samples/${wgs_dir}/${wgs_file}"
    bcftools isec "samples/${rna_dir}/${rna_file}" "samples/${wgs_dir}/${wgs_file}" -p results/${print_dir}/${sample}_rna_wgs

done < $1

## Quality filter: HaplotypeCaller vs HaplotypeCaller
mkdir /expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/samples/placenta_quality_hc_wgs/
mkdir /expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/samples/placenta_quality_hc_rnaseq/

for sample in /expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/samples/placenta_hc_wgs/*.filt.vcf.gz; do
    base=$(basename "${sample/}" .vcf.gz);
    echo -e "$base";
    bcftools filter -i 'GQ>20 && QD>5 && MQ>40 && ReadPosRankSum>-3 && MQRankSum>-10 && FS<200 && QUAL>50' $sample -o /expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/samples/placenta_quality_hc_wgs/$base.vcf.gz;
done

for sample in /expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/samples/placenta_hc_rnaseq/*-gatk-haplotype-annotated.vcf.gz; do
    base=$(basename "${sample/}" .vcf.gz);
    echo -e "$base";
    bcftools filter -i 'GQ>20 && QD>5 && MQ>40 && ReadPosRankSum>-3 && MQRankSum>-10 && FS<200 && QUAL>50' $sample -o /expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/samples/placenta_quality_hc_rnaseq/$base.vcf.gz;
done < $1


rna_dir="placenta_quality_hc_rnaseq"
wgs_dir="placenta_quality_hc_wgs"
print_dir="placenta_quality_hc-rna_placenta_quality_hc-wgs"

while IFS= read -r line; do
    sample=$(echo "$line" | awk '{print $1 "_" $2}')
    rna_file=$(echo "$line" | awk '{print $3}')
    wgs_file=$(echo "$line" | awk '{print $4}')

    echo "Sample: $sample"
    echo "RNA: ${rna_dir}/${rna_file}"
    echo "WGS: ${wgs_dir}/${wgs_file}"

    bcftools index "samples/${rna_dir}/${rna_file}"
    bcftools index "samples/${wgs_dir}/${wgs_file}"
    bcftools view -f 'PASS' "samples/${rna_dir}/${rna_file}" -o "samples/${rna_dir}/${rna_file}.PASS.vcf.gz"
    bcftools index "samples/${rna_dir}/${rna_file}.PASS.vcf.gz"
    bcftools view -f 'PASS' "samples/${wgs_dir}/${wgs_file}" -o "samples/${wgs_dir}/${wgs_file}.PASS.vcf.gz"
    bcftools index "samples/${wgs_dir}/${wgs_file}.PASS.vcf.gz"

    bcftools isec "samples/${rna_dir}/${rna_file}.PASS.vcf.gz" "samples/${wgs_dir}/${wgs_file}.PASS.vcf.gz" -p results/${print_dir}/${sample}_filt_rna_filt_wgs

done < $1


for item in "hc-rna_hc-wgs" "hc-rna_mutect2-wgs" "placenta_quality_hc-rna_placenta_quality_hc-wgs"
do
    comp_path="/expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/results/${item}/"
    output="/expanse/lustre/projects/csd728/m7huang/wgs_rna_validation/results/results_table/${item}.txt"
    echo $comp_path
    
    printf "Folder\tRNA-seq\tWGS\tShared\tPerc_RNA\tPerc_WGS\n" > "$output"
    
    for folder in "$comp_path"/*/; do
        [ -d "$folder" ] || continue

        name=$(basename "$folder")

        rna=$(grep -vc '^#' "$folder/0000.vcf")
        wgs=$(grep -vc '^#' "$folder/0001.vcf")
        shared=$(grep -vc '^#' "$folder/0002.vcf")

        perc_rna=$(echo "scale=5; $shared * 100 / ($shared + $rna)" | bc -l)
        perc_wgs=$(echo "scale=5; $shared * 100 / ($shared + $wgs)" | bc -l)

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$name" "$rna" "$wgs" "$shared" "$perc_rna" "$perc_wgs" >> "$output"
    done
done


