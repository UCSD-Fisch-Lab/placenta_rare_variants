## samples that need combining
## GSM8817240: SRR32499798, SRR32499799, SRR32499800, SRR32499801, SRR32499802
## GSM8817241: SRR32499793, SRR32499794, SRR32499795, SRR32499796, SRR32499797
## GSM8817242: SRR32499788, SRR32499789, SRR32499790, SRR32499791, SRR32499792
## GSM8817243: SRR32499786, SRR32499814, SRR32499815, SRR32499816, SRR32499817

source /home/m7huang/miniconda3/etc/profile.d/conda.sh

FOLDER=$1
echo $FOLDER

cd /expanse/lustre/projects/csd728/m7huang/single_cell_paik/data/

## concatenation
#cat SRR32499798_4.fastq SRR32499799_4.fastq SRR32499800_4.fastq SRR32499801_4.fastq SRR32499802_4.fastq > GSM8817240_4.fastq
#cat SRR32499798_3.fastq SRR32499799_3.fastq SRR32499800_3.fastq SRR32499801_3.fastq SRR32499802_3.fastq > GSM8817240_3.fastq
#cat SRR32499793_4.fastq SRR32499794_4.fastq SRR32499795_4.fastq SRR32499796_4.fastq SRR32499797_4.fastq > GSM8817241_4.fastq
#cat SRR32499793_3.fastq SRR32499794_3.fastq SRR32499795_3.fastq SRR32499796_3.fastq SRR32499797_3.fastq > GSM8817241_3.fastq
#cat SRR32499788_4.fastq SRR32499789_4.fastq SRR32499790_4.fastq SRR32499791_4.fastq SRR32499792_4.fastq > GSM8817242_4.fastq
#cat SRR32499788_3.fastq SRR32499789_3.fastq SRR32499790_3.fastq SRR32499791_3.fastq SRR32499792_3.fastq > GSM8817242_3.fastq
#cat SRR32499786_4.fastq SRR32499814_4.fastq SRR32499815_4.fastq SRR32499816_4.fastq SRR32499817_4.fastq > GSM8817243_4.fastq
#cat SRR32499786_3.fastq SRR32499814_3.fastq SRR32499815_3.fastq SRR32499816_3.fastq SRR32499817_3.fastq > GSM8817243_3.fastq


############### Start: STARSolo Alignment ######################

# STARSolo: contains STAR conda installation, 16 CPUs
# Genome Alignment for 10X Chromium V3 3' Paired-End
# in: R1 FASTQ (with barcode), R2 FASTQ
# Chromium Single Cell 3′ Library Kit (10X Genomics)

conda activate STARSolo_2
SINGLE_CELL_cDNA="/expanse/lustre/projects/csd728/m7huang/single_cell_paik/data/${FOLDER}_4.fastq"
SINGLE_CELL_CB_UMI="/expanse/lustre/projects/csd728/m7huang/single_cell_paik/data/${FOLDER}_3.fastq"

## check the sequencing manifest for FASTQ number
#SINGLE_CELL_cDNA="/expanse/lustre/projects/csd728/m7huang/single_cell_paik/data/${FOLDER}_3.fastq"
#SINGLE_CELL_CB_UMI="/expanse/lustre/projects/csd728/m7huang/single_cell_paik/data/${FOLDER}_2.fastq"

STAR --runThreadN 16 --soloType CB_UMI_Simple --soloCBwhitelist /home/m7huang/software/3M-february-2018.txt --genomeDir /expanse/lustre/projects/csd728/m7huang/STAR_arcgenome --readFilesIn ${SINGLE_CELL_cDNA} ${SINGLE_CELL_CB_UMI} --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --limitOutSJcollapsed 5000000 --soloCellFilter EmptyDrops_CR --twopassMode Basic --limitSjdbInsertNsj 1500000 --soloCellReadStats Standard --outFileNamePrefix /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/$FOLDER/

################ End: STARSolo Alignment #######################

####### Start: Split bam by barcode for genome coverage ########
conda activate umi_tools
mkdir /expanse/lustre/projects/csd728/m7huang/single_cell_paik/cell_barcodes/${FOLDER}
mkdir /expanse/lustre/projects/csd728/m7huang/single_cell_paik/cell_barcodes/${FOLDER}/cells

cd /expanse/lustre/projects/csd728/m7huang/single_cell_paik/cell_barcodes/${FOLDER}/cells

python /home/m7huang/software/split_bam_by_barcode.py /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${FOLDER}/Aligned.sortedByCoord.out.bam /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${FOLDER}/Solo.out/Gene/filtered/barcodes.tsv /expanse/lustre/projects/csd728/m7huang/single_cell_paik/cell_barcodes/${FOLDER}/cells/ .
######## End: Split bam by barcode for genome coverage #########

###### Start: Monopogen (run w/128GB Memory, 22 threads) #######
path="/expanse/lustre/projects/csd728/m7huang/single_cell_afshar/Monopogen/"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

data_path="/expanse/lustre/projects/csd728/m7huang/single_cell_paik/monopogen_results"
cd $data_path

sample=$1
echo $sample

conda activate samtools
samtools index /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Aligned.sortedByCoord.out.bam

echo "$sample,/expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Aligned.sortedByCoord.out.bam" > ${data_path}/bam_list/${sample}.bam.lst

conda activate monopogen_2
python ${path}/src/Monopogen_2.py preProcess -b ${data_path}/bam_list/${sample}.bam.lst -o ${data_path}/results/bm_${sample} -a ${path}/apps -t 22
python ${path}/src/Monopogen_2.py germline  \
    -a ${path}/apps -r /expanse/lustre/projects/csd728/m7huang/single_cell_afshar/GRCh38_chr.region.lst \
    -p /expanse/lustre/projects/csd728/m7huang/single_cell_afshar/Monopogen/resource/20201028_3202_phased/ -t 22 -m 3 \
    -g /home/m7huang/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa -s all -o ${data_path}/results/bm_${sample}

awk -F '\t' 'NR==FNR { keep[$1]; next } BEGIN { OFS=","; print "cell","id" } FNR>1 && ($1 in keep) { print $1, $2 }' /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Solo.out/Gene/filtered/barcodes.tsv /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Solo.out/Gene/CellReads.stats > /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Solo.out/Gene/Filtered_CellReads.csv

python  ${path}/src/Monopogen_2.py  somatic  \
    -a  ${path}/apps  -r  /expanse/lustre/projects/csd728/m7huang/single_cell_afshar/GRCh38_chr.region.lst  -t 22 \
    -i  ${data_path}/results/bm_${sample}  -l  /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Solo.out/Gene/Filtered_CellReads.csv -s featureInfo     \
    -g  /home/m7huang/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa
python  ${path}/src/Monopogen_2.py  somatic  \
    -a  ${path}/apps  -r  /expanse/lustre/projects/csd728/m7huang/single_cell_afshar/GRCh38_chr.region.lst  -t 22 \
    -i  ${data_path}/results/bm_${sample}  -l  /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Solo.out/Gene/Filtered_CellReads.csv -s cellScan \
    -g  /home/m7huang/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa
python  ${path}/src/Monopogen_2.py  somatic  \
    -a  ${path}/apps  -r  /expanse/lustre/projects/csd728/m7huang/single_cell_afshar/GRCh38_chr.region.lst  -t 22 \
    -i  ${data_path}/results/bm_${sample}  -l  /expanse/lustre/projects/csd728/m7huang/single_cell_paik/samples_2/${sample}/Solo.out/Gene/Filtered_CellReads.csv -s LDrefinement \
    -g  /home/m7huang/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa

exit 0
