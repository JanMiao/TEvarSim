####################################################
### complete workflow for the short-read and long-read based tool benchmarking
####################################################

#####################################################
### install TEvarSim
conda create -n tevarsim -c bioconda mason pbsim3 python=3.8
conda activate tevarsim
pip install TEvarSim

### install other tools
conda install -c bioconda bwa samtools bcftools

# Dir=/lustre/scratch/jiamiao/tevarsim/
DIR=`pwd`
#cd $DIR
mkdir data short long

#####################################################
### download required files
cd data
# genome sequence of GRCh38
wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
# repeatmasker annotation for GRCh38
wget https://www.repeatmasker.org/genomes/hg38/rmsk4.0.5_rb20140131/hg38.fa.out.gz
gunzip hg38.fa.out.gz
# human TE library
wget https://raw.githubusercontent.com/JanMiao/TEvarSim/refs/heads/main/data/human_TE.fa

# split genome sequence into chromosomes
for chr in {1..22}
do
samtools faidx hg38.fa chr${chr} > chr${chr}.fa
done

#####################################################
### benchmark for short reads

# load tevarsim environment
conda activate tevarsim

# TErandom module and genome simulation
cd ${DIR}/short
for chr in {1..22}
do
tevarsim TErandom --consensus ../data/human_TE.fa --knownDEL ../data/hg38.fa.out --CHR chr${chr} --nTE 100 --ins-ratio 0.6 --seed 1 --outprefix Sim${chr}
tevarsim Simulate --ref ../data/chr${chr}.fa --bed Sim${chr}.bed --num 1 --pool Sim${chr}.fa --seed 1 --outprefix Sim_chr${chr}
done
# combine to a single genome
cat Sim_chr*fa > genome.fa
bcftools concat -Ov -o genome.vcf Sim_chr*vcf
# remove intermediate files
rm Sim_chr*fa Sim_chr*vcf Sim*.bed Sim*fa
# short read simulation
tevarsim Readsim --type short --genome genome.fa --depth 30
zcat chr*R1.fq.gz > read_R1.fq.gz
zcat chr*R2.fq.gz > read_R2.fq.gz
# remove intermediate files
rm chr*R1.fq.gz chr*R2.fq.gz

### METL
REF=${DIR}/data/hg38.fa
# map
bwa index ${REF}
bwa mem -t 16 ${REF} read_R1.fq.gz read_R2.fq.gz | samtools sort -o sorted.bam
samtools index sorted.bam

# MELT
MELT_DIR=${DIR}/MELTv2.2.2 # MELT_DIR should be modified according to your installation path
# MELT_DIR=/home/jiamiao/software/MELTv2.2.2
MELT_JAR=${MELT_DIR}/MELT.jar
# preprocess
java -jar ${MELT_JAR} Preprocess -h ${REF} -bamfile sorted.bam
# single sample
ls ${MELT_DIR}/me_refs/Hg38/*zip > meltTransposonFileList.txt 
java -jar ${MELT_JAR} Single \
-h ${REF} \
-bamfile sorted.bam \
-n ${MELT_DIR}/add_bed_files/Hg38/Hg38.genes.bed \
-c 30 -t meltTransposonFileList.txt -w melt


#####################################################
### benchmark for long reads
cd ${DIR}/long
tevarsim Readsim --type long --genome ../short/genome.fa --depth 30

minimap2 -aYx map-ont -t 16 ${REF} reads.fq.gz | samtools view -bhS - | samtools sort -o sorted.bam - 
samtools index sorted.bam

# xTEA for long reads
echo "sample" >> sample_id.txt
echo "sample ${DIR}/long/sorted.bam" >> long_read_bam_list.txt

# you may activate your conda environment by: conda activate xtea
xTEA_DIR=${DIR}/xTea # xTEA_DIR should be modified according to your installation path
xtea_long -i `pwd`/sample_id.txt -b `pwd`/long_read_bam_list.txt -p `pwd` \
-o submit_jobs.sh --rmsk ${xTEA_DIR}/lib/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out \
-r ${REF} --cns ${xTEA_DIR}/lib/consensus/LINE1.fa \
--rep ${xTEA_DIR}/lib/ --xtea ${xTEA_DIR}/xtea_long \
-f 31 -y 15 -n 10 -m 32 --slurm -t 2-0:0:0 # maybe -p should be used to specify your partition

cp sample/run_xTEA_pipeline.sh run_xTEA_pipeline.sh
# you may need to modify the run_xTEA_pipeline.sh script according to your environment
sbatch run_xTEA_pipeline.sh


### Accuracy evaluation
# metl
tevarsim Compare -T genome.vcf -P metl/ALU.final_comp.vcf -e Alu --INSonly --nHap 1 --predType VCF -O metl_alu --truthID Hap1 --predID sorted
tevarsim Compare -T genome.vcf -P metl/LINE1.final_comp.vcf -e L1 --INSonly --nHap 1 --predType VCF -O metl_L1 --truthID Hap1 --predID sorted
tevarsim Compare -T genome.vcf -P metl/SVA.final_comp.vcf -e SVA --INSonly --nHap 1 --predType VCF -O metl_SVA --truthID Hap1 --predID sorted
tevarsim Compare -T genome.vcf -P metl/HERVK.final_comp.vcf -e ERVK --INSonly --nHap 1 --predType VCF -O metl_HERVK --truthID Hap1 --predID sorted

# xtea
awk '{print $1, $2, $2+1, 0}' xTea/classified_results.txt.merged_ALU.txt | sort -k1,1 -k2,2n > xtea_alu.bed
awk '{print $1, $2, $2+1, 0}' xTea/classified_results.txt.merged_LINE1.txt | sort -k1,1 -k2,2n > xtea_l1.bed
awk '{print $1, $2, $2+1, 0}' xTea/classified_results.txt.merged_SVA.txt | sort -k1,1 -k2,2n > xtea_sva.bed
awk '{print $1, $2, $2+1, 0}' xTea/classified_results.txt.merged_HERV.txt | sort -k1,1 -k2,2n > xtea_herv.bed
tevarsim Compare -T genome.vcf -P xtea_alu.bed -e Alu --INSonly --nHap 1 --predType BED -O xtea_alu --truthID Hap1 --predID sorted
tevarsim Compare -T genome.vcf -P xtea_l1.bed -e L1 --INSonly --nHap 1 --predType BED -O xtea_l1 --truthID Hap1 --predID sorted
tevarsim Compare -T genome.vcf -P xtea_sva.bed -e SVA --INSonly --nHap 1 --predType BED -O xtea_sva --truthID Hap1 --predID sorted
tevarsim Compare -T genome.vcf -P xtea_herv.bed -e ERVK --INSonly --nHap 1 --predType BED -O xtea_herv --truthID Hap1 --predID sorted
