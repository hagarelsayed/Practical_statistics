


### Images
<img scr="./Images/Screenshot%20(1177).png" width = "100" >  
<img scr="./Images/Screenshot (1178).png" height = "100" >  
<img scr="./Images/Screenshot (1179).png" >  
 #You can also specify height but only one of them to maintain the aspect ratio. 

<img scr="Lesson3_Descriptiive Part2/Images/Screenshot (1189).png" >  

![alt text](./Images/Screenshot%20(1177).png)

![Alt text](relative/path/to/img.jpg?raw=true "Title")
![Alt text](./Images/Screenshot (1179).png?raw=true "Title")



That's right - the units of the variance are the square of the original units of your data.

## Standard Deviation

Important Final Points
The variance is used to compare the spread of two different groups. A set of data with higher variance is more spread out than a dataset with lower variance. Be careful though, there might just be an outlier (or outliers) that is increasing the variance when most of the data are actually very close.

When comparing the spread between two datasets, the units of each must be the same.

When data are related to money or the economy, higher variance (or standard deviation) is associated with higher risk.

The standard deviation is used more often in practice than the variance because it shares the units of the original dataset.

Use in the World
The standard deviation is associated with risk in finance, assists in determining the significance of drugs in medical studies, and measures the error of our results for predicting anything from the amount of rainfall we can expect tomorrow to your predicted commute time tomorrow.

These applications are beyond the scope of this lesson as they pertain to specific fields, but know that understanding the spread of a particular set of data is extremely important to many areas. In this lesson, you mastered the calculation of the most common measures of spread.

#### Questions:
1. For each of the below: If the statement is true, mark the box next to the statement.

If two datasets have the same variance, they will also have the same standard deviation.

That is correct! Besides the mean return of an investment, we should also consider the spread associated with the return. But just because the standard deviation associated with each investment is the same, this does not mean the max you could make for each investment is the same.

2. That's right! Because the return is the same year over year for Investment 1, it has 'no spread' or a standard deviation of 0. This smaller standard deviation is associated with smaller risk. Understanding the spread of values we could earn is just as important as understanding the expected return (mean return).

#### Useful Insight
The above example is a simplified version of the real world but does point out something useful that you may have heard before. Notice if you were not fully invested in either Investment 1 or fully invested in Investment 2, but instead, you were diversified across both investment options, you could earn more than either investment individually. This is the benefit of diversifying your portfolio for long-term gains. For short-term gains, you might not need or want to diversify. You could get lucky and hit short-term gains associated with the upswings (12%, 10%, or 7%) of Investment 2. However, you might also get unlucky, and hit a down term and earn nothing or even lose money on your investment using this same strategy.


```s


# downloading human whole exome ref
#wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

#downloading the WES data with ~40X coverage
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/052/SRR16939752/SRR16939752_1.fastq.gz -o $fastqs/SRR16939752_WES_of_Jianwei_Dong_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/052/SRR16939752/SRR16939752_2.fastq.gz -o $fastqs/SRR16939752_WES_of_Jianwei_Dong_2.fastq.gz

#saving paths for simpler code
ref="./ref/Homo_sapiens.GRCh38.cds.all.fa"
fastqs="./process/fastqs"
bam="./process/bams"
sam="./process/sams"
qc="./process/qc/stats"
sample="SRR16939752_WES_of_Jianwei_Dong"
samples="./samples.txt"
downsampling="./downsampling.txt" ##This file is manually filled by the % of reads needed for the final fastq file to be aligned (as the cov. of the original fastq is 48.9x, the percentages corresponding to 2x, 5x, 10x, 20x, were calculated to be 0.0408,0.102,0.204,0.4)
depth_file="./process/qc/stats/samples_depth.txt"

#bwa index $ref
##Algin original sample and produce depth to calculate % of read subsampling
#bwa mem $ref $fastqs/$sample'_1'.fastq.gz $fastqs/$sample'_2'.fastq.gz > $sam/$sample.sam
#samtools view -hbo $bam/$sample.bam $sam/$sample.sam
#samtools sort $bam/$sample.bam -o $bam/$sample.sorted.bam
#samtools index $bam/$sample.sorted.bam
samtools depth $bam/$sample.sorted.bam | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' >> $qc/$sample.depth.txt
##Create subsampled fastq files
x=$(cat $downsampling)

for fraction in $x
do
seqtk sample -s100 $fastqs/$sample'_1.fastq.gz' ${fraction} > $fastqs/$sample'_'${fraction}'_1.fastq'
seqtk sample -s100 $fastqs/$sample'_2.fastq.gz' ${fraction} > $fastqs/$sample'_'${fraction}'_2.fastq'
gzip $fastqs/$sample'_'${fraction}'_1.fastq'
gzip $fastqs/$sample'_'${fraction}'_2.fastq'
done
##Align subsampled fastq files and produce depth stats
for fraction in $x
do
bwa mem $ref $fastqs/$sample'_'${fraction}'_1.fastq.gz' $fastqs/$sample'_'${fraction}'_2.fastq.gz' > $sam/$sample'_'${fraction}.sam
samtools view -hbo $bam/$sample'_'${fraction}.bam $sam/$sample'_'${fraction}.sam
samtools sort $bam/$sample'_'${fraction}.bam -o $bams/$sample'_'${fraction}.sorted.bam
samtools index $bam/$sample'_'${fraction}.sorted.bam
#samtools depth $bam/$sample'_'${fraction}.sorted.bam >> $qc/$sample'_'${fraction}.depth.txt
done

##Produce the final sample list for stats loop
cd $bam
ls *.sorted.bam >> /media/mahmoud/storage/Backup/NU_Diploma/NGS2/Project2/sample_list.txt
cd ..
cd ..
pwd
s=$(cat sample_list.txt)
##Calculate depth comparison file 
for y in $s
do
echo "*******************************************" >> $qc/depth.txt
echo ${y} >> $qc/depth.txt
echo "*******************************************" >> $qc/depth.txt
samtools depth $bam/${y} | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' >> $qc/overall_depth.txt
echo "*******************************************" >> $qc/depth.txt
echo ${y} >> $qc/depth.txt
echo "*******************************************" >> overall_QC.txt
samtools flagstat $bam/${y} >> overall_QC.txt
done
##depth.txt is the final output for comparison of depth produced after subsampling
```
######################################################################

```s
############################
#BCFtools Variant calling###
############################

for f in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*.dedup.bam
do 
	l=$(basename $f | cut -d"." -f1)
	r=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/GCA_000001405.28_GRCh38.p13_genomic.adj.fna
	bcftools mpileup -Ou -f $r $f |bcftools call -Ov -mv > /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/bcftools/$l.bcftools.vcf
done

# GCA_000001405.28_GRCh38.p13_genomic.fna should be indexed using samtools faidx

######################################################################################3



#################################
#Imp note      #################
################################

#notes about using the reference: after many trials and errors, having a hard time of not finding any kind of variants with variant caller. we found that the name of chromosemes in the reference file GCA_000001405.28_GRCh38.p13_genomic.fna is "chromosome" and the name of the name of chromosomes in all other file is "chr". may be due to using different version of the index files and the reference file. so what will we do now is changing of the chromosome name in the reference file 

gunzip -c GCA_000001405.28_GRCh38.p13_genomic.fna.gz |awk -F '[ ,]' '{if ($1 ~ ">" ) print ">chr"$5; else print $0 }' > GCA_000001405.28_GRCh38.p13_genomic.adj.fna
#then we need to reindex the file again with samtools faidx and generate the sdf file for vcfeval tool again


#####################################################################################################33
#GATK Variant calling
######################


#the files of the sumpsampling are already aliged and sorted. so we will add read group to them.

#1-Add ReadGroups for each file

##############################
 picard_path=$CONDA_PREFIX/share/picard-*
for b in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*sorted.bam
do
	echo $b
	l=$(basename $b | cut -d"_" -f1,2)
	echo $l
	java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=$b O=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/$l.sorted.RG.bam RGID=1 RGSM=$l RGLB=lib1 RGPL=ILLUMINA RGPU=unit1

done


###########################
#2-mapping QC
# we did this step before in the 2nd group activty

##############################
#3-Mark duplicate using MarkDuplicates in Picard tools
for sample in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*.sorted.RG.bam
do
  l=$(basename $sample | cut -d"." -f1,2,3)
   java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/$l.dedup.bam METRICS_FILE=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/$l.metrics.txt;
done

#############################
#4-QC for errors in BAM files by ValidateSamFile

for s in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*.dedup.bam
do
java -Xmx2g -jar $picard_path/picard.jar ValidateSamFile I=$s MODE=SUMMARY

done

#No errors found
##########################################################3
#indexing#
# 1-samples#
##########

picard_path=$CONDA_PREFIX/share/picard-*
for sample in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*.dedup.bam
do
	java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$sample
done

# 2- Reference#
############
#we used the avilable index of the whole genome HG38
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz


gunzip -k GCA_000001405.28_GRCh38.p13_genomic.fna.gz  
samtools faidx GCA_000001405.28_GRCh38.p13_genomic.adj.fna # it will be used in bcftools and GATK

#Most GATK tools additionally require that the main FASTA file be accompanied by a dictionary file ending in .dict and an index file ending in .fai, because it allows efficient random access to the reference bases.

#Creates a sequence dictionary for a reference sequence.
gatk CreateSequenceDictionary R=GCA_000001405.28_GRCh38.p13_genomic.adj.fna O=GCA_000001405.28_GRCh38.p13_genomic.adj.dict

#############################
#5-Download known varinats
#wget 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/'

#for the known variant we will use HG002_NA24385_son/latest/GRCh38, the last modified was in 2020

wget 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_all.vcf.gz'
wget 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_all.vcf.gz.tbi'


#############################

#6-Recalibrate Bases BQSR
#under reconstruction...........


#indexing of the known variant

sed -i '6783465d;4878083d;4847292d;4205570d;1243522d' HG002_GRCh38_1_22_v4.2.1_all.vcf  # removing malformed records for indexing
gatk IndexFeatureFile -I HG002_GRCh38_1_22_v4.2.1_all.vcf


for sample in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*.dedup.bam
do
	r=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/GCA_000001405.28_GRCh38.p13_genomic.adj.fna
	v=/home/ngs/Documents/NGS2/group.activity/3rd_group_activity/goldstandard/HG002_GRCh38_1_22_v4.2.1_all.vcf

	n=$(basename $sample| cut -d"_" -f1)

	gatk --java-options "-Xmx2G" BaseRecalibrator -R $r -I $sample --known-sites $v -O/ home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.report


	gatk --java-options "-Xmx2G" ApplyBQSR -R $r -I $sample -bqsr /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.report -O /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.bam --add-output-sam-program-record --emit-original-quals

	gatk --java-options "-Xmx2G" BaseRecalibrator -R $r -I /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.bam --known-sites $v -O /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.report2

	gatk AnalyzeCovariates -before /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.report -after /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.report2 -plots /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/$n.bqsr.pdf
done



#A USER ERROR has occurred: Illegal argument value: Positional arguments were provided ',home/ngs/Documents/NGS2/group.activity/3rd_group_activity/BQSR/trial.bqsr.report}' but no positional argument is defined for this tool.

#we tried many times and we couldn't do BQSR due to the previus error. sowe just escaped it. 



####################################################################

#7-variant calling using gatk HaplotypeCaller.


for sample in /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/sunsampling-align/*.RG.dedup.bam
do
	
	r=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/GCA_000001405.28_GRCh38.p13_genomic.adj.fna


	n=$(basename $sample| cut -d"_" -f1)
	echo $sample

gatk --java-options "-Xmx2G" HaplotypeCaller \
-R $r \
-I $sample \
--pcr-indel-model NONE \
-O /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/gatk-vcf/$n.gatk.vcf.gz


done
 # I did this step for just 3 bam file due to lack of computational resources 2X,5X and 10X

###########################################################################
#8-VCFeval


v=/home/ngs/Documents/NGS2/group.activity/3rd_group_activity/goldstandard/HG002_GRCh38_1_22_v4.2.1_all.vcf.gz
r=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/GCA_000001405.28_GRCh38.p13_genomic.fna

#creat SDF file for the known variants


rtg format --output /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/GCA_000001405.28_GRCh38.p13_genomic.SDF $r

s=/home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/000001405.28_GRCh38.p13_genomic.SDF 
##############################################3
#indexing#
#########
rtg index $v # indexing of the known variant

for sample in /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/bcftools/*ools.vcf
do

	rtg bgzip $sample
	rtg index $sample.gz
done


#################################
#for gatk vcf files
for sample in /home/ngs/Documents/NGS2/group.activity/3nd_Group_activity/gatk-vcf/*.vcf.gz 


	n=$(basename $sample| cut -d"_" -f1)
	echo $sample

 
rtg vcfeval -b $v -c $sample -t $s -o /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/vcfeval/eval.$n

echo /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/vcfeval/eval_$n

done

##############################

#for bcftools vcf files
for sample in /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/bcftools/*ools.vcf.gz
do



	n=$(basename $sample| cut -d"_" -f1)
	echo $sample

 
rtg vcfeval -b $v -c $sample -t $s -o /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/vcfeval/eval.$n

echo /home/ngs/Documents/NGS2/group.activity/3rd_group_activity/vcfeval/eval_$n

done



#Error: Duplicate sequence names detected in SDF: /home/ngs/Documents/NGS2/group.activity/2nd_Group_activity/ref-HG38/GCA_000001405.28_GRCh38.p13_genomic.adj.SDF

###########################################################################
#we reached the previuos step and we didn't have enough time to proceed.
``` 



```s
## Code to check for the matching header names between bwa index, and NCBI reference  
conda create -y --name ngs1 python=3.6
conda activate ngs1
conda install -c bioconda -y bwa

workdir=$(pwd)
mkdir -p $workdir/fqData && cd $workdir/fqData
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/061/SRR14182261/SRR14182261_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/061/SRR14182261/SRR14182261_2.fastq.gz"

wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz
tar -xvf *.tar.gz

R1="SRR14182261_1.fastq.gz"
R2="SRR14182261_2.fastq.gz"

conda install -c bioconda -y bwa
/usr/bin/time -v bwa mem -t 8 -M GCA_000001405.15_GRCh38_no_alt_analysis_set.fna $R1 $R2 > SRR14182261.sam

conda install -y samtools
samtools view -hbo SRR14182261.bam SRR14182261.sam
samtools sort SRR14182261.bam -o SRR14182261.sorted.bam

conda install -c bioconda java-jdk=8.0.112
conda install -c bioconda picard
picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0
conda install -c bioconda gatk4


wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
#gunzip -c GCA_000001405.28_GRCh38.p13_genomic.fna.gz |awk -F '[ ,]' '{if ($1 ~ ">" ) print ">chr"$5; else print $0 }' > GCA_000001405.28_GRCh38.p13_genomic.adj.fna
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt
conda install -c conda-forge dos2unix
dos2unix GCA_000001405.28_GRCh38.p13_assembly_report.txt


wget 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_all.vcf.gz'
wget 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_all.vcf.gz.tbi'


grep "^0 chr" GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.ann | awk '{print $2}' > index.header
grep -v "^@" SRR14182261.sam | awk '{print $3}' | sort | uniq > sam.headers
gunzip -c GCA_000001405.28_GRCh38.p13_genomic.fna.gz | grep "^>" > fna.headers ## NCBI IDs
grep "^>" GCA_000001405.28_GRCh38.p13_genomic.adj.fna > adj.fna.headers  ## xxxx
gunzip -c HG002_GRCh38_1_22_v4.2.1_all.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq > vcf.headers   ## chr1-22

grep -v "^#" GCA_000001405.28_GRCh38.p13_assembly_report.txt | awk -F"\t" '{print $5,$10}' > GenBank_UCSC.map
cat index.header | grep -Fwf - GenBank_UCSC.map > GenBank_UCSC.index_map

conda install -c bioconda seqkit
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz
cat GenBank_UCSC.index_map | awk '{print $1}' > GenBank.targets
cat GCA_000001405.28_GRCh38.p13_genomic.fna | awk '{print $1}' > temp.fna
seqkit grep -n -f GenBank.targets temp.fna > GCA_000001405.28_GRCh38.p13_genomic.sel.fna && rm temp.fna
awk 'FNR==NR{a[">"$1]=$2;next}{if(a[$1])print ">"a[$1];else print $1}' GenBank_UCSC.index_map GCA_000001405.28_GRCh38.p13_genomic.sel.fna > GCA_000001405.28_GRCh38.p13_genomic.sel.newIds.fna


```




### Salma Solution
```s
Each group should take the time to read and discuss this paper "Comparing Variant Call Files by vcfeval in RTG Tools" (https://www.biorxiv.org/content/10.1101/023754v2.full.pdf). Use GATK and BCF-tools for single sample variant calling for each of the 5 alignments generated in the previous assembly. Use vcfeval to compare the outputs

We have 5 sorted bam files from the previous activity which have five different depths corrosponding to 2x,5x,10x,20x, and ~48x denoted as fractions of subsetted reads; 0.0408, 0.1022, 0.204, 0.4 and withoutfractioning(whole) which was 48x. 

SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.bam
SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.bam
SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.bam
SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.bam
SRR16939752_WES_of_Jianwei_Dong.sorted.bam

#note all commands were run for each individual file and no loops were used as each command -for each file- takes too long, so as to decrease the chances of 
incomplete/failed runs the commands were done for each file individually
##1. first we need to add read group information to the bams as this wasn included
java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups \
       I=SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.dedup.bam \
       O=SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.RG.dedup.bam \
       RGID=1 \
       RGSM=SRR16939752_WES_of_Jianwei_Dong_0.0408 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1

java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups \
       I=SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.dedup.bam \
       O=SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.RG.dedup.bam \
       RGID=1 \
       RGSM=SRR16939752_WES_of_Jianwei_Dong_0.1022 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1

java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups \
       I=SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.bam \
       O=SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.RG.bam \
       RGID=1 \
       RGSM=SRR16939752_WES_of_Jianwei_Dong_0.204 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1

java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups \
       I=SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.dedup.bam \
       O=SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.RG.dedup.bam \
       RGID=1 \
       RGSM=SRR16939752_WES_of_Jianwei_Dong_0.4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1

java -Xmx2g -jar $picard_path/picard.jar AddOrReplaceReadGroups \
       I=SRR16939752_WES_of_Jianwei_Dong.sorted.bam \
       O=SRR16939752_WES_of_Jianwei_Dong.sorted.RG.bam \
       RGID=1 \
       RGSM=SRR16939752_WES_of_Jianwei_Dong\
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1

2. next we need to remove duplicate reads since this is exome sequencing so we expect while preparing the library lots of duplicate reads were created
picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0

java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates \ 
INPUT=SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.RG.bam \ 
OUTPUT=SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.RG.dedup.bam \ 
METRICS_FILE=SRR16939752_WES_of_Jianwei_Dong_0.0408.metrics.txt

java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates \ 
INPUT=SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.RG.bam \ 
OUTPUT=SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.RG.dedup.bam \ 
METRICS_FILE=SRR16939752_WES_of_Jianwei_Dong_0.1022.metrics.txt

java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates \ 
INPUT=SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.RG.bam \ 
OUTPUT=SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.RG.dedup.bam \ 
METRICS_FILE=SRR16939752_WES_of_Jianwei_Dong_0.204.metrics.txt

java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates \
INPUT=SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.bam \ 
OUTPUT=SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.RG.dedup.bam \ 
METRICS_FILE=SRR16939752_WES_of_Jianwei_Dong_0.4.metrics.txt

java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates \
INPUT=SRR16939752_WES_of_Jianwei_Dong.sorted.RG.bam \ 
OUTPUT=SRR16939752_WES_of_Jianwei_Dong.sorted.RG.dedup.bam \ 
METRICS_FILE=SRR16939752_WES_of_Jianwei_Dong.metrics.txt

3. now index bam files and ref for GATK 
java -Xmx60g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.RG.dedup.bam
java -Xmx60g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.RG.dedup.bam
java -Xmx100g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.RG.dedup.bam
java -Xmx150g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.RG.dedup.bam
java -Xmx150g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR16939752_WES_of_Jianwei_Dong.sorted.RG.dedup.bam

index ref:
java -Xmx150g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.cds.all.fa O=Homo_sapiens.GRCh38.cds.all.dict
samtools faidx Homo_sapiens.GRCh38.cds.all.fa

4. Now do GATK variant calling using haplotype caller

gatk --java-options "-Xmx150G" HaplotypeCaller \
  -R ../Homo_sapiens.GRCh38.cds.all.fa -I SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.RG.dedup.bam  \
  --sites-only-vcf-output TRUE \
  -O GATK_0.0408.vcf

gatk --java-options "-Xmx150G" HaplotypeCaller \
  -R ../Homo_sapiens.GRCh38.cds.all.fa -I SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.RG.dedup.bam  \
  --sites-only-vcf-output TRUE \
  -O GATK_0.1022.vcf

gatk --java-options "-Xmx150G" HaplotypeCaller \
  -R ../Homo_sapiens.GRCh38.cds.all.fa -I SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.RG.dedup.bam  \
  --sites-only-vcf-output TRUE \
  -O GATK_0.204.vcf

gatk --java-options "-Xmx150G" HaplotypeCaller \
  -R ../Homo_sapiens.GRCh38.cds.all.fa -I SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.RG.dedup.bam  \
  --sites-only-vcf-output TRUE \
  -O GATK_0.4.vcf

gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R ../Homo_sapiens.GRCh38.cds.all.fa -I SRR16939752_WES_of_Jianwei_Dong.sorted.RG.dedup.bam  \
  --sites-only-vcf-output TRUE \
  -O GATK.vcf

5.Variant calling using bcftools mpileup and call
for file in *.sorted.RG.dedup.bam;do
  samtools index $file
done

index genome
bwa index Homo_sapiens.GRCh38.cds.all.fa


bcftools mpileup -Ou -f ../Homo_sapiens.GRCh38.cds.all.fa \
SRR16939752_WES_of_Jianwei_Dong_0.0408.sorted.RG.dedup.bam |\
bcftools call -Ov -mv > bcftools_0.0408.vcf

bcftools mpileup -Ou -f ../Homo_sapiens.GRCh38.cds.all.fa \
SRR16939752_WES_of_Jianwei_Dong_0.1022.sorted.RG.dedup.bam |\
bcftools call -Ov -mv > bcftools_0.1022.vcf

bcftools mpileup -Ou -f ../Homo_sapiens.GRCh38.cds.all.fa \
SRR16939752_WES_of_Jianwei_Dong_0.204.sorted.RG.dedup.bam |\
bcftools call -Ov -mv > bcftools_0.204.vcf

bcftools mpileup -Ou -f ../Homo_sapiens.GRCh38.cds.all.fa \
SRR16939752_WES_of_Jianwei_Dong_0.4.sorted.RG.dedup.bam |\
bcftools call -Ov -mv > bcftools_0.4.vcf

bcftools mpileup -Ou -f ../Homo_sapiens.GRCh38.cds.all.fa \
SRR16939752_WES_of_Jianwei_Dong.sorted.RG.dedup.bam |\
bcftools call -Ov -mv > bcftools.vcf

#If done in form of loop
#for file in *.sorted.RG.dedup.bam;do
#  sample=${file%.sorted.RG.dedup.bam} \
#  bcftools mpileup -Ou -f ../Homo_sapiens.GRCh38.cds.all.fa \
#  $file | bcftools call -Ov -mv > $sample.bcftools.vcf
#done

6. now compare called variant sets using vcfeval

7. first convert ref format to sdf 
rtg format -o Homo_sapiens.GRCh38.cds.all_SDF Homo_sapiens.GRCh38.cds.all.fa

8. Prepare vcfs for vcfeval
rtg bgzip *.vcf
rtg index -f vcf *.vcf.gz

9.compare variant sets (using 2 different variant callers) using vcfeval
#vcfs called from bams of highest depths
rtg vcfeval -b GATK.vcf.gz -c bcftools.vcf.gz -t Homo_sapiens.GRCh38.cds.all_SDF -o whole
#using vcfs called from bams of lowest depths
rtg vcfeval -b GATK_0.0408.vcf.gz -c bcftools_0.0408.vcf.gz -t Homo_sapiens.GRCh38.cds.all_SDF -o d_0.0408
#comparing the called sets from vcfs of highest and lowest depths using each of the 2 aligners
rtg vcfeval -b GATK.vcf.gz -c GATK_0.0408.vcf.gz -t Homo_sapiens.GRCh38.cds.all_SDF -o bigvssmall_GATK
rtg vcfeval -b bcftools.vcf.gz -c bcftools_0.0408.vcf.gz -t Homo_sapiens.GRCh38.cds.all_SDF -o bigvssmall_bcftools

#some statistics on varaints called from each caller and from different depths
(base) salma@metagenomics-bio:~/workdir/bams/RG/vcf/gunzip$ grep -v "^#" bcftools_0.0408.vcf | wc -l
21431
(base) salma@metagenomics-bio:~/workdir/bams/RG/vcf/gunzip$ grep -v "^#" GATK_0.0408.vcf | wc -l
28895

#using a loop
(base) salma@metagenomics-bio:~/workdir/bams/RG/vcf/unzip$ for i in *.vcf;do echo $i;grep -v "^#" $i | wc -l; done
bcftools_0.0408.vcf
21431
bcftools_0.1022.vcf
28008
bcftools_0.204.vcf
31106
bcftools_0.4.vcf
34288
bcftools.vcf
37162
GATK_0.0408.vcf
28895
GATK_0.1022.vcf
45536
GATK_0.204.vcf
54888
GATK_0.4.vcf
62426
GATK.vcf
70105
#we can see that with increasing depth, more variant records were observed and we also note more variants were called by GATK (much more than bcftools); this could be due to different representation of variants so a comparison tool such as vcfeval is needed to see how many variants are true positives by unifying the way the variants get represented

##since we get "run failed" error with vcfeval when using these vcfs -but with other data such as the tutorial's -dogchr5 vcfs we generated from different bams which worked fine with vcfeval, we don't have results to show 

###we tested vcfeval on variants called by bcftools and gatk on the dog bams that we worked with in the tutorial and it worked fine as below:

rtg format -o human_REF_SDF human_REF.fasta

rtg format -o dog_chr5_SDF dog_chr5.fa
rtg vcfeval -b GATK.vcf.gz -c bcftools.vcf.gz -t dog_chr_SDF -o eval


Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None               5782           5791       1699       1415     0.7732       0.8034     0.7880


It makes sense, that ~5000 variants are shared between the variants in both vcfs and some are false positives (seen in the callset 'bcftools' but not gatkvcf) and some are false nagtives (missed in bcftoolsvcf but seen in GATKvcf)

To compare this with a similar tool such as bcftools norm then intersect -unifyng variant representation we did:
first we checked just by intersecting both vcfs as they are represented
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf$ bcftools isec bcftools.vcf.gz GATK2.vcf.gz -p dir
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf$ cd dir
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/dir$ ls
0000.vcf  0001.vcf  0002.vcf  0003.vcf  README.txt
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/dir$ wc -l 0000.vcf
2271 0000.vcf
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/dir$ wc -l 0001.vcf
1894 0001.vcf
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/dir$ wc -l 0003.vcf
5279 0003.vcf
##again ~5000 variants shared

then we tried normalizing variant representation 

to use bcftools norm: we did
bgzip file.vcf
tabix file.vcf.gz

bcftools norm -f dog_chr5.fa ./vcf/gunzip/GATK2.vcf -o ./vcf/gunzip/GATK2.norm.vcf
bcftools norm -f dog_chr5.fa ./vcf/gunzip/bcftools.vcf -o ./vcf/gunzip/bcftools.norm.vcf

#here we notice 
with GATK.vcf there was no realigning when bcftools norm was used but for the bcftools.vcf 1272 lines of variants were realigned (changed)
Lines   total/split/realigned/skipped:  7115/0/0/0
Lines   total/split/realigned/skipped:  7490/0/1272/0


(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/gunzip/norm$ bcftools isec bcftools.norm.vcf.gz GATK2.norm.vcf.gz -p dir
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/gunzip/norm$ cd dir
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/gunzip/norm/dir$ ls
0000.vcf  0001.vcf  0002.vcf  0003.vcf  README.txt
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/gunzip/norm/dir$ wc -l 0000.vcf
1439 0000.vcf
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/gunzip/norm/dir$ wc -l 0001.vcf
1062 0001.vcf
(hackbio) salma@metagenomics-bio:~/workdir/GATK_tutorial/vcf/gunzip/norm/dir$ wc -l 0002.vcf
6117 0002.vcf

##so now when we do bcftools isec the realigned vcf files have more in common (because earlier some same variants were just represented differently)


##not sure why vcfeval shows less common variants than when using bcftools norm then bcftool isec as supposedly vcfeval maximizes common variants 



```

```s
picard_path="/media/mahmoud/storage/NU/NGS2/GA3/Project2/picard-2.26.9"
bams="./process/bams"
sample_list="./process/sams/sample_list.txt"
ref="./ref"
sams="./process/sams"
vcf="./process/vcf"
cd $sams
ls *.sam >> temp_sample_list.txt
cat  temp_sample_list.txt | sed 's/.sam//' >> sample_list.txt
cd ..
cd ..
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=$ref/Homo_sapiens.GRCh38.cds.all.fa O=$ref/ref.dict
#read -p "********************************************************DICT FILE IS DONE********************************************************************" p
samtools faidx $ref/Homo_sapiens.GRCh38.cds.all.fa
#read -p "********************************************************SAM INDEXING IS COMPLETE********************************************************************" s
rtg RTG_MEM=60G format -o ./ref/reference.sdf ./ref/Homo_sapiens.GRCh38.cds.all.fa
x=$(cat $sample_list)
for sample in $x
do
echo "********************************************************************************"
echo $sample
echo "********************************************************************************"
#read -p "enter to proceed to sorting" b
samtools view -hbo $bams/$sample.bam $sams/$sample.sam
samtools sort $bams/$sample.bam -o $bams/$sample.sorted.bam
#read -p "enter to proceed to PICARD" b
java  -Xmx60g -jar $picard_path/picard.jar MarkDuplicates INPUT=$bams/$sample.sorted.bam OUTPUT=$bams/$sample.dedup.bam METRICS_FILE=$bams/$sample'_metrics'.txt
#read -p "***********************************STEP 1 DONE *******************************************" z
java -Xmx60g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$bams/$sample.dedup.bam
#read -p "***********************************STEP 2 DONE *******************************************" v
java -Xmx60g -jar $picard_path/picard.jar ValidateSamFile I=$bams/$sample.dedup.bam MODE=SUMMARY >> ./process/bams/$sample.txt
gatk --java-options "-Xmx60G" HaplotypeCaller -R $ref/Homo_sapiens.GRCh38.cds.all.fa -I $bams/$sample.dedup.bam --emit-ref-confidence GVCF --pcr-indel-model NONE -O $vcf/$sample.gvcf
gatk --java-options "-Xmx60G" GenotypeGVCFs -R $ref/Homo_sapiens.GRCh38.cds.all.fa -V $vcf/$sample.gvcf --max-alternate-alleles 3 -O $vcf/$sample.vcf
bgzip $vcf/$sample.vcf
rtg index $vcf/$sample.vcf.gz
bcftools mpileup -Ou -f $ref/Homo_sapiens.GRCh38.cds.all.fa $bams/$sample.dedup.bam | bcftools call -Ov -mv > $vcf/bcftools/$sample.vcf
bgzip $vcf/bcftools/$sample.vcf
rtg index $vcf/bcftools/$sample.vcf.gz
rtg vcfeval --baseline $vcf/$sample.vcf.gz --calls $vcf/bcftools/$sample.vcf.gz --output ./process/qc/$sample --template ./ref/reference.sdf/
done

```
```s
#There are  sorted bam files for processing located at  /home/ngs/Downloads/Data
#create directory for gatk results
mkdir gatk_result && cd gatk_result
#marking duplicates
for sample in /home/ngs/Downloads/Data/*.bam
do
name=$(basename $sample | cut -d "." -f1)
java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt
done
#create directory for bcftools result 
cd ../
mkdir bcftools_result && cd bcftools_result
ref=$"/home/ngs/Downloads/bwa_align/bwa_index/Homo_sapiens.GRCh38.dna.chromosome.22.fa"
for sample in /home/ngs/Downloads/Data/gatk_result/*.dedup.bam
do
name=$(basename $sample | cut -d "." -f1)
bcftools mpileup -Ou -f $ref $sample | bcftools call -Ov -mv > $name.vcf
done
#indexing vcf files for subsequent use by vcfeval 
for sample in *.vcf
do
bgzip $sample
tabix -p vcf $sample.gz
done
#download a set of known variant for BQSR
cd ../gatk_result
wget 'ftp://ftp.ensembl.org/pub/release-105/variation/vcf/homo_sapiens/homo_sapiens_chr22.vcf.gz' -O homo_sapiens_chr22.vcf
#indexing sample
for sample in *.dedup.bam
do
java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$sample
done
#indexing reference file
ln -s /home/ngs/Downloads/bwa_align/bwa_index/Homo_sapiens.GRCh38.dna.chromosome.22.fa
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.chromosome.22.fa O=Homo_sapiens.GRCh38.dna.chromosome.22.dict
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.22.fa
#indexing known variants file
gatk IndexFeatureFile -I homo_sapiens_chr22.vcf.gz
#apply BQSR for samples
for sample in *.dedup.bam 
do
name=$(basename $sample | cut "." -d -f1)
gatk --java-options "-Xmx2G" BaseRecalibrator -R Homo_sapiens.GRCh38.dna.chromosome.22.fa \
-I $sample --known-sites homo_sapiens_chr22.vcf.gz -O $name.report
gatk --java-options "-Xmx2G" ApplyBQSR -R Homo_sapiens.GRCh38.dna.chromosome.22.fa \
-I $sample -bqsr $name.report -O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done
#variant calling with GATK HaplotypeCaller
for sample in *.bqsr.bam
do
name=$(basename $sample | cut -d "." -f1)
gatk --java-options "-Xmx2G" HaplotypeCaller -R Homo_sapiens.GRCh38.dna.chromosome.22.fa \
-I $sample --pcr-indel-model NONE -O $name.vcf
done
#indexing vcf for subsequent use by vcfeval 
for file in *.vcf
do 
bgzip $file 
tabix -p vcf $file.gz
done
#create a subset of high confiedence variants to use as baseline for comparison
cd ../Data
grep "^#" HG006_GRCh38_1_22_v4.2.1_benchmark.vcf > homo_sapiens_chr22.vcf
grep -v "^#" HG006_GRCh38_1_22_v4.2.1_benchmark.vcf | grep "^chr22" | sed 's/^chr22/22/' >> homo_sapiens_chr22.vcf
bgzip homo_sapiens_chr22.vcf
tabix -p vcf homo_sapiens_chr22.vcf.gz
#comparing out put of both tools using homo_sapiens_chr22.vcf subset for benchmarking
rtg format -o Homo_sapiens.GRCh38.dna.chromosome.22 Homo_sapiens.GRCh38.dna.chromosome.22.fa
ref=$"/home/ngs/Downloads/Data/Homo_sapiens.GRCh38.dna.chromosome.22"
for sample in ~/Downloads/Data/bcftools_result/*.vcf.gz
do
name=$(basename $sample | cut -d "." -f1 | cut -d "_" f2)
rtg vcfeval -b homo_sapiens_chr22.vcf -c $sample -t $ref -f QUAL -o bcf_eval_$name
done
for sample in ~Downloads/Data/gatk_result/*.vcf.gz
do
name=$(basename $sample | cut -d "." -f1 | cut -d "_" f2)
rtg vcfeval -b homo_sapiens_chr22.vcf -c $sample -t $ref -f QUAL -o gatk_eval_$name
done
#the result is 5 directories for each sample  with the name gatkeval_5x and bcfeval_5x for example
#plotting result of gatk and bcftools for each subsample
for i in original 20x 10x 5x 2x
do
rtg rocplot *$i/weighted_roc.tsv.gz --png $i
done




