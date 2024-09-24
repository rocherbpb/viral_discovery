BASE_DIR=/path/to/fastq/input/directory
HOST_GEN=/path/to/Host/genome/bowtie/ref/db
#
mkdir fastp kneaddata
#
for sample in $(cat sample_name.list); do
#
fastp -i ${BASE_DIR}/${sample}_R1_001.fastq.gz \
-I ${BASE_DIR}/${sample}_R2_001.fastq.gz \
-o fastp/${sample}.R1.fastq.gz \
-O fastp/${sample}.R2.fastq.gz \
--thread $NSLOTS --low_complexity_filter
#
kneaddata -i fastp/${sample}.R1.fastq.gz -i fastp/${sample}.R2.fastq.gz \
-db ${HOST_GEN} \
--output kneaddata/${sample} --bypass-trim -t $NSLOTS --run-fastqc-end
#
gzip kneaddata/${sample}/${sample}.R1_kneaddata_paired_1.fastq
gzip kneaddata/${sample}/${sample}.R1_kneaddata_paired_2.fastq
rm kneaddata/${sample}/*.fastq
#
done
