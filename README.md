### Raw data adaptor trimmed and quality filtered 
```sh
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
```

### De-novo assembly using five assemblers
#### Megahit Assembly with min contigs of 1000 bp
```sh
megahit -1 kneaddata/${sample}/${sample}.R1_kneaddata_paired_1.fastq.gz \
-2 kneaddata/${sample}/${sample}.R1_kneaddata_paired_2.fastq.gz \
-o megahit/${sample} -t $NSLOTS --min-contig-len 1000
```
#### Trinity Assembly with min contigs of 1000 bp
```sh
Trinity --seqType fq \
--left kneaddata/${sample}/${sample}.R1_kneaddata_paired_1.fastq.gz \
--right kneaddata/${sample}/${sample}.R1_kneaddata_paired_2.fastq.gz \
--CPU $NSLOTS --max_memory 760G --output trinity/${sample}_trinity --min_contig_length 1000
```
#### Meta-Spades Assembly
```sh
spades.py --pe1-1 kneaddata/${sample}/${sample}.R1_kneaddata_paired_1.fastq.gz \
--pe1-2 kneaddata/${sample}/${sample}.R1_kneaddata_paired_2.fastq.gz \
--meta -o spades/${sample} --threads $NSLOTS --memory 990 --only-assembler
```
###### Filter Meta-Spades assembly for >1000 bp contigs
```sh
seqkit seq -m 1000 spades/${sample}/contigs.fasta > spades/${sample}/contigs_1000.fasta
```
#### Meta-Viral-Spades Assembly
```sh
spades.py --pe1-1 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_1.fastq.gz \
--pe1-2 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_2.fastq.gz \
--metaviral -o metavSpades/${sample} --threads $NSLOTS --memory 990 --only-assembler
```
###### Filter Meta-Viral-Spades assembly for >1000 bp contigs
```sh
seqkit seq -m 1000 metavSpades/${sample}/contigs.fasta > metavSpades/${sample}/contigs_1000.fasta
```
#### RNA-Viral-Spades Assembly
```sh
spades.py --pe1-1 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_1.fastq.gz \
--pe1-2 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_2.fastq.gz \
--rnaviral -o rnavSpades/${sample} --threads $NSLOTS --memory 990 --only-assembler
```
###### Filter RNA-Viral-Spades assembly for >1000 bp contigs
```sh
seqkit seq -m 1000 rnavSpades/${sample}/contigs.fasta > rnavSpades/${sample}/contigs_1000.fasta
```
#### Concatenate contigs from the three assemblers
```sh
cat trinity/${sample}_trinity.Trinity.fasta megahit/kneaddata/${sample}/final.contigs.fa spades/${sample}/contigs_1000.fasta > de_novo_merge/${sample}_merge.fasta
```
#### Cluster contigs from all assemblies and select one representative from each cluster using CD Hit
```sh
cd-hit-est -i de_novo_merge/${sample}_merge.fasta -o de_novo_merge/${sample}_CD_HIT_c95.fasta -c 0.95 -n 10 -d 0 -M 0 -T $NSLOTS
```
