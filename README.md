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
#### Meta-Viral-Spades Assembly
```sh
spades.py --pe1-1 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_1.fastq.gz \
--pe1-2 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_2.fastq.gz \
--metaviral -o metavSpades/${sample} --threads $NSLOTS --memory 990 --only-assembler
```
#### RNA-Viral-Spades Assembly
```sh
spades.py --pe1-1 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_1.fastq.gz \
--pe1-2 kneaddata/${sample}_clean/${sample}.R1_kneaddata_paired_2.fastq.gz \
--rnaviral -o rnavSpades/${sample} --threads $NSLOTS --memory 990 --only-assembler
```
### Clean sequence headers and limit sequence length to >1000bp
#### Trinity
```sh
seqkit seq -i trinity/${sample}_trinity.Trinity.fasta | seqkit seq -m 1000 > trinity/${sample}_trinity_1000.Trinity.fasta
```
#### Megahit
```sh
seqkit seq -i megahit/kneaddata/${sample}/final.contigs.fa | seqkit replace -p '(.+)' -r 'megahit_${1}' | seqkit seq -m 1000 > megahit/kneaddata/${sample}/final_1000.contigs.fasta
```
#### Meta-Spades
```sh
seqkit replace -p '(NODE_\d+_length_\d+)_.+' -r 'Spades_${1}' spades/${sample}/contigs.fasta | seqkit seq -m 1000 > spades/${sample}/contigs_1000.fasta
```
#### Meta-Viral-Spades
```sh
seqkit replace -p '(NODE_\d+_length_\d+)_.+' -r 'metavSpades_${1}' metavSpades/${sample}/contigs.fasta | seqkit seq -m 1000 > metavSpades/${sample}/contigs_1000.fasta
```
#### RNA-Viral-Spades
```sh
seqkit replace -p '(NODE_\d+_length_\d+)_.+' -r 'rnavSpades_${1}' rnavSpades/${sample}/contigs.fasta | seqkit seq -m 1000 > rnavSpades/${sample}/contigs_1000.fasta
```
### Concatenate assemblies, cluster and select longest representative
#### Concatenate
```sh
cat \
trinity/${sample}_trinity_1000.Trinity.fasta \
megahit/kneaddata/${sample}/final_1000.contigs.fasta \
spades/${sample}/contigs_1000.fasta \
metavSpades/${sample}/contigs_1000.fasta \
rnavSpades/${sample}/contigs_1000.fasta \
> de_novo_merge/${sample}_merge.fasta
```
#### Cluster and sequence selection using CD Hit
```sh
cd-hit-est -i de_novo_merge/${sample}_merge.fasta -o de_novo_merge/${sample}_CD_HIT_c95.fasta -c 0.95 -n 10 -d 0 -M 0 -T $NSLOTS
```
### Viral classification of sequences
#### DeepMicroClass analysis and extract Eukaryotw and Prokaryote Viruses
```sh
DeepMicroClass predict -i de_novo_merge/${sample}_CD_HIT_c95.fasta -o DeepMicroClass/${sample}
DeepMicroClass extract --tsv DeepMicroClass/${sample}/${sample}_CD_HIT_c95.fasta_pred_one-hot_hybrid.tsv --fasta de_novo_merge/${sample}_CD_HIT_c95.fasta  --class EukaryoteVirus --output DeepMicroClass/${sample}/${sample}_EukaryoteVirus.fasta
DeepMicroClass extract --tsv DeepMicroClass/${sample}/${sample}_CD_HIT_c95.fasta_pred_one-hot_hybrid.tsv --fasta de_novo_merge/${sample}_CD_HIT_c95.fasta  --class ProkaryoteVirus --output DeepMicroClass/${sample}/${sample}_ProkaryoteVirus.fasta
```
#### GeNomad analysis with conservative classification setting
```sh
genomad end-to-end --cleanup --conservative --threads $NSLOTS de_novo_merge/${sample}_CD_HIT_c95.fasta genomad/${sample} /scratch/wrbu/databases/genomad_db
```
#### Drop GeNomad viral sequences already present in DeepMicroClass classification
```sh
cat DeepMicroClass/${sample}/${sample}_EukaryoteVirus.fasta DeepMicroClass/${sample}/${sample}_ProkaryoteVirus.fasta | seqkit seq -n -i > DeepMicroClass/${sample}/${sample}_DeepMicroClass_Virus.ID
cat genomad/${sample}/${sample}_CD_HIT_c95_summary/${sample}_CD_HIT_c95_virus.fna | seqkit grep -f DeepMicroClass/${sample}/${sample}_DeepMicroClass_Virus.ID --invert-match \
> genomad/${sample}/${sample}_CD_HIT_c95_summary/${sample}_DeepMicroClass_rem.fasta
```
### Check quality
```sh
checkv end_to_end genomad/conserve/${sample}/${sample}_CD_HIT_c95_summary/${sample}_DeepMicroClass_rem.fasta checkv/${sample}_genomad -d /scratch/wrbu/databases/checkv/checkv-db-v1.5 -t $NSLOTS
checkv end_to_end DeepMicroClass/${sample}/${sample}_EukaryoteVirus.fasta checkv/${sample}_DeepMicroClassEukV -d /scratch/wrbu/databases/checkv/checkv-db-v1.5 -t $NSLOTS
checkv end_to_end DeepMicroClass/${sample}/${sample}_ProkaryoteVirus.fasta checkv/${sample}_DeepMicroClassProkV -d /scratch/wrbu/databases/checkv/checkv-db-v1.5 -t $NSLOTS
```
#### Extract high-quality of complete viral genome sequences
```sh
grep -e "High-quality" -e "Complete" checkv/${sample}_genomad/quality_summary.tsv | awk '{print $1}' > checkv/${sample}_genomad/HQ_viruses.ID
cat checkv/${sample}_genomad/viruses.fna | seqkit grep -f checkv/${sample}_genomad/HQ_viruses.ID > checkv/${sample}_genomad/HQ_viruses.fasta
grep -e "High-quality" -e "Complete" checkv/${sample}_DeepMicroClassEukV/quality_summary.tsv | awk '{print $1}' > checkv/${sample}_DeepMicroClassEukV/HQ_viruses.ID
cat checkv/${sample}_DeepMicroClassEukV/viruses.fna | seqkit grep -f checkv/${sample}_DeepMicroClassEukV/HQ_viruses.ID > checkv/${sample}_DeepMicroClassEukV/HQ_viruses.fasta
grep -e "High-quality" -e "Complete" checkv/${sample}_DeepMicroClassProkV/quality_summary.tsv | awk '{print $1}' > checkv/${sample}_DeepMicroClassProkV/HQ_viruses.ID
cat checkv/${sample}_DeepMicroClassProkV/viruses.fna | seqkit grep -f checkv/${sample}_DeepMicroClassProkV/HQ_viruses.ID > checkv/${sample}_DeepMicroClassProkV/HQ_viruses.fasta
done
```

