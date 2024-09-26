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
done
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
#### Meta-Spades Assembly NB this is a long and extremely demanding assembly, consider assemnling seperately
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
#### DeepMicroClass analysis and extract Eukaryote and Prokaryote Viruses
```sh
DeepMicroClass predict -i de_novo_merge/${sample}_CD_HIT_c95.fasta -o DeepMicroClass/${sample}
DeepMicroClass extract --tsv DeepMicroClass/${sample}/${sample}_CD_HIT_c95.fasta_pred_one-hot_hybrid.tsv --fasta de_novo_merge/${sample}_CD_HIT_c95.fasta  --class EukaryoteVirus --output DeepMicroClass/${sample}/${sample}_EukaryoteVirus.fasta
DeepMicroClass extract --tsv DeepMicroClass/${sample}/${sample}_CD_HIT_c95.fasta_pred_one-hot_hybrid.tsv --fasta de_novo_merge/${sample}_CD_HIT_c95.fasta  --class ProkaryoteVirus --output DeepMicroClass/${sample}/${sample}_ProkaryoteVirus.fasta
```
#### GeNomad analysis with conservative classification setting
```sh
genomad end-to-end --cleanup --conservative --threads $NSLOTS de_novo_merge/${sample}_CD_HIT_c95.fasta genomad/${sample} /scratch/wrbu/databases/genomad_db
```
#### Drop GeNomad viral sequences already present in DeepMicroClass classifications
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
#### Extract high-quality or complete viral sequences
```sh
grep -e "High-quality" -e "Complete" checkv/${sample}_genomad/quality_summary.tsv | awk '{print $1}' > checkv/${sample}_genomad/HQ_viruses.ID
cat checkv/${sample}_genomad/viruses.fna | seqkit grep -f checkv/${sample}_genomad/HQ_viruses.ID > checkv/${sample}_genomad/HQ_viruses.fasta
grep -e "High-quality" -e "Complete" checkv/${sample}_DeepMicroClassEukV/quality_summary.tsv | awk '{print $1}' > checkv/${sample}_DeepMicroClassEukV/HQ_viruses.ID
cat checkv/${sample}_DeepMicroClassEukV/viruses.fna | seqkit grep -f checkv/${sample}_DeepMicroClassEukV/HQ_viruses.ID > checkv/${sample}_DeepMicroClassEukV/HQ_viruses.fasta
grep -e "High-quality" -e "Complete" checkv/${sample}_DeepMicroClassProkV/quality_summary.tsv | awk '{print $1}' > checkv/${sample}_DeepMicroClassProkV/HQ_viruses.ID
cat checkv/${sample}_DeepMicroClassProkV/viruses.fna | seqkit grep -f checkv/${sample}_DeepMicroClassProkV/HQ_viruses.ID > checkv/${sample}_DeepMicroClassProkV/HQ_viruses.fasta
```
### Taxonomically classify high-quality/complete viral sequences
```sh
for classifer in DeepMicroClassEukV DeepMicroClassProkV genomad; do
# Diamond read classification
diamond blastx --db /scratch/wrbu/databases/diamond/nr --out checkv/${sample}_${classifer}/HQ_viruses --outfmt 100 \
-q checkv/${sample}_${classifer}/HQ_viruses.fasta \
--threads $NSLOTS -b20 --evalue 1e-6 -F 15 --range-culling --top 10
# Re-formatting for Megan software
daa-meganizer --in checkv/${sample}_${classifer}/HQ_viruses.daa --classify --mapDB /scratch/wrbu/databases/megan/megan-map-Feb2022.db --threads $NSLOTS --minSupport 1 --minPercentIdentity 40 --maxExpected 1.0E-6 --lcaAlgorithm longReads --lcaCoveragePercent 51 --longReads --readAssignmentMode readCount --propertiesFile /home/bourkeb/megan/MEGAN.vmoptions --only none 
done
```
### Extract reads based on classifications
```sh
  echo "+ Starting ${sample}"
for classifer in DeepMicroClassEukV DeepMicroClassProkV genomad; do
## Make a classification data folder
mkdir -p ${sample}_${classifer}/classify_virus/reads
  ###  
  # Takes the megan output (.rma) and filters for virus classifications (via CZID_taxonID.list file) and produces a Krona graph
  ###
  
  ## take the diamond alignment outputs and create a classification file
  daa2info -i ${sample}_${classifer}/HQ_viruses.daa -r2c Taxonomy > ${sample}_${classifer}/classify_virus/${sample}_classify.txt
  
  ## This first part cleans the read header names in the "pre-processed" fastq file (and the output is used to make a list of reads names with read lengths
  seqkit seq -i ${sample}_${classifer}/HQ_viruses.fasta | seqkit fx2tab --length --name --header-line > ${sample}_${classifer}/classify_virus/${sample}_length.txt
  
  ## This step adds a column read lengths to the classification file to create a "${sample}_table.txt" file
  awk -v OFS='\t' 'NR==FNR {a[$1] = $2; next} $1 in a {print $1, $2, a[$1]}' ${sample}_${classifer}/classify_virus/${sample}_length.txt ${sample}_${classifer}/classify_virus/${sample}_classify.txt > ${sample}_${classifer}/classify_virus/${sample}_table.txt
  
  ## This step filters the "${sample}_table.txt" file based on a virus TaxonID list creating a virus classification table file. These last three lines of code could be merged into one.
  awk 'NR == FNR {a[$1]; next} $2 in a {print}' Virus_TaxonID.list ${sample}_${classifer}/classify_virus/${sample}_table.txt > ${sample}_${classifer}/classify_virus/${sample}_virus.txt

  ###
  # Creates a table of virus species names and associated fasta read names
  ###
  
  awk '{print $2}' ${sample}_${classifer}/classify_virus/${sample}_virus.txt \
    | taxonkit lineage -r -L --data-dir /home/bourkeb/TaxonKit\
    | taxonkit reformat -I 1 -F -S -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" --data-dir /home/bourkeb/TaxonKit\
    | cut -f 9 \
    | csvtk add-header -t -n "species" \
    | csvtk pretty -t | tail -n+3 | paste ${sample}_${classifer}/classify_virus/${sample}_virus.txt - | cut -f1,4 > ${sample}_${classifer}/classify_virus/${sample}_virus_reads.txt 
   sed -i -e 's/ /_/g' -e 's/\//_/g' ${sample}_${classifer}/classify_virus/${sample}_virus_reads.txt  


  ###
  # Splits the previous file by virus species name
  ###

  for x in $(awk {'print $2'} ${sample}_${classifer}/classify_virus/${sample}_virus_reads.txt | sort -u); do
    egrep "[[:space:]]${x}$" ${sample}_${classifer}/classify_virus/${sample}_virus_reads.txt | awk {'print $1'} > ${sample}_${classifer}/classify_virus/reads/${sample}_${x}.txt
    echo "reads/\"${sample}_${x}.txt\" created"
  done


  ###
  # Uses the previous file to pre-processed fasta files to pull out the read names for each virus file
  ###
    for x in ${sample}_${classifer}/classify_virus/reads/${sample}_*.txt; do
      BASE=$(basename $x .txt)
      seqtk subseq ${sample}_${classifer}/HQ_viruses.fasta ${x} > ${sample}_${classifer}/classify_virus/reads/${BASE}.fasta
    done

  echo "- Ending ${sample}_${classifer}"
  echo
  echo
done
#
done
```
