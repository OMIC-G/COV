# SARS-COV-2 genomic surveillance in Galicia

The main goal of genomic surveillance of SARS-CoV-2 is to report the geographical spread and temporal patterns of SARS-CoV-2 clades and to detect resistance mutations.

Surveillance of respiratory viruses has a multidisciplinary character, and must integrate clinical, epidemiological and virological data in an organized way to report to the ECDC by TESSY (The European Surveillance System).

As a part of the RELECOV (Red de Laboratorios de secuenciación de SARS-CoV-2) and for coordination with CNM-ISCIII (Centro Nacional de Microbiología, Instituto de Salud Carlos III), several activities are part of the routine of a sequencing laboratory.

Here we present an integrated SARS-CoV-2 analysis pipeline including :
1. Sequencing and bioinformatic analysis
2. Submission to GISAID of consensus sequences
3. Report of clades to autonomous community Health Authorities (SERGAS)

This is the result of a continuous collaborative effort of the OMIC-G network (Red de Laboratorios para la aplicación de Ómicas a la Microbiología Clínica en Galicia)

## Dependencies

The following software was used for this pipeline.

* [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
* [samtools](https://www.htslib.org/)
* [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html)
* [nextclade](https://github.com/nextstrain/nextclade)
* [multiqc](https://multiqc.info/)
* [covCLI](https://gisaid.org/)


## Sample files and folder structure

The starting point of this analysis are .fastq files stored in the folder *fastq*

The following folders are used through the whole analysis:

* *fastq*: with the original .fastq.gz files
* *gisaid*: for uploading the final .fasta and .csv files to GISAID
* *mid*: temporary files created in the process and .fa files of each sample
* *out*: final files created in the analysis: report, fasta, csv, ...
* *refs*: .fasta reference files and primers bed file.


## 1. Sequence alignment

We align the sequences in each fast.gz file with `bwa-mem2` and study the coverage of each one with `samtools`. 

```
# Aligning:
bwa-mem2 mem -t 8 covid-bwa fastq/${SAMPLE}_L001_R1_001.fastq.gz fastq/${SAMPLE}_L001_R2_001.fastq.gz > mid/${SAMPLE}.sam

# Converting from .sam to .bam:
samtools view -bS mid/${SAMPLE}.sam -o mid/${SAMPLE}.bam

# Sorting:
samtools sort mid/${SAMPLE}.bam -o mid/${SAMPLE}.sorted.bam

# Triming primers:
ivar trim -i mid/${SAMPLE}.sorted.bam -b refs/SARs-CoV-2_v5.3.2_400.primer.bed -p mid/${SAMPLE}.trimmed -e -m 32

# Re-sorting:
samtools sort mid/${SAMPLE}.trimmed.bam -o mid/${SAMPLE}.sorted.bam

# Cleaning:
if [ -f "mid/${SAMPLE}.sorted.bam" ]
then
	rm mid/${SAMPLE}.sam
	rm mid/${SAMPLE}.bam
	rm mid/${SAMPLE}.trimmed.bam
fi
```

As a result, in the _mid_ folder we get all the .bam files of the samples with this name structure:
`${SAMPLE}.sorted.bam`


## 2. Generation of .fasta consensus of each sample

We use `samtools` and `ivar` in order to get the .fasta file of the consensus of each sample, using the previously obtained .sorted.bam files
```
samtools mpileup -aa -A -d 0 -Q 0 mid/${SAMPLE}.sorted.bam | ivar consensus -t 0.01 -p mid/${SAMPLE} -i ${SAMPLE}
```

As a result, in the _mid_ folder we get the .fa files with the same name structure:
`${SAMPLE}.fa`

Besides, we create a single .fasta file for all samples. We'll need these files later to query nextclade. 
```
cat mid/*.fa > out/all_SAMPLES.fasta
```

## 3. Pangolin typing

Clades clasification with `pangolin` database. It's a good idea to run `pangolin update` before.
```
pangolin out/all_SAMPLES.fasta --outfile out/pango.csv
```

## 4. Sumarize quality and coverage data

We obtain coverage data at 10X and 100X for each sample and store the results in a .csv file called _coverage.csv_. Later, we will use this file as part of the final report.

```
samtools coverage -H mid/${SAMPLE}.sorted.bam
samtools mpileup mid/${SAMPLE}.sorted.bam
```


## 5. Geting information of each sample from nextclade

We use the `nextclade CLI` from Nextstrain to obtain data for each sample. For this, we need to provide the .fasta file of all the samples

```
nextclade -i out/all_SAMPLES.fasta -c out/nextclade.csv
```

The results, stored as .csv files in the _out_ folder, also will be used as a part of the final report.


## 6. Variant calling

We get data for minority variants, so we can make histogram plots in the final report based on this data.

```
samtools mpileup -A -d 0 --reference refs/reference.fasta \
    -Q 0 mid/${SAMPLE}.sorted.bam | ivar variants -p out/tsv/${SAMPLE} \
    -t 0.03 -r refs/reference.fasta -m 10 -g refs/reference.gff
```


## 6. Final report

The final report is made with an R script, summarizing all data obtained previously plus quality data obtained with `FastQC`, `Qualimap` and `multiqc` tools.

```
samtools stats ${SAMPLE} > mid/${SAMPLE}.stats

fastqc fastq/*.fastq.gz --outdir=mid/fastqc/

qualimap multi-bamqc -d mid/qualimap_list.txt \
    -gff refs/reference.gff -outdir out/qualimap/ -r mid/*.sorted.bam

multiqc --force -o out -n "multiqc.html" -i "Report" \
    -b "<a href='report.html'>Global Report</a>" mid
```


## Consensus sequences submission to GISAID

The objective of this part is to ease the process to send the metadata and .fasta files to *GISAID* database. It's done with an R script and with the functionality of the `covCLI` application from GISAID. 

In order to use the `covCLI` command line utility, a GISAID user and a client_id it's needed.


### Metadata and fasta files

You need to complete the `gisaid_cov_template.csv` with the metadata of the samples (date, patient age, gender, Lab ID, ...). In our case, we use an .ods file from the LIS and process it with R to fill the template.csv.

The .fasta file with the sequence of each sample is done by concatenating the individual .fa files in the _mid_ folder. The names of each sequence in the final .fasta file, should match the names in the template.csv. Again, we use R to do the work.


### Uploading to GISAID

The files generated in the previous step (.fasta and .csv) are uploaded to GISAID with the `covCLI` utility. 

```
covCLI upload --username XXXX --password YYYY --clientid ZZZZ \
    --metadata metadata.csv --fasta sequences.fasta \
    --dateformat YYYYMMDD --log result.log
```

In the _result.log_ generated you'll have the accession_id of each sample uploaded to GISAID.


## Additional mutations

In this step we aim to obtain the additional mutations (the ones not defined in each subtype) for each sample.

The first step is to upload the .fasta file to https://clades.nextstrain.org and process the samples choosing as reference _BA.2_. Then export the .csv file with all the information of the samples. We store this .csv file as _nextclade.csv_ in the _out_ folder.

The process to obtain additional mutations is made in R using the web service from Cov-Spectrum

```{r}
# Load the nextclade.csv obtained from clades.nextstraing.org:
clades <- read.csv("out/nextclade.csv", sep = ";")

# Get defining mutations for each lineage
getMutations <- function (lin){
  refe <- httr::GET(paste0("https://lapis.cov-spectrum.org/open/v1/sample/aa-mutations?pangoLineage=",lin,"&aaMutations=,&downloadAsFile=false&dataFormat=csv"))
  refe <- data.frame(httr::content(refe, as="parsed", show_col_types=FALSE))
  refe$mutation
}

lineages <- unique(clades$Nextclade_pango)
mutas <- data.frame(lineages)
rownames(mutas) <- mutas$lineages
mutas$mutations <- apply(mutas, 1, getMutations)

# Obtainig the additional mutations of each variant
tbl <- matrix(ncol = 3, nrow = 0)
for(i in 1:nrow(clades)){
  muestra <- clades[i,]$seqName
  lin <- clades[i,]$Nextclade_pango
  mut <- explode(clades[i,]$aaSubstitutions, sep=",")
  ref <- mutas[lin, "mutations"]
  adi <- mut[!(toupper(mut) %in% toupper(ref[[lin]]))]
  row <- c(muestra, lin, paste(adi, collapse=","))
  tbl <- rbind(tbl, row)
}
tbl <- data.frame(tbl)
colnames(tbl) <- c("virus_name", "NEXTCLADE_LINEAGE","ADDITIONAL_MUTATIONS")

# We keep all the S gene mutations and remove most of the mutations in other regions:
tbl$ADDITIONAL_MUTATIONS <- gsub(pattern = "([ENM]:[A-Z]\\d+[A-Z],)|(ORF\\d[a-b]?:[A-Z]\\d+[A-Z],)", "", tbl$ADDITIONAL_MUTATIONS)
```

## Final report to SAÚDE PÚBLICA DE GALICIA

In this final step, an .xlsx file is created with the data requested by Saúde Pública de Galicia. 
Again, we use R to write the .xlsx file including:
- accession_id numbers from the covCLI log
- patient data and clade from the .ods from the LIS
- GISAID sample names from the .csv sent to GISAID
- additional mutations


