# EBAME 2019:
# Metagenomics of a gueuze type lambic beer from Sussex.


## Introduction

Today we aim to investigate the community composition of a blended (gueuze) lambic beer brewed in the South Downs, UK. Traditional styles of beer and wine are commonly fermented using domesticated strains of yeasts (_Saccharomyces_ _cerevisiae_). While such strains show phenotypic diversity within and between different beers, wines, geographic regions and even breweries, most fermentations use a defined monoculture of _S. cervevisiae_ to give reproducible results. To maintain reproducibility, fermentation is undertaken in controlled anaerobic environments that limit exposure to contaminants from the environment which can lead to spoilage. 

Lambic beer is unlike traditional ales from the UK as the fermentation process is reliant on the natural seeding of the prepared fermentable (termed wort) by yeasts and bacteria from the local environment. Following the heating and boiling of grain to produce a range of fermentable sugars, the hot wort needs to be cooled before it is suitable for fermentation. There are many ways to achieve this under sterile conditions in modern breweries, however, one traditional method is to use "Coolships" which consisted of large open metal vats with a high surface to volume ratio. While highly effective in rapidly cooling the wort to fermentable temperatures, the open nature of coolships provides an opportunity for the natural microflora of the local area to enter the wort before it is sealed in barrels for fermentation. The beer we are examining today has been seeded and fermented with wild microflora originating from the Firle region of the south downs, UK.  


![alt text](https://github.com/BadgerRob/Staging/blob/master/AllagashCoolship_1200.jpg "Coolship")  
_Cooling wort in a coolship_  

![alt text](https://github.com/BadgerRob/Staging/blob/master/image3.jpeg "Firle Microflora")  
_The ingredients of the beer_  

![alt text](https://github.com/BadgerRob/Staging/blob/master/image4.jpeg "beer")
_The beer_  

## Setup  

Organise into groups of 3 - 4.  

## Data  

Sample fast5 files:
Sample quick fastq files:
Sample high accuracy fastq files:

## Basecalling

Nanopore sequencing results in fast5 files that contain raw signal data termed "squiggles". This signal needs to be processed into the `.fastq` format for onward analysis. This is undertaken through a process called 'basecalling'. The current program released for use by Oxford Nanopore is called `Guppy` and can be implemented in both GPU and CPU modes. Two forms of basecalling are available, 'fast' and 'high-accuracy' (HAC). HAC basecalling implements a 'flipflop' basecalling algorithm which is highly computationally intensive and thus slower than the fast basecalling method. Compare the two methods on the subset of fast5 files.  

```
guppy_basecaller -r --input_path path/to/fast5/ --save_path /path/to/fastq/ --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 4 --num_callers 2

```


```
guppy_basecaller -r --input_path path/to/fast5/ --save_path /path/to/fastq/ --config dna_r9.4.1_450bps_hac.cfg  --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 4 --num_callers 2

```

`guppy_basecaller`      : calls guppy  
`-r`                    : recursive  
`--input-path`          : path to fast5 dir/  
`--save-path`           : path to output fastq files  
`--qscore filtering`    : enable quality score filtering (optional)  
`--min_qscore`          : qscore filtering value (optional, default 7)  
`cpu_threads_per_caller`: number of threads to run.  
`--num_callers`         : number of basecallers to run.  
`--config`              : configuration file to run (HAC).  

Reads are output as `.fastq` files containing 4000 reads and quality data per file. Sequences contain a `@` followed by a header then sequence. This is separated from quality data by `+`. 

(Optional) If you want to you can watch in real time how many sequences are being written you can change to the directory where your fastq files are being written (/pass):

```

watch -n 10 'find . -name "*.fastq" -exec grep 'read' -c {} \; | paste -sd+ | bc'

```

`watch`                    : invoke 'watch' command : tells the computer to watch for something  
`- n 10`                   : every 10 seconds  
`find .`                   : look in every directory below where you are   
`-name "*.fastq"`          : target of find which will find every file with 'fastq' in its name  
`-exec`                    : execute the following command on find .  
`grep`                     : command to look for stuff within each file found by 'find . -name'  
`'read'`                   : the thing grep is looking for in the header of each sequence (1 per sequence)  
`-c` = count               : count the number of the occurrence 'read' in all files  
`{} \; | paste -sd+ | bc'` : paste the output from grep to basic calculator to display a count on the screen 


### Observations

How long did the different runs take?  
How do the identities differ with simple blast searching NCBI?  

## Read QC
Count the number of fastq reads in the Guppy pass dir.

```

cat pass/*.fastq | grep 'read=' - -c

```

Create a single file for onward analysis

```

cat path/to/pass/*.fastq > workshop.reads.fastq

```

Nanoplot


## Sorting and filtering reads 

FastqSample


## Taxonomic identification using Kraken2.

Kraken and Kraken2 provides a means to assign taxonomic identification to reads using a k-mer based indexing against a reference database. We provide a slimline reference database compiled for this workshop as well as the minikraken2 database. (ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz) Other databases such as the Loman labs [microbial-fat-free](https://lomanlab.github.io/mockcommunity/mc_databases.html) and [maxikraken](https://lomanlab.github.io/mockcommunity/mc_databases.html) are also available. 

### Optional extra information

Custom reference databases can be created using `kraken2-build --download-library`, `--download-taxonomy` and `--build` [commands](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases). Mick Wattson has written [Perl scripts](https://github.com/mw55309/Kraken_db_install_scripts) to aid in customisation. Example of customisation of databases can be found [here](http://porecamp.github.io/2017/metagenomics.html).




```

kraken2 --db path/to/kraken2_workshop_db/ --threads 8 --report path/to/output/report.txt path/to/workshop.reads.fastq > path/to/output

```

`kraken2`   :Call kraken2  
`--db`      :database name flag  
`--threads` :number of threads to use  
`--report`  : generate a user friendly report of taxonomy  

### Observations

Have a look at both the output and `report.txt` files using head and more to get a first look at the sample. Use `head` and `more` bash commands. Does this look correct for a lambic beer? Anything odd in the sample?

Try the assembly again using the minikraken2 database and see how your database can affect your results.

##Visualization of output

While scrolling through the kraken2 outputs can be fun and somewhat alarming, it does not produce a user friendly interpretation of the results. Here we present two methods to 'eye ball' your raw read diversity.

### Krona

Krona produces an interactive `.html` file based on your `--report` file. While not fully integrated with kraken2, the use of the report file gives an overall approximation of your sample diversity based on individual reads. Try this on the two kraken outputs.

```

ktImportTaxonomy -q 2 -t 3 report.txt -o kraken_krona_report.html

```

`ktImportTaxonomy`  : call KronaTools Import taxonomy  
`q 1 -t 3`          : For compatibility with kraken2 output  
`report.txt.`       : Kraken2 report.txt file  
`-o`                : HTML output  

Copy the html files to your local machine and open in a browser.

```

scp USERNAME@IP:/path/to/report.txt ~/Desktop

```

It should look something like this:  
![alt text](https://github.com/BadgerRob/Staging/blob/master/Krona.png "Krona report")


### Pavian

Pavian is an R based program that is useful to produce Sankey plots. It can be run on your local machine if you have R [installed](https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/installr.html). You may need to install `r-base-dev`. To set up Pavian open an R terminal and enter the following.

```
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

```

To run Pavian enter into R terminal:

```

pavian::runApp(port=5000)

```
You can now access Pavian at http://127.0.0.1:5000 in a web browser if it does not load automatically.  

Alternatively a shiny route is available.

```

shiny::runGitHub("fbreitwieser/pavian", subdir = "inst/shinyapp")

```

You should now be presented with a user interface to which you can browse for your report files on your local machine.

![alt text](https://github.com/BadgerRob/Staging/blob/master/pavin_snap.png "Pavian input")

Once you have loaded your file, navigate to the "sample" tab and try interacting with the plot. It should look something like this:

![alt text](https://github.com/BadgerRob/Staging/blob/master/kraken.png)

### Observations
How to the different databases affect your results?  
How does read depth affect your results?  
How does basecalling mode affect your results?  

##Assembly

```
minimap2 -x ava-ont -t 8 workshop.reads.fastq workshop.reads.fastq | gzip -1 > workshop.paf.gz

```

```
miniasm -f workshop.reads.fastq qorkshop.paf.gz > workshop.contigs.gfa
```

gfa to .fasta
```

awk '/^S/{print ">"$2"\n"$3}' workshop.contigs.gfa | fold > workshop.contigs.fasta

```

## Polishing with racon

```
minimap2 -t 8 -x map-ont workshop.contigs.fasta workshop.reads.fastq | gzip -1 > workshop.reads_to_assembly.paf.gz
```

Polish with racon

```

racon -t 12 workshop.reads.fastq workshop.reads_to_assembly.paf.gz workshop.contigs.fasta > workshop.contigs.racon.fasta

```
kraken2 on assembly
```
kraken2 --db path/to/kraken2_workshop_db/ --threads 8 --report workshop.contigs.racon.txt workshop.contigs.racon.fasta > workshop.contigs.kraken

```


Flye assembly

```

flye --nano-raw path/to/workshop.reads.fastq -g 1g -o flye_workshop/ -t 8 --meta

```



