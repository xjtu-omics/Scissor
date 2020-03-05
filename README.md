## Scissor: a flexible complex genome rearrangement simulator for short and long reads

Scissor is a flexible tool to simulate self-defined genome rearrangements, including chromothripsis and chromoplexy. Meanwhile, Scissor produces sequence dotplots of altered region for visualization and validation.

### Install and dependency

A python3.6 environment (suggested to install by Anaconda) with the following dependencies:

- intervaltree, pyfaidx, matplotlib, numpy. These packages are used to simulate variation genome.
- wgsim, pbsim. These two are used to simulate short and long reads, respectively.

All the above mentioned packages can be installed through conda command.

### General usage

#### Scissor SIM

SIM allows users to simulate self-defined complex genome rearrangement based on a template reference genome and produce a variation genome contain these complex events.

##### Rearrangement file

SIM requires a tab separated file with 4 columns, without header or any other lines. The columns are:

- Column 1: A simple name of the event. We use ***type1***, ***type2*** and etc. to name each events.
- Column 2: A comma separated symbolic representation of the reference allele. And please add ***#*** at the end of your representation. For example, ***A,B,C,D,#***.
- Column 3: A comma separated symbolic representation of the alteration towards the reference allele representation. 
  - We use ***^*** to indicate a inversion of the original sequence. For example, we can create a deletion-inversion by using ***A,C^,D***.
  - If you use a new symbol that dose not appear in the reference representation, such as ***E***. Then, the program will produce ***E*** based on either cut-paste or copy-paste from other places of the genome (intra and inter).
- Column 4: An integer to indicate the number of such rearrangement type you want to create.

An example rearrangement file based on the above description can be:

```
type1	A,B,C,D,#	A,C^,D	3
type2	A,B,C,D,#	A,C,A,D^	4
type3	A,B,C,D,#	A,C,B,E,B^	1
```

##### Chromosome size file

```
samtools faidx template_reference.fa
cut -f1,2 template_reference.fa.fai > template.sizes.tsv
```

##### Run SIM

Given the rearrangement file and chromosome size file prepared, we can run SIM with the following parameters:

- A template reference genome (-g). 
- Dimension of the template reference genome (-s).
- Some regions to exclude during the simulation (-x).
- Rearrangement file (-t).
- Output directory (-o).

The program will randomly assign size to each symbol, ranging from 500bp to 10Kbp. This can be set through -l and -u options.

#### Scissor short/long

Short/long read sequencing of the variation genome based on wgsim and pbsim.

##### Config file

A tab separated configuration file with four column for short read sequencing.

- Column1: chromosome.
- Column2: 0.
- Column3: dimension of this chromosome.
- Column4: variant allele fraction of the rearrangement. Determine the percentage of reads from variation genome and reference genome.

An example file:

```
chr20	0	64444167	50
```

##### Run short/long

```
# Short read sequencing
Scissor short -g tempalte_reference.fa -v variation_genome.fa -o output_dir/ -f config_file.txt

# Long read sequencing
Scissor long -g tempalte_reference.fa -v variation_genome.fa -o output_dir/ -f config_file.txt
```

### Use cases

Before we go to certain cases, it has to be aware that CGRs is a cluster of SVs, belonging to one haplotype. Therefore, Scissor only manipulate sequence segments on same haplotype. But you can provide two haplotypes in FASTA format and use one of the haplotypes for CGR simulation.

For example in the ***example_type.txt***, we define two types of CGRs as follows to simulate:

```
type1	A,B,C,D,#	C^,B,A,C,E	2
type2	A,B,C,D,#	A,C^,A^,B,D	1
```

##### Prepare template genome

1. Provide a normal reference genome. Scissor will treat this as a haploid genome by default, and implants CGRs.
2. Provide diploid genome named with ***h1*** and ***h2***. Scissor will use one haplotype to implant CGRs. Haplotypes ***h1*** and ***h2*** can be hacked by simple SVs, and ensure the SV regions are added to the exclude region file. *Please refer to VISOR if you want to simulate haplotype-aware simple SVs.*

Create the chromosome size file ***chrom.sizes.tsv*** for genome to hack. 

```
samtools faidx template.fa
cut -f1,2 template.fa > chrom.sizes.tsv
```

##### Run sim

By default, the excluded regions contain gaps, centromere, telomere and etc. The default file provided by SVelter (https://github.com/mills-lab/svelter/tree/master/Support).

```
cd /path/to/work_dir/
Scissor sim -g ./input/template.fa -s ./input/chrom.sizes.tsv -t ./input/example_type.txt -x ./input/exclude.bed -o ./output/
```

The output directory will have:

- ***variation_genome.fa***: hacked genome with CGRs.
- ***.png***: Dotplot view of each CGRs.
- ***gr_info.bed:*** Detailed information of CGRs, including ***REF*** and ***ALT*** position of each segments, the size of each segment as well as the start and end of the corresponding event. 

##### Run sequencing

In principle, once you have obtained the variation genome in FASTA format. You can also use modules from other tools for sequencing, such as VISOR SHORtS and LASeR. Here we introduce the default sequencing module provided by Scissor.

Prepare the configuration file. 

```
cd ./output/
find . -name "*.fa" -exec samtools faidx {} \;
cut -f1,2 *.fai > var_chrom.sizes.tsv
# Find the maximum dimenstion for accurate reads count.
cat var_chrom.sizes.tsv | sort  | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' > maxdims.tsv
# Allele fraction is set to be 100 in the case.
awk 'OFS=FS="\t"''{print $1, "0", $2, "100"}' maxdims.tsv > scissor_config.txt
```

Start sequencing with required I/O.

```
cd ./output/ && mkdir long_reads && mkdir long_reads
# Sequence short reads
Scissor short -g ./input/template.fa -v ./output/variation_genome.fa -o ./output_dir/short_reads/ -f ./output/config_file.txt

# Sequence long reads
Scissor long -g ./input/template.fa -v ./output/variation_genome.fa -o ./output_dir/long_reads/ -f ./output/config_file.txt
```

### Contact

Jiadong Lin: jiadong324@gmail.com

Songbo Wang: songbowang125@163.com