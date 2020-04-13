<img src="https://github.com/jiadong324/Scissor/blob/master/Scissor.png" alt="scissor_logo" style="zoom:0.01%;" align=center/>

## Scissor: a flexible complex genome rearrangement simulator for short and long reads

Scissor is a flexible tool to simulate self-defined genome rearrangements, including chromothripsis and chromoplexy. Meanwhile, Scissor produces sequence dotplots of altered region for visualization and validation.

### Install and dependency

A python3.6 environment (suggested to install by Anaconda) with the following dependencies:

- Scissor sim dependencies: intervaltree, pyfaidx, matplotlib, numpy
- Scissor short/long dependencies: wgsim, pbsim, bwa, ngmlr, samtools

All the above mentioned packages can be installed through conda command. 

Afterwards, install Scissor from source.

```
git clone https://github.com/jiadong324/Scissor
cd Scissor
python setup.py install
```

### General usage

#### Scissor SIM

SIM allows users to simulate self-defined complex genome rearrangement based on a template reference genome and produce a variation genome contain these complex events.

##### Rearrangement file

SIM requires a tab separated file with 4 columns, without header or any other lines. The columns are:

- Column 1: A simple name of the event. We use ***type1***, ***type2*** and etc. to name each events.
- Column 2: A comma separated symbolic representation of the reference allele. Please add ***#*** at the end of your representation. For example, ***A,B,C,D,#***.
- Column 3: A comma separated symbolic representation of the alteration. **Note** that the first and last symbol are used as anchors for the breakpoint definition, so please manipulate or add symbols between them.
  - We use ***^*** to indicate a inversion of the original sequence. For example, we can create a deletion-inversion by using ***A,C^,D***.
  - If you use a new symbol that dose not appear in the reference representation, such as ***E***. Then, the program will produce ***E*** based on either cut-paste or copy-paste from other places of the genome (intra and inter).
- Column 4: An integer to indicate the number of such rearrangement type you want to create.

##### Run SIM

Given the rearrangement file and chromosome size file prepared, we can run SIM with the following parameters:

- A template genome to implant CGRs (-g). 
- Dimension of the template genome (-s).
- Some regions to exclude during the simulation (-x).
- Rearrangement file (-t).
- Output directory (-o).

The program will randomly assign size to each symbol, ranging from 500bp to 10Kbp. This can be set through -l and -u options. 

**Note**: Scissor only implants CGRs to chromosomes involved in the dimension file.

#### Scissor short/long

Short/long read sequencing of the variation genome based on wgsim and pbsim.

##### Config file

A tab separated configuration file with four column for short read sequencing. This file can be easily obtained from .fai file. 

- Column1: chromosome.
- Column2: start position.
- Column3: end position.
- Column4: variant allele fraction of the rearrangement. Determine the percentage of reads from variation genome and reference genome.

##### Run short/long

```
# Short read sequencing
Scissor short -g reference.fa -v alt_genome_folder/ -o output_folder/ -f scissor_config.txt -n sample_name

# Long read sequencing
Scissor long -g reference.fa -v alt_genome_folder/ -o output_folder/ -f scissor_config.txt -n sample_name
```

Scissor uses NGMLR for the alignment. You can also align sequenced reads under your output directory with other aligners.

### Use cases

Before we go to certain cases, it has to be aware that CGRs is a cluster of SVs, belonging to one haplotype. Therefore, Scissor only manipulate sequence segments on same haplotype. There are three FASTA files will be used in Scissor workflow:

- template.fa: sequence where CGRs are implanted.
- alt.h1.fa: sequence with CGRs.
- reference.fa: normal reference genome. It can be used as the template sequence.

#### Scissor common workflow

For example in the ***example_type.txt***, we define two types of CGRs as follows to simulate:

```
type1	A,B,C,D,#	A,C^,B,C,E,D	2
type2	A,B,C,D,#	A,C^,A^,B,D	1
```

##### Prepare template genome

1. Provide a normal reference genome. Scissor will use this as a haploid genome by default, and implants CGRs.
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
Scissor sim -g ./input/template.fa -s ./input/chrom.sizes.tsv -t ./input/example_type.txt -x ./input/exclude.bed -o ./path/to/work_dir/
```

The output directory will have:

- ***alt.h1.fa***: hacked genome with CGRs (h1 is default haplotype).
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

Start sequencing with required I/O. (The hacked FASTA file is under ./output/ directory)

```
cd ./output/ && mkdir short_reads && mkdir long_reads
# Sequence short reads
Scissor short -g ./input/reference.fa -v ./ -o ./output/short_reads/ -f ./output/scissor_config.txt

# Sequence long reads
Scissor long -g ./input/reference.fa -v ./ -o ./output_dir/long_reads/ -f ./output/scissor_config.txt
```

#### Scissor with VISOR

VISOR provides an efficient way to simulate haplotype-aware simple SVs and it is able to simulate different clones. Scissor is able to combine these two functions in its workflow, leading to full spectrum variation simulation.

##### Case 1: Full spectrum variants simulation

First, VISOR HACK module produces two haplotypes ***h1.fa*** and ***h2.fa*** with  ***hack.h1.bed*** and ***hack.h2.bed*** respectively (Please refer VISOR document for details).

Next, we use ***h1.fa*** and ***h2.fa*** as input to Scissor, and please specify the haploid with ***-i*** option to implant CGRs (h1 by default). With the required I/O, Scissor runs:

```
cd /path/to/work_dir/ && mkdir output
mv h1.fa ./output && mv h2.fa ./output/alt.h2.fa

cd output/ 
samtools faidx h1.fa
cut -f1,2 h1.fa.fai > chrom.sizes.tsv

# adding exising SV region to exclude file.
cat exclude.bed hack.h1.bed hack.h2.bed | sortBed > exclude_new.bed 

Scissor sim -g ./output/ -s ./output/chrom.sizes.tsv -t ./output/example_type.txt -x ./output/exclude_new.bed -o ./path/to/work_dir/output/
```

Next we sequence the diploid genome:

```
mkdir short_reads/ && mkdir long_reads

Scissor short -g ./path/to/reference.fa -v /path/to/work_dir/output/ -o ./output_dir/short_reads/ -f ./output/scissor_config.txt

Scissor long -g ./path/to/reference.fa -v /path/to/work_dir/output/ -o ./output_dir/long_reads/ -f ./output/scissor_config.txt
```

##### Case 2: Variants at different clones

First, get your variation genome implanted with CGRs, either haploid or diploid and put them in a single folder. Next, we can use VISOR modules with multiple FASTA folder as input to simulate variants at different clones.

### Contact

Jiadong Lin: jiadong324@gmail.com

Songbo Wang: songbowang125@163.com