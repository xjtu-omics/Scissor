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

##### Outputs

***variation_genome.fa***: The variation genome hacked by the rearrangements.

***.png***: Dotplots image for each rearrangements.

#### Scissor SHORT/LONG

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

##### Run SHORT/LONG

```
# Short read sequencing
Scissor short -g tempalte_reference.fa -v variation_genome.fa -o output_dir/ -f config_file.txt

# Long read sequencing
Scissor short -g tempalte_reference.fa -v variation_genome.fa -o output_dir/ -f config_file.txt
```

### Use cases



### Contact

Jiadong Lin: jiadong324@gmail.com

Songbo Wang: songbowang125@163.com