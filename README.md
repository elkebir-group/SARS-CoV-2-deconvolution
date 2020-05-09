# Characterization of SARS-CoV-2 viral diversity within and across hosts

In light of the current COVID-19 pandemic, there is an urgent need to accurately infer the evolutionary and
transmission history of the virus to inform real-time outbreak management, public health policies and mitigation strategies. Current phylogenetic and phylodynamic approaches typically use consensus sequences,
essentially assuming the presence of a single viral strain per host. This code takes as input the variant allele frequencies of the mutations in multiple samples collected from infected hosts and finds the genotypes and proportions of the strains present in the hosts.

![Overview of Strain Deconvolution Problem](deconvolution.png)

In the STRAIN DECONVOLUTION problem, we are given the variant allele frequency (VAF) matrix F, containing the VAF of every mutation in each sample, and a number k of strains to be inferred. Our goal is to infer the genotype matrix B and mixture matrix U such that F ≈ BU, thus elucidating strains that occur within and across COVID-19 hosts along with their sample-specific proportions.

## Contents

  1. [Compilation instructions](#compilation)
     * [Dependencies](#dep)
     * [Compilation](#comp)
  2. [Usage instcructions](#usage)
     * [I/O formats](#io)
     * [Simulation](#simulate)
     * [Deconvolution](#gradient)
     * [Exposure](#exposure)

<a name="compilation"></a>
## Compilation instructions

<a name="dep"></a>
### Dependencies

Deconvolution solver is written in C++11 and requires a modern C++ compiler
(GCC >= 4.8.1, or Clang). In addition it has the following dependencies

* [CMake](http://www.cmake.org/) (>=3.1)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)
* [Eigen](http://eigen.tuxfamily.org/) (>= 3.3.7)
* [LBFGS++](https://lbfgspp.statr.me/) (>= 1.8.17)
* [GUROBI](https://www.gurobi.com/) (>= 8.1.0)

<a name="comp"></a>
### Compilation

To compile execute the following commands from the root of the
repository

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

In case CMake fails to detect LEMON, run the following command with adjusted paths:

    $ cmake -DLIBLEMON_ROOT=~/lemon

The compilation results in the following files in the `build' directory

EXECUTABLE       | DESCRIPTION
-----------------|-------------
`simulate`       | simuate perfect mixture data with missing entries for input
`gradient`       | solve strain deconvolution problem using penalty thresholding approach
`exposure'       | for a given set of strains, find the strain proportions in the sample

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats

The STRAIN DECONVOLUTION input is text based. There are two input files, the reference allele count file and the alternate variant allele count file.
Each line in both files correpsonds to a position in the genome.
Both the files should be tab-seperated.
The first 10 columns of the file have information about the mutations and its presence in the SRA and GISAID sequences. Only the position information is used in the algorithm to identify the mutation and the rest of the columns can be left empty.
The first 10 columns in the input files are --- '\<position\> \<reference allele\> \<variant allele\> \<gene\> \<N/S (nonsynonymous/synonymous)\> \<Amino acid change\> \<number of subclonal SRA samples\> \<number of clonal SRA samples\> \<total number of SRA samples\> \<number of consensus sequences\>'.

For the rest of the columns, each column is a sample and the entry in the reference file is the number of reads in the sample that have the reference allele and the entry in the variant file is the number of reads in the sample that have the variant allele.

For the EXPOSURE problem, there is an additional input of the genotype matrix B. This is also a tab-seperated value file in which the 
first column are the positions of the mutations. The subsequent columns correspond to the strains and an entry is 0 if the strain does not have the mutation and 1 is the strain has the mutation.

<a name="simulate"></a>

###  Sankoff Labeling (`simulate`)

    Allowed options:
      -h [ --help ]                  produce help message
      -k [ --strains ] arg (=10)     number of strains
      --missing arg (=0)             missing mutation rate
      -m [ --samples ] arg (=50)     number of samples
      -K [ --expstrains ] arg (=3)   number of expected strains per sample
      -N [ --expmutations ] arg (=5) number of expected mutations per strain
      -n [ --mutations ] arg (=100)  number of mutations
      -d [ --depth ] arg (=1000)     number of reads
      -s [ --seed ] arg (=0)         random number generator seed
      -o [ --output ] arg (=out)     output prefix


An example execution:

    $ ./simulate -s 5 --missing 0.1 -k 25 -o simuated_instance

<a name="gradient"></a>
### Deconvolution (`gradient`)

    Allowed options:
      -h [ --help ]               produce help message
      -k [ --strains ] arg        number of strains
      -l [ --lambda ] arg (=1.05) rate of increase of lambda
      --lambdaInit arg (=1)       initial value of lambda
      -B [ --initB ] arg          genotype matrix for initialization
      -U [ --initU ] arg          mixture matrix for initialization
      -m [ --maxIter ] arg (=100) maximum number of iterations
      -e [ --eps ] arg (=0.01)    termination condition
      -f [ --filter ]             filtering flag
      -s [ --seed ] arg (=0)      random number generator seed
      -T [ --threads ] arg (=1)   number of threads
      --input arg                 input files (ref and alt read counts)
      -o [ --output ] arg (=out)  output prefix

An example execution:

    $ ./gradient -s 0 -k 25 -o gradient_output -f simulated_instance_ref.tsv simulated_instance_alt.tsv --lambda 1.005 --eps 0 --maxIter 200
    
<a name="exposure"></a>
### Expsoure (`exposure`)

    Allowed options:
      -h [ --help ]              produce help message
      -k [ --strains ] arg       maximum number of strains allowed
      -B [ --initB ] arg         genotype matrix for initialization
      -T [ --threads ] arg (=1)  number of threads
      --input arg                input files (ref and alt read counts)
      -o [ --output ] arg (=out) output prefix

An example execution:

    $ exposure -k 25 simulated_instance_ref.tsv simulated_instance_alt.tsv -B gradient_output_B -o exposure_output
