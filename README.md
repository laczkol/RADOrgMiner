# RADOrgMiner

**RADseq Organellar DNA Miner and Genotyper**

If you find this pipeline useful, please do not forget to credit our work by citing this paper:

````
Laczkó, L., Jordán, S., Sramkó, G. (2021) A magnet to draw a bright needle out from the haystack – RADOrgMiner, an automated pipeline to genotype organellar reads from RADseq data, submitted
````

## Overview

The RADOrgMiner pipeline is designed to genotype sequence tags found in reduced genomic complexity libraries that originate from the organellar genome(s). It relies on aligning short sequences to a reference genome then, after defining loci by a minimal read depth (organellar genome is expected to have a higher read depth) calls the sequence of each locus. The pipeline can analyze multiple individuals simultaneously and the main outputs are fasta files that contain the alignments of loci. 

Running the analysis consists of two main steps: 1) mapping the reads to a reference genome and 2) calling haplotypes. After short read alignment, the pipeline emits some basic statistics about the read depth distribution along the reference sequence. This information helps to narrow down haplotype calling to those loci with low missingness over the dataset and are of extrachromosomal origin. Parameters of both steps can be fine-tuned using the command line. For the complete list of options see '**Details and example**' or the help menu of the pipeline that can be checked by starting the shell script RADOrgMiner.sh without specifying any additional parameters.

This pipeline is distributed without any warranty. Use it at your own risk and please, double-check if everything turned out well.

Any feedback is welcome. 

## Installation and dependencies

Downloading the RADOrgMiner.sh script or cloning this repository to the desired directory (by typing `git clone https://github.com/laczkol/RADOrgMiner.git`) in the command line is all that should be done to use the pipeline, if all dependencies are installed correctly and are added to the $PATH. The script was tested on Debian 10.1, but should run on any GNU/Linux-based distribution with dependencies installed. 

- makeblastdb and blastn from the [ncbi-blast+ package](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for self-blasting to identify inverted repeats

- [bedtools](https://bedtools.readthedocs.io/en/latest/) to mask the inverted repeat and create alignment intervals for haplotype calling

- [bwa](http://bio-bwa.sourceforge.net/) to align the short reads to a reference genome

- [samtools](http://www.htslib.org/) (please, make sure its version is at least 1.10) and [HTSlib](http://www.htslib.org/) to convert and subset the data files

- [freebayes](https://github.com/freebayes/freebayes) for haplotype calling

- [GNU Parallel](https://www.gnu.org/software/parallel/) to speed up the pipeline by splitting tasks and using threads in parallel

- vcf2fasta from the [vcflib](https://github.com/vcflib/vcflib) for file conversion

- [vcftools](http://vcftools.sourceforge.net/) to filter for missingness of loci 

- [muscle](https://www.drive5.com/muscle/) to align the loci

- [R](https://www.r-project.org/) and [Rscript](https://rdrr.io/r/utils/Rscript.html) for visualization

- [AMAS.py](https://github.com/marekborowiec/AMAS) to concatenate loci and check summary statistics

- [awk](https://www.gnu.org/software/gawk/manual/gawk.html) is an effective programming language that is used at many points of the pipeline to filter data and convert strings

- [GNU Core Utilities](https://www.gnu.org/software/coreutils/) is practically the spine of the pipeline. It is used in the majority of `bash` scripts and should be preinstalled on most of the GNU/Linux-based operating systems. `RADOrgMiner.sh` will not run at all without this.

  These dependencies can be installed with a package manager (e.g. `apt` or `conda`) or can be downloaded from GitHub and should be fairly easy to install. The presence of dependencies are checked before each run, but if some of them are missing the pipeline is not stopped. 

  If you use these tools within our pipeline, please do not forget to credit the authors of the above-mentioned softwares by citing them.

## Details and example

The help menu of RADOrgMiner.sh contains all the parameters necessary to set to run the pipeline. It can be invoked by just typing `RADOrgMiner.sh` in the terminal without specifying any parameters.

To run the pipeline, the two mandatory parameters to specify are the reference genome (`-r` or `--reference-genome`) and the population map (`-popmap` or `--population-map`):

````text
        This help menu is shown because you did not specify the reference genome and population map, the two mandatory parameters for this pipeline.
        -r --reference-genome
                                                        The file that contains the reference genome that should be used for both aligning the reads haplotype calling.
        -popmap --population-map
                                                        Tab-delimited text file with two columns: 1. identifier of samples without the file extension; 2. population identifier. The names in the first column must exactly match the file names that contain the reads. The population identifier will be used by freebayes. Please make sure IDs do not contain dots. This file can control which samples are being analysed.
````

The reference genome is used for aligning and haplotype calling, whereas the main goal to use a population map is to control which sample files should be included in the run. The population map should be a tab-delmited text file with sample_id<tab>population. File extensions should not be included in the first column, please specify only the sample names that should be analyzed. The second column can truly reflect the real populations, but, if desirable, any grouping can be specified (e.g. species, subspecies). Freebayes will use the grouping as a prior. Please make sure that empty lines and extra white spaces are not included in the population map file. Sample names also should not contain white spaces.

The working directories by default are assumed to be the current directory, but optionally can be specified by the following parameters:

````
        -s --samples-directory
                                                        The directory where fastq files are stored. The extension of the files can be .fq, .fastq,. fq.gz and fastq.gz. When paired-end reads are used, read pairs must be indicated by .1 and .2; e.g. SAMPLE_ID.1.fq.gz and SAMPLE_ID.2.fq.gz, which is the default of Stacks' process_radtags. [default is the current directory]

        -o --output-directory
                                                        Output directory to store all output files. Use absolute path. [default is the current directory]
````

Input files can have a `fq` or `fastq` extension and can be gzipped (`fq.gz` or `fastq.gz`).  The naming convention of paired-end reads follows the output of the `process_radtags` component of [Stacks](https://catchenlab.life.illinois.edu/stacks/) ([Rochette et al. 2019](https://doi.org/10.1111/mec.15253)). The pairs should be indicated by `.1` and `.2` before the file extension. If you specify an input or output directory, please use an absolute path (e.g. /home/user/dummy_reads) instead of a relative path (e.g. ../dummy_reads).

The main options for the first step of the pipeline, the alignment of reads are the following:

````
        -mask --mask-reference 
                                                        Valid options are yes or no (without quotes). Indicates if the inverted repeat of the chloroplast should be masked with N. Automatically removes only the largest duplication. [default no]

        -align --align-reads
                                                        Valid options are yes or no (without quotes). Indicates if reads stored in --samples-directory and named by first column of --population-map should be aligned to reference. [default yes]

        -np --number-of-processors
                                                        The number of processors (threads) that will be used. [default 1]

        -type --sequencing-type
                                                        Valid options are SE (single-end) or PE (paired-end). In case of PE read pairs should be separated to SAMPLE_ID.1.fq.gz and SAMPLE_ID.2.fq.gz. [default SE]
````

Masking could be important for plastomes, as inverted repeats (IR) (~30 kbp) are known to be present in many plastid genomes and might result in ambiguous alignments that will not be used in the haplotype calling step. Automatic masking removes the longest repeat after self-blasting and masks one of the IR regions with `N`-s. After masking, loci located in the IR are expected to have a roughly 2× coverage relative to the single-copy regions. The need for automatic masking should be explicitly specified by setting `-mask yes`. The default is to use the reference sequence as is.

Alignment of the reads is turned on by default. When only calling haplotypes of previously aligned samples, it might be useful to turn it off to reduce run time and avoid repeated tasks by setting `-align no`. The short read alignment can be fine-tuned by specifying the minimum seed length,  the matching-score, mismatch-penalty and gap-open penalty. These values can be set by the following parameters:  

````
        -bwa_k --min-seed-length
                                                        Minimum seed length for the Burrows-Wheeler aligner. [default 19]

        -bwa_A --matching-score
                                                        Matching-score for the Burrows-Wheeler aligner. [default 1]

        -bwa_B --mismatch-penalty
                                                        Mismatch-penalty for the Burrows-Wheeler aligner. [default 4]

        -bwa_O --gap-open-penalty
                                                        Gap-open penalty for the Burrows-Wheeler aligner. [default 6]
````



If not specified otherwise, the default values of `bwa` are used. Please, see the manual of `bwa` if you need further information about these parameters.

The run time can be decreased by using more CPU cores than the default 1. The number of threads to use **during the whole run** can be set by `-np` or `--number-of-processors`. 

The pipeline can use both forward and reverse reads. For single-end reads the default to set `-type SE` can be unchanged. To use both read-pairs, `-type PE` should be specified. Please, see naming convention of input files above.

The second main step of the pipeline is the calling of haplotypes that the following parameters can control:

````
        -call --call-haplotypes
                                                        Valid options are yes or no (without quotes). Indicates if haplotypes should be extracted to vcf and fasta format. The fasta format contains the entire sequences of the loci. It is suggested to first align reads to the reference, then investigate the read depth of loci. The bed files produced after alignment also report the bed coverage. 

        -mbc --min-bed-coverage 
                                                        Minimal coverage of a locus to be included in the output haplotypes. If the bed coverage of a locus is higher in any sample than this value, it will be a subject of analysis in all the samples. As organellar DNA can be overrepresented, it can be a really high number (e.g. 500). When trying different thresholds it is not needed to align the reads to the reference again, only the genotype call step should be redone. [default 3]  

        -sp --subsampling-prop
                                                        Float value between 0 and 1. To use less memory when calling haplotypes, each bed locus can be downsampled using samtools view by this proportion. 0.1 means: use 10% percent of all the reads found in a bed locus. The effect of downsampling is not tested properly, so use at your own risk and double check the results. [default 1.0]

        -mbs --min-bases-sequenced
                                                        Defines the minimum number of bases sequenced in a sample to be included in the genotype calling. [default 0]

        -mrc --min-read-count
                                                        Defines the minimum number of reads needed in a sample to be included in the genotype calling. [default 1]
````

The calling of haplotypes is turned off by default. To turn it on `-call yes` should be specified. The amount of loci is controlled by the minimal coverage of a locus (`-mbc` or `--min-bed-coverage`). If any locus of any individual included in the pipeline has higher coverage than this value, the locus will be included in the haplotype calling for all samples. Alignment intervals of each sample are merged at this step to get the total length of overlapping loci in the intervals used for the haplotype calling step. Theoretically, organellar reads can have a much higher read depth than the nuclear loci; thus, the value of `-mbc` may be really high. It is suggested to check the read depth distribution of samples after aligning the reads to the reference. The default value is 3, which probably will result in a scattered final dataset.

The number of samples can be narrowed down by specifying the minimum number of bases sequenced in any individual (`-mbs` or `--min-bases-sequenced`) and the minimal read count for any sample (`-mrc` or `--min-read-count`). Samples that fail these criteria will be omitted from the haplotype calling. The default for the minimum number of bases is 0, and for the minimum number of reads is 1, which parameters will only exclude samples with zero reads aligned to the reference.

The effect of setting different subsampling proportions (`-sp`) is currently not extensively tested. If subsampling is turned on, the given proportion of reads of each loci will be downsampled randomly for the haplotype calling. If you set this, please try more different values to check for consistency and double-check the results to ensure no real variability is lost. Also, pay attention to not to downsample the coverage lower than the minimal read depth that `freebayes` even considers analyzing. This parameter is here for testing reasons only.

The pipeline will carry out a subsetting of the dataset, but it does not involve the downsampling of read depth, just the exclusion of alignments outside the alignment intervals (loci) set by `-mbc`. By default, read depth of the desired loci will be left unchanged. This is done for an easier handling of loci in later steps.

The haplotype calling by freebayes can be (similarly to the alignment step) fine-tuned by setting the following parameters:

````
        -p --fb-ploidy
                                                        Ploidy level used by freebayes to call haplotypes. For organellar DNA it should be 1. [default 1]

        -q --fb-min-base-quality
                                                        Minimal base quality to include in the analysis. [default 20]

        -Q --fb-mismatch-base-quality-threshold
                                                        Minimal quality of a mismatched base to be accepted. [default 15]

        -m --fb-min-mapping-quality
                                                        Minimal mapping quality for a read to be included in haplotype calling. [default 30]

        -F --fb-min-alternate-fraction
                                                        Minimal supporting fraction of reads for an alternate allele. [default 0.4]

        -mc --fb-min-coverage
                                                        Minimal coverage to accept when calling haplotypes. [default 5]

        -n --fb-best-alleles
                                                        The number of most probable alleles to consider. [default 5]

````

As organellar loci are expected to be haploid, the `-p` (or `--fb-ploidy`) value should be left unchanged. This parameter can be set to a different value for use cases other than originally intended.

Other parameter values used by `freebayes` are set for a high accuracy and fairly short run time but can be altered via the command line as shown above. If you need further details on the parameters, please see the manual of `freebayes`.

A final filtering of the loci can be done by setting the following parameters:

````

        -minlen --min-locus-length
                                                        Minimal length of loci to include in the final dataset. [default 100]

        -miss --max-missing
                                                        Maximal missingness when filtering sites. [default 0.8]
````

The `-minlen` (or `--min-locus-length`) specifies the minimal length of any locus to be included in the final `muscle` alignment. This is by default set to 100 to exclude partial or ambiguous alignments from the dataset, but can be altered, especially if loci are expected to be relatively long (e.g. for ezRAD). 

The maximal missingness (`-miss` or `--max-missing`), similarly for SNP datasets, controls for the missingness of loci. Missingness is checked on the entire length of loci by `vcftools` and sites with a higher missingness than the specified will be excluded. Loci will be filtered for their minimal length after filtering for missingness.

To illustrate how the pipeline works this example uses dummy reads simulated from plastid genomes of *Arabidopsis thaliana* and *Arabidopsis lyrata*. These species are widely used model systems and numerous assembled plastomes are available in public databases. For this example run, 3-3 accessions of both species from the NCBI Nucleotide database were randomly picked. PE reads were simulated with [radinitio](http://catchenlab.life.illinois.edu/radinitio/) ([Rivera-Colón et al. 2021](https://doi.org/10.1111/1755-0998.13163)) using the sdRAD protocol with the PstI enzyme (cutsite: C^TGCAG) and requiring a 20× coverage with a mean insert size of 400bp. Some reads simulated using the same parameters that originate from the Chromosome 4 of *A. thaliana* (NC_003075) were added to the organellar read set to illustrate how the separation of the read sets work, but will not be analyzed separately in this example. The run will use the reference plastome of *A. thaliana* (NC_000932). The population map (saved as popmap) includes six samples and groups them by species:

````
lyr_KU764768_TGCA       a_lyrata
lyr_KX886355_TGCA       a_lyrata
thal_MK380719_TGCA      a_thaliana
thal_MK380720_TGCA      a_thaliana
thal_MK380721_TGCA      a_thaliana
lyr_NC_034379_TGCA      a_lyrata
````

The example assumes that all necessary files are stored in the current (working) directory:

````
popmap                      thal_MK380720_TGCA.1.fq.gz  lyr_KX886355_TGCA.2.fq.gz
NC_000932.fa                thal_MK380719_TGCA.2.fq.gz  lyr_KX886355_TGCA.1.fq.gz
thal_MK380721_TGCA.2.fq.gz  thal_MK380719_TGCA.1.fq.gz  lyr_KU764768_TGCA.2.fq.gz
thal_MK380721_TGCA.1.fq.gz  lyr_NC_034379_TGCA.2.fq.gz  lyr_KU764768_TGCA.1.fq.gz
thal_MK380720_TGCA.2.fq.gz  lyr_NC_034379_TGCA.1.fq.gz
````

To test the basic functionality of the pipeline, these files can be found in the directory `dummy_reads` on this repository. Sample names refer to the NCBI accession of the plastome that was used to generate the example reads.

 First, the example reads will be aligned to the reference by the following command:

````bash
RADOrgMiner.sh --mask-reference yes -r NC_000932.fa -align yes -call no -np 6 -popmap popmap -type PE
````

This command using six CPU cores (i7-4910MQ) should finish under six seconds with a maximum memory usage of 58 Mb. Parameters set for the run, the stage of analysis, and the alignment statistics created by `samtools coverage -m` are shown the `stdout`, which can be redirected to a text file to store it as a log. The alignment statistics for the first sample file in the example dataset should look like the following:

````
There are 8599 sites in the samples lyr_KU764768_TGCA with read depth > 0
NC_000932.1 (154.5Kbp)
>  54.24% │        █                               │ Number of reads: 680
>  48.22% │        █                               │ 
>  42.19% │        █                               │ Covered bases:   8.6Kbp
>  36.16% │        █                               │ Percent covered: 5.566%
>  30.13% │        █                ▃              │ Mean coverage:   0.646x
>  24.11% │    ▇   █    ▅           █     ▄        │ Mean baseQ:      40
>  18.08% │    █   █    █       ▄   █     █        │ Mean mapQ:       60
>  12.05% │▆   █   █    █       █   █     █        │ 
>   6.03% │█   █   █    █      ▃█   █     █        │ Histo bin width: 3.9Kbp
>   0.00% │█   █   █    █      ██   █     █        │ Histo max bin:   60.269%
          1       38.6K     77.2K     115.8K     154.5K 
Per site read depth can be found in working_directory/aligned/lyr_KU764768_TGCA.depth

````

As the masking of the reference was requested, the pipeline should have created the `blastn` database and the output of self-blasting (`NC_000932_selfblast.fmt6`  in this example) in the output directory. The range, that was masked prior the short read alignment, can be checked by opening the file ending with `IR_boundary` (`NC_000932_IR_boundary` in this example, it should be 128215-154478). The index file used for `bwa mem` are also stored here. The directory `unaligned` within the output directory contains all read pairs as `fq.gz` files that did not align to the reference genome. Using real data these can be processed further with any pipeline designed to process raw RADseq reads. In this example all `fq.gz` files stored in this directory should contain identical (n=500) reads. These are added to the example only to illustrate how the pipeline works. 

The alignments subject to genotype the loci of organellar origin can be found in the directory `aligned`. Alignments are stored here as `bam` files together with basic statistics about the alignment process. The read depth of each genomic position of each sample is exported to a file with the name of the sample with the extension `depth`. This information is summarized in the file `ind_depths.tsv`. At the end of the alignment step an R script is generated to visualize the read depth distribution of each sample, and the mean read depth of genomic positions across the dataset. The plots showing this information can be exported to `png` by running `Rscript /working_directory/aligned/plot_depth.R`, as suggested by the message output on the `stdout`. 

These plots should be checked to properly set the `-mbc` value. It should be noticed that not all loci overlap between the two species analyzed, resulting in potential missing data in the final dataset. This is expected for RADseq, as, if the cut site of the restriction enzyme is a subject of mutations, the experiment will result in allele dropout. Another table worth checking is the `coverage_table.tsv` that contains the statistics produced by `samtools coverage` for each individual.

The second step of the genotyping is the calling of haplotypes. After checking of the statistics of the alignments, an `-mbc` value of 15 seems to retrieve the sequence of all loci. The haplotype calling can be run with the following command:

````bash
RADOrgMiner.sh --mask-reference yes -r NC_000932.fa -align no -call yes -mbc 15 -np 6 -popmap popmap -type PE 
````

If you specify an output directory, please use the same for both the alignment and haplotype call steps. The script assumes that the directory `aligned`, that contains the files needed for haplotype calling, is located **within** the output directory. Using six CPU cores, haplotype calling of this dataset should take approximately 40 seconds with a maximal memory usage of 47 Mb. 

Some information worth checking can be found in the directory `aligned` . The `bam` files that contain `_subset` after the sample name only contain the alignments of the desired loci. The reference_subset.fa, which is subsetted from the original reference sequence by the file subset_alignments.bed, shows the reference's genomic regions that had the minimal specified number of reads aligned. Annotated haplotypes of each locus in vcf format can be found in the directory `aligned/vcf_loci`.

The main output of the pipeline can be found in the output directory, called `concat_loci.fa` that is the concatenated sequence of all loci filtered out from the alignment of organellar reads in standard fasta format. The range of each loci in this alignment can be retrieved from `concat_loci.parts`. The alignment of individuals' haplotype sequences can be found in the directory `fasta_loci`.  These files are in standard fasta format and are named to identify the location of loci easily. In this case, the analysis should yield 7 individual loci. It might be desirable to include only polymorphic in downstream analyses that. The variability of loci analyzed by `AMAS.py summary` could be checked by opening `summary.txt` located in `fasta_loci` with any spreadsheet editor software. Alternatively, the specified columns of this table can be checked by command-line utilities, e.g. running `cut -f 1,2,3,7,9 fasta_loci/summary.txt`, that should output the following:

````bash
 Alignment_name  No_of_taxa      Alignment_length        No_variable_sites       Parsimony_informative_sites
loc_NC_000932_1183-1700_0-517.fa        6       517     12      12
loc_NC_000932_121345-122411_0-1066.fa   6       1029    7       7
loc_NC_000932_16906-17868_0-962.fa      6       962     22      21
loc_NC_000932_32398-33473_0-1075.fa     6       1075    8       8
loc_NC_000932_33592-34562_0-970.fa      6       970     7       4
loc_NC_000932_52616-53705_0-1089.fa     6       1063    10      10
loc_NC_000932_97333-98493_0-1160.fa     6       1141    1       1
````

This analysis shows that all loci were polymorphic and also contain informative sites. If a subset of the loci is desirable, it can be easily achieved by using `AMAS.py` or any software that is designed for viewing and manipulating DNA sequences in fasta format.



