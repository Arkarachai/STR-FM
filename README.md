# *STR-FM*, a short tandem repeat profiling using a flank-based mapping approach

## User manual and guide
We designed the STR profiling pipeline as a collection of tools which can be executed in both commandline or via a GUI on Galaxy. The easiest way to use STR-FM pipeline is to via Galaxy platform. Current, we have all tools in Galaxy main toolshed (See installation of STR-FM tools from toolshed below) and in Galaxy test website (STR-FM: microsatellite analysis).

## Overview

Our tools in ‘str_fm’ can be used to: 

**(1) profile STRs from short read data with STR-FM pipeline** (tools: ‘STR detection’, ‘Read name modifier’, ‘Fetch bases flanking’, ‘Combine mapped faux paired-end reads’, ‘Check STR motif compatibility between reference and read STRs’, ‘Select uninterrupted STRs’)

This pipeline needs several tools on Galaxy to complete the process. It can be customized with different mapper or STRs detection algorithm. Either single-end or paired-end sequencing data can be utilized; for paired-end read data, each read is treated separately. The core of the pipeline consists of the following three procedures 

First, STR-FM runs a short-read STR detection tool using a string comparison algorithm (see publication details). The algorithm can detect exact (pure, or uninterrupted) STRs (mono- through hexanucleotide STRs greater than or equal to two repeats), incomplete motifs (e.g., ATATATA), interrupted STRs (e.g., AAAATAAAAA), or multiple STRs in a read. Reads that do not have sufficient upstream or downstream sequences flanking the STRs are discarded (we used a threshold of 20 bp on each side of an STR). Each read is split into two “pseudoreads,” containing the upstream and downstream flanks surrounding the STR. 

Second, these are mapped to the reference genome using a standard paired-end read-mapping algorithm, e.g., BWA, Bowtie, or Bowtie2, treating each pair of flanking sequences as a faux paired-end read. 

Finally, STR-FM runs a profiler tool, which groups all reads with STRs that are mapped to the same location in the reference genome. As a result, an array of all STR lengths from the reads mapping to a particular STR-containing locus is generated.

**(2) genotype STRs with error correction** (tool ‘Correct genotype for STR errors’)

This pipeline needs only one of our tools to complete process. It will take STR-profile file and sequencine error rates file as inputs. The program will calculate the maximum likelihood of genotype for each STR locus in STR-profile file. Then it will report the mostly likely genotype and the log odds ratio between their probabilities, which can be interpreted as a confidence of genotyping (the more this value deviates from 0, the more confidence we have in this genotype).

**(3) estimate the minimum informative read depth from error rates** (tools: ‘Generate all possible combination of STR length profile’, ‘Evaluate the probability of the allele combination to generate read profile’, ‘Combine read profile probabilities’)

This pipeline needs other tools on Galaxy to complete the process. This pipeline will generate all possible read profiles from sequencing error spectrum, select the profiles that can distinguish heterozygote from homozygote, calculate the probability to produce such profiles from sequencing error spectrum, and report the probability that a certain sequence depth can distinguish heterozygote from homozygote under a given sequencing error rates (see publication details). We recommend that you should try to run with less than 10x depth for initial trial.

**(4) convert informative read depth to locus-specific and genome-wide sequencing depth** (tool ‘Convert informative read depth to sequencing depth’).  

This pipeline needs only one of our tools to complete process. It will convert *informative read depth* to *locus-specific sequencing depth* (given read length) and *genome-wide sequencing depth* (given confidence intervals).


## Description of tools

The short description for each tool is provided below.

1. “STR detection” = Detect STRs from short reads (FASTQ), reference genome (FASTA), or alignments (SAM)
2. “Read name modifier” = Change space in read name to ‘_’ to prevent read name truncation by mapping tools
3. “Fetch bases flanking” = Generate two FASTQ files containing flanking bases around STRs for mapping as faux paired-end reads
4. “Combine mapped faux paired-end reads” = For each mapped faux paired-end reads, infer STR sequence in reference genome between the two mapped ends of the pair
5. “Check STR motif compatibility between reference and read STRs” = Check if two STRs have the same motif
6. “Select uninterrupted STRs” = Select STRs that do not contain an interruption
7. “Correct genotype for STR errors” = Build error correction model from pre-defined error rates and identify most likely genotype of the input data
8. “Generate all possible combination of STR length profile” = Use STR error spectrum to generate all possible combinations of read profile at each read depth
9. “Evaluate the probability of the allele combination to generate read profile” = Calculate the probability of a given genotype to generate read profiles (instead of finding most likely genotype like tool number 7)
10. “Combine read profile probabilities” = Sum the probability of the given allele combinations to generate read profile at certain read depth
11. “Convert informative read depth to sequencing depth” = Calculate ‘locus-specific’ and ‘genome-wide’ sequencing depth from the given informative read depth
The detailed description for each tool is embedded within the tool.

## Citing *STR-FM*
Fungtammasan A, Ananda G, Hile SE, Su MS, Sun C, Harris R, Medvedev P, Eckert K, Makova KD. 2015. Accurate Typing of Short Tandem Repeats from Genome-wide Sequencing Data and its Applications, Genome Research

## Installation of STR-FM tools from toolshed


The installation can be done as follows


1 Install and set configuration of local Galaxy 

1.1 Download and install Galaxy (https://wiki.galaxyproject.org/Admin/GetGalaxy). Galaxy works on both Unix and Mac OS.

1.2 From your Galaxy directory, add your E-mail as admin E-mail to the Galaxy configuration file. Depending on the Galaxy version, this file can be either universe_wsgi.ini or config/galaxy.ini (https://wiki.galaxyproject.org/Admin/Interface)

1.3 Set directory for tool dependencies (step 2 in https://wiki.galaxyproject.org/Admin/Tools/AddToolFromToolShedTutorial). 

1.4 Run local Galaxy from the command line by running ‘sh run.sh’ from your Galaxy directory. 

1.5 Open your Galaxy from your browser at address http://localhost:8080 (https://wiki.galaxyproject.org/Admin/GetGalaxy)

1.6 Register using your admin E-mail in the ‘User’ tab on the top.

1.7 Refresh your browser


2 Install tools and dependencies

2.1 From your local galaxy, click ‘Admin’ tab on the top.

2.2 On the left panel, click ‘Search and browse tool sheds’ under ‘Tool sheds’. ‘Accessible Galaxy tool sheds’ will appear on main panel.

2.3 Click on ‘Galaxy main tool shed’ and select ‘Browse valid repositories’. (https://wiki.galaxyproject.org/Admin/Tools/AddToolFromToolShedTutorial)

2.4 Type ‘str_fm in search box and click enter.

2.5 The ‘suite_str_fm_0_1’ repository that has ‘arkarachai-fungtammasan’ as the owner will appear. The user may click on this repository name and click ‘Preview and install’. The ‘Install to Galaxy’ button will appear on upper right corner. This button allows the user to install all our tools and workflows -- pipelines containing tools for specific purpose such as STR profiling from short read sequencing data, microsatellite detection of the reference genome, and estimating minimum informative read depth. None of our tools have any dependencies. However, some of the other tools that used in our workflows (e.g. SAM flag filter, unique element selection, etc.) are not included in the standard Galaxy installation. For the user’s convenience, we included all dependency tools for the workflows in this repository. Therefore, installing ‘suite_str_fm_0_1’ will be sufficient to operate all workflows we provided. 

2.6 After clicking on ‘Install to Galaxy’ and ‘Install’ button in confirmation page, all our tools, workflows, and test datasets will be downloaded to your local Galaxy. After the download is completed, all our tools will be available on your local Galaxy. If the user wants to use the workflows that we suggested (i.e. STR profiling from short read sequencing data, microsatellite detection of the reference genome, and estimating minimum informative read depth), please proceed to step 3.

2.7 Refresh your browser


3 Install workflows

3.1 Click on the ‘Admin’ tab at the top again.

3.2 On the right panel, click ‘Manage installed tool shed repositories’ under ‘Server’. ‘Installed tool shed repositories’ will appear on main panel.

3.3 Click to open ‘str_fm’ repository. 

3.4 Scroll down to ‘Workflows’ section and select the workflow that you want to install. The SGV graphic of the workflow will appear.

3.5 Click on the ‘Repository Actions’ on the upper right corner and select ‘Import workflow to Galaxy’. If success, the ‘Workflow <workflow name> imported successfully’ will appear. Once the workflow is imported to your Galaxy, you can view and modify it from ‘Workflow’ tab on the top. 
