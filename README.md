# Variant calling from NGS data of two accessions of Lablab purpureus

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/enezermjema/Mini-project-group-06/master)

## Introduction

Lablab purpureus is a bean (family Fabaceae) commonly known as “lablab” which is native to Africa and widely cultivated in East Africa. It is called “Njahi” or black beans in Kenya where it is an important part of the daily diet. Lablab is, however, still an orphan crop with limited genomics and genetics resources. Recently, three genomes assemblies were separately created for two lablab accessions by the African Orphan Crops Consortium (AOCC) and ILRI (Feed and Forage Development program and Bioinformatics Community of Practice). These include an Illumina-based and an Oxford Nanopore-based assemblies for accession 147D, and a nanopore-based assembly for accession 147X. These assemblies can be used as templates for identifying molecular markers to fast-track the genetic mapping of economically important traits in lablab. 

## Objective
* To create a variant calling pipeline for Lablab purpureus.
* To convert the pipeline into a portable and reproducible workflow using Snakemake language.

## Data

Raw sequence data called SeqData_147D illumina based and SeqDataX Oxford Nanopore-based were provided. Also two reference files that is an Illumina-based and an Oxford Nanopore-based assemblies for accession 147D were provided. The data can be downloaded from [here](https://hpc.ilri.cgiar.org/~jbaka/EANBiT-RT2020-project6/)

## Setting up conda environment

. **Miniconda3**: A bootstrap version of Anaconda for python packages, to check if conda is installed

`conda -v`

If the version is not displayed then; then install (a 64bit system)

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh `

`sh Miniconda3-latest-Linux-x86_64.sh`

Follow on-screen instructions until the installation is complete

**NOTE: When asked to add conda_init , enter YES**

Add conda to PATH

`source ~/.bashrc `

If the installation is successful, you should see a list of installed packages with

`conda list`

If the command cannot be found, add conda to PATH environment manually, open the bashrc file and add the export PATH command to the end of the file and save it.

`sudo nano ~/.bashrc`

`export PATH=~/miniconda3/bin:$PATH`

## Usage:

Create a virtual environment called *variant* 

`conda env create --name variant --file Config.yml`

Activate the environment 

`source activate variant`


## Pipeline

Snakemake workflow language was used so as to have the pipeline portable and reproducible

[pipeline link]<https://github.com/enezermjema/Mini-project-group-06/blob/master/Snakefile>

## Software documentation in wiki pages following this link

[wiki link] <https://github.com/enezermjema/Mini-project-group-06/wiki>


## Contributors
* Immaculate Nahereza
* Jane Njeri
* Winfred Gatua 
* Nsangi Olga Tendo
* Davis Kiberu 
* Nsubuga Moses 
* Eneza Yoel
