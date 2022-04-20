# Allele-specific RNA-seq workflow

Nextflow workflow to process raw RNA-seq data for allele-specific expression analysis. It takes compressed fastq files as input and generates normalized signal tracks and counts of total ("Gall") and allelic reads. Also, it saves log files that can be parsed to calculate the number of reads filtered out in each step.

## Author

Yuvia Alhelí PÉREZ-RICO

Affiliation: European Molecular Biology Laboratory | EMBL

## Dependencies

To use this workflow, you will need to install the following programs, the indicated version is the one that I have used to analyze my data:

- nextflow (20.04.1)
- FastQC (v0.11.8)
- Trim Galore (0.6.3)
- STAR (2.7.2b)
- samtools (1.9)
- SNPsplit (0.3.4)
- Picard Tools (2.20.8)
- featureCounts (v2.0.1)
- R (3.2.2)
- deepTools (3.1.3)

Note: I suggest to install all programs in a conda environment.

## How to use the workflow

Thank you for your interest in using the workflow!

After downloading the main script, the configuration file and the additional script within the bin folder, you will need to prepare the following files and change the paths in the configuration file accordingly:

- An N-masked STAR index using SNPs of 2 strains.
- Gene annotations in GTF format.
- A compressed file with the SNPs of the 2 strains for SNPsplit.
- A comma-separated file indicating the sample name and paths to the fastq files (an example is available in the 'docs' folder). The name of the fastq files is relevant for line 235 of the main script.

If you are not using SLURM, then do additional modifications to the configuration file considering the workload manager that you use. Finally, write a simple bash script that will be submitted to the cluster to activate the conda environment and start the main nextflow job:

`source /home/user/miniconda2/bin/activate /home/user/conda-envs/ASE_RNA-seq`

`nextflow run allelic_RNA-seq.nf`

Please, keep in mind that this workflow was written for the processing of mouse paired-end data, therefore, you will need to do some changes to the main script before using it with data of other species, for example, changing the strain names. Also, there is a parameter called 'TXline' that I use to skip some chromosomes for the normalization of signal tracks, as I know that they tend to get duplicated in the cell lines that are commonly handled in our lab or to ease the visualization of differentiated and undifferentiated cells. Simply, type 'no' if you do not want to skip chromosomes.

## Acknowledgements

The workflow is a nextflow adaptation of a bash script written by Samuel Collombet, former postdoctoral fellow in the Heard group.

This repository is part of a project that has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 882771.

## License

Licenced under the GNU general public license.

