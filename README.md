# Instability of the pseudoautosomal boundary in house mice

Supporting scripts for this manuscript:

> Morgan AP, _et al_. (2019) Instability of the pseudoautosomal boundary in house mice. bioRxiv [doi:10.1101/561951](https://doi.org/10.1101/561951).

Files are organized into one directory for each of the main arms of the analysis.

* **`./data`**: raw genotype data for intercross progeny; a kinder version is provided on [Figshare](https://gsajournals.figshare.com/s/28ec4963b5de2740e475)
* **`./genetic_map`**: inference of crossovers (and from that, a genetic map) in and near the PAR from genotype data in intercross progeny
* **`./chipseq`**: re-analysis of published ChIP-seq experiments in which we used sequence variants to assign read (pairs) to parental haplotypes

Code provided here is not "one-click" reproducible -- the ChIP-seq analyses in particular require BAM files from those experiments and a large [VCF file](ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/) from the Sanger Mouse Genomes Project. However, the material in this repository should be sufficient for interested readers to evaluate the analyses performed and judge their quality.
