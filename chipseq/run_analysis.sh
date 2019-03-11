# Assign ChIPseq reads from F1 male spermatocytes to maternal and paternal haplotypes,
#	then calculate lightly-smoothed coverage on each haplotype across the X chromosome.
# This uses H3K4me3 data from Baker et al. (2015) PLoS Genet, and SSDS data from Khil et al. (2016) Genes Dev

classify_reads.py \
	-v /proj/pmdvlab/data/variants/v5.snps.CC8.vcf.gz \
	-p CAST_EiJ C57BL_6NJ \
	-r chrX \
	-b aln/FB_1.chip.bam \
	-o allele_specific/H3K4me3/FB.u.bam && \
samtools sort allele_specific/H3K4me3/FB.u.bam >allele_specific/H3K4me3/FB.bam && \
samtools index allele_specific/H3K4me3/FB.bam

reads_by_tag.py -t XX -v m allele_specific/H3K4me3/FB.bam | \
samtools view -bhS >allele_specific/H3K4me3/FB.m.bam && \
samtools index allele_specific/H3K4me3/FB.m.bam

reads_by_tag.py -t XX -v p allele_specific/H3K4me3/FB.bam | \
samtools view -bhS >allele_specific/H3K4me3/FB.p.bam && \
samtools index allele_specific/H3K4me3/FB.p.bam

bamCoverage -b allele_specific/H3K4me3/FB.m.bam -of bigwig -o allele_specific/H3K4me3/FB.m.bw -r chrX:0:171031299 -bs 100 --smoothLength 500
bamCoverage -b allele_specific/H3K4me3/FB.p.bam -of bigwig -o allele_specific/H3K4me3/FB.p.bw -r chrX:0:171031299 -bs 100 --smoothLength 500

classify_reads.py \
	-v /proj/pmdvlab/data/variants/v5.snps.CC8.vcf.gz \
	-p C57BL_6NJ CAST_EiJ \
	-r chrX \
	-b aln/BF_1.chip.bam \
	-o allele_specific/H3K4me3/BF.u.bam && \
samtools sort allele_specific/H3K4me3/BF.u.bam >allele_specific/H3K4me3/BF.bam && \
samtools index allele_specific/H3K4me3/BF.bam

reads_by_tag.py -t XX -v m allele_specific/H3K4me3/BF.bam | \
samtools view -bhS >allele_specific/H3K4me3/BF.m.bam && \
samtools index allele_specific/H3K4me3/BF.m.bam

reads_by_tag.py -t XX -v p allele_specific/H3K4me3/BF.bam | \
samtools view -bhS >allele_specific/H3K4me3/BF.p.bam && \
samtools index allele_specific/H3K4me3/BF.p.bam

bamCoverage -b allele_specific/H3K4me3/BF.m.bam -of bigwig -o allele_specific/H3K4me3/BF.m.bw -r chrX:0:171031299 -bs 100 --smoothLength 500
bamCoverage -b allele_specific/H3K4me3/BF.p.bam -of bigwig -o allele_specific/H3K4me3/BF.p.bw -r chrX:0:171031299 -bs 100 --smoothLength 500


classify_reads.py \
	-v /proj/pmdvlab/data/variants/v5.snps.CC8.vcf.gz \
	-p CAST_EiJ C57BL_6NJ \
	-r chrX \
	-b ssds_bams/FB.bam \
	-o allele_specific/SSDS/FB.u.bam && \
samtools sort allele_specific/SSDS/FB.u.bam >allele_specific/SSDS/FB.bam && \
samtools index allele_specific/SSDS/FB.bam

reads_by_tag.py -t XX -v m allele_specific/SSDS/FB.bam | \
samtools view -bhS >allele_specific/SSDS/FB.m.bam && \
samtools index allele_specific/SSDS/FB.m.bam

reads_by_tag.py -t XX -v p allele_specific/SSDS/FB.bam | \
samtools view -bhS >allele_specific/SSDS/FB.p.bam && \
samtools index allele_specific/SSDS/FB.p.bam

bamCoverage -b allele_specific/SSDS/FB.m.bam -of bigwig -o allele_specific/SSDS/FB.m.bw -r chrX:0:171031299 -bs 100 --smoothLength 500
bamCoverage -b allele_specific/SSDS/FB.p.bam -of bigwig -o allele_specific/SSDS/FB.p.bw -r chrX:0:171031299 -bs 100 --smoothLength 500

classify_reads.py \
	-v /proj/pmdvlab/data/variants/v5.snps.CC8.vcf.gz \
	-p C57BL_6NJ CAST_EiJ \
	-r chrX \
	-b ssds_bams/BF.bam \
	-o allele_specific/SSDS/BF.u.bam && \
samtools sort allele_specific/SSDS/BF.u.bam >allele_specific/SSDS/BF.bam && \
samtools index allele_specific/SSDS/BF.bam

reads_by_tag.py -t XX -v m allele_specific/SSDS/BF.bam | \
samtools view -bhS >allele_specific/SSDS/BF.m.bam && \
samtools index allele_specific/SSDS/BF.m.bam

reads_by_tag.py -t XX -v p allele_specific/SSDS/BF.bam | \
samtools view -bhS >allele_specific/SSDS/BF.p.bam && \
samtools index allele_specific/SSDS/BF.p.bam

bamCoverage -b allele_specific/SSDS/BF.m.bam -of bigwig -o allele_specific/SSDS/BF.m.bw -r chrX:0:171031299 -bs 100 --smoothLength 500
bamCoverage -b allele_specific/SSDS/BF.p.bam -of bigwig -o allele_specific/SSDS/BF.p.bw -r chrX:0:171031299 -bs 100 --smoothLength 500
