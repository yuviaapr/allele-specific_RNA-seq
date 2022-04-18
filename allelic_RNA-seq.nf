#!/usr/bin/env nextflow

/*
 Gene-centered processing of allele-specific RNA-seq data
 Authors: 
	- Yuvia A. PEREZ RICO <yuvia.perez-rico@embl.de>
*/

log.info "        allele-specific RNA-seq - version 0.0.2        "
log.info "#######################################################"
log.info "Sample description file	= ${params.sampleInfo}"
log.info "Reads per split file		= ${params.chunkSize}"
log.info "Mouse strain 1		= ${params.G1}"
log.info "Mouse strain 2		= ${params.G2}"
log.info "Output directory		= ${params.outDir}"
log.info "Path to genome index		= ${params.genomeDirPath}"
log.info "Gene annotations		= ${params.annotations}"
log.info "SNP reference file		= ${params.snpFile}"
log.info "Script to calculate TPM	= ${params.rsTPM}"
log.info "Data from TX lines		= ${params.TXline}"
log.info "Number of threads		= ${params.numCPUs}"
log.info "\n"

/*
 Validate input parameters
*/

if( !(params.chunkSize instanceof Number) ){
	exit 1, "Invalid chunk size = ${params.chunkSize}"
}

if( !(params.numCPUs instanceof Number) ){
	exit 1, "Invalid number of CPUs = ${params.numCPUs}"
}

if( !(params.G1 in ['C57BL-6J', 'CAST-EiJ', 'PWK-PhJ', '129S1-SvImJ', 'FVB-NJ', '129P2-OlaHsd']) ){
	exit 1, "Invalid strain 1 name = ${params.G1}"
}

if( !(params.G2 in ['C57BL-6J', 'CAST-EiJ', 'PWK-PhJ', '129S1-SvImJ', 'FVB-NJ', '129P2-OlaHsd']) ){
	exit 1, "Invalid strain 2 name = ${params.G2}"
}

if( !(params.TXline in ['yes', 'no']) ){
	exit 1, "Invalid value to indicate if data was generated using TX lines = ${params.TXline}"
}

/*
 Validate input files
*/

sdFile = file(params.sampleInfo)
if( !sdFile.exists() ){
	exit 1, "The specified sample description file does not exist = ${params.sampleInfo}"
}
log.info "Checking sample description file = $sdFile"

featuresFile = file(params.annotations)
if( !featuresFile.exists() ){
	exit 1, "The specified feature annotation file does not exist = ${params.annotations}"
}
log.info "Checking feature annotation file = $featuresFile"

varFile = file(params.snpFile)
if( !varFile.exists() ){
	exit 1, "The specified SNP annotation file does not exist = ${params.snpFile}"
}
log.info "Checking SNP annotations file = $varFile"

calTPM = file(params.rsTPM)
if( !calTPM.exists() ){
	exit 1, "The specified R script file does not exist = ${params.rsTPM}"
}
log.info "Checking R script to calculate TPM values = $calTPM"

gnmDir = file(params.genomeDirPath)
if( !gnmDir.exists() ){
	exit 1, "The specified genome annotation directory does not exist = ${params.genomeDirPath}"
}
log.info "Checking genome annotation directory = $gnmDir"

resDir = file(params.outDir)
if( !resDir.exists() && !resDir.mkdirs() ){
	exit 1, "The specified directory to save results cannot be created = ${params.outDir}\n Check file system access permission"
}
log.info "Checking results directory = $resDir"

/*
 Program versions
*/

process get_program_versions{
	publishDir "${resDir}/software", mode: 'move'

	output:
	file('programs_version.txt') into programs_version

	"""
	echo nextflow ${nextflow.version} > tmp_version.txt
	fastqc -v >> tmp_version.txt
	echo trim_galore \$(trim_galore -v | awk '/version/{print\$2}') >> tmp_version.txt
	echo STAR \$(STAR --version) >> tmp_version.txt
	samtools --version | grep samtools >> tmp_version.txt
	echo SNPsplit \$(SNPsplit --version | awk '/Version/{print\$2}') >> tmp_version.txt
	picard SortSam --version 2>&1 | awk '{sub(/-SNAPSHOT/,"");print"picard "\$1}' >> tmp_version.txt
	featureCounts -v 2>&1 | grep feature >> tmp_version.txt
	bamCoverage --version >> tmp_version.txt
	sort tmp_version.txt > programs_version.txt
	"""

}

/*
 Create channels with the fastq files for processing and quality control
*/

// Reads1

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.sample, file(row.reads1)) }
	.set { samples_r1 }

// Reads2

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.sample, file(row.reads2)) }
	.set { samples_r2 }

// Both reads

Channel
	.fromPath(sdFile)
	.splitCsv(header:true)
	.map{ row -> tuple(row.sample, file(row.reads1), file(row.reads2)) }
	.set { samples_quality }

/*
 Step 0. Quality check
*/

process fastq_quality{
	publishDir "${resDir}/qc/fastqc", mode: 'copy'

	input:
	set val(name), file(reads1), file(reads2) from samples_quality

	output:
	file("${name}_fastqc") into fastqc_results

	"""
	mkdir ${name}_fastqc
	fastqc -o ${name}_fastqc -q -t ${params.numCPUs} ${reads1} ${reads2}
	"""

}

/*
 Step 1. Split fastq files
*/

process split_Reads1{
	label 'fastq_splitting'

	input:
	set val(name), file(reads) from samples_r1

	output:
	set val(name), file('*.gz') into samples_r1_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads1_"
	"""

}

process split_Reads2{
	label 'fastq_splitting'

	input:
	set val(name), file(reads) from samples_r2

	output:
	set val(name), file('*.gz') into samples_r2_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads2_"
	"""

}

samples_r1_split
	.combine(samples_r2_split, by: 0)
	.transpose()
	.set { samples }

/*
 Step 2. Trimming
*/

process trim_reads{
	publishDir "${resDir}/qc/trimming", mode: 'copy', pattern: "${name}/*report.txt"

	input:
	set val(name), file(reads1), file(reads2) from samples

	output:
	set val(name), file("${name}/*{1,2}.fq.gz") into trimmed_samples
	file("${name}/*report.txt") into trim_stats

	"""
	trim_galore --paired ${reads1} ${reads2} -o ${name}
	"""

}

/*
 Step 3. Mapping
*/

process read_mapping{
	publishDir "${resDir}/qc/mapping", mode: 'copy', pattern: '*Log.final.out'

	input:
	set val(name), file(trimmed_reads) from trimmed_samples

	output:
	set val(name), file('*.bam') into split_mapping
	file('*Log.final.out') into mapping_stats

	"""
	file=\$(echo ${trimmed_reads[0]} | sed s/_1_sequence.fastq.gz_val_1.fq.gz//)
	STAR --readFilesIn ${trimmed_reads} --outFileNamePrefix \$file. --outTmpDir ${params.tmpOutDir}/\$file \
	--readFilesCommand zcat --runThreadN ${params.numCPUs} --genomeDir ${gnmDir} --sjdbGTFfile ${featuresFile} \
	--sjdbOverhang 99 --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.06 \
	--alignIntronMax 500000 --alignMatesGapMax 500000 --alignEndsType EndToEnd \
	--outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted
	"""

}

/*
 Step 4. Remove singletons
*/

process singleton_filter{
	publishDir "${resDir}/qc/singleton_filter", mode: 'copy', pattern: '*stats'
	label 'filter_bams'

	input:
	set val(name), file(mapped_reads) from split_mapping

	output:
	set val(name), file('*.bam') into split_mapping_noSE
	file('*stats') into singleton_stats

	"""
	samtools flagstat ${mapped_reads} > ${mapped_reads.baseName}_rmSE.stats
	samtools view -b -f 0x2 ${mapped_reads} > ${mapped_reads.baseName}_rmSE.bam
	samtools flagstat ${mapped_reads.baseName}_rmSE.bam > ${mapped_reads.baseName}_rmSE_post.stats
	"""

}

/*
 Step 5. Remove mitochondrial reads
*/

process chrM_filter{
	publishDir "${resDir}/qc/chrM_filter", mode: 'copy', pattern: '*stats'
	label 'filter_bams'

	input:
	set val(name), file(PE_mapped_reads) from split_mapping_noSE

	output:
	set val(name), file('*.bam') into split_mapping_noMit
	file('*stats') into chrM_stats

	"""
	samtools view -h ${PE_mapped_reads} | grep -v chrM | samtools view -bS - > ${PE_mapped_reads.baseName}_rmChrM.bam
	samtools flagstat ${PE_mapped_reads.baseName}_rmChrM.bam > ${PE_mapped_reads.baseName}_rmChrM.stats
	"""

}

/*
 Step 6. Genotype assignation
*/

process SNPsplit{
	publishDir "${resDir}/qc/SNPsplit", mode: 'copy', pattern: "${name}/*.txt"

	input:
	set val(name), file(filtered_reads) from split_mapping_noMit

	output:
	set val(name), file("${name}/*genome1.bam") into split_genome1
	set val(name), file("${name}/*genome2.bam") into split_genome2
	set val(name), file("${name}/*unassigned.bam") into split_unassigned
	file("${name}/*report.txt") into SNPs_report_stats
	file("${name}/*sort.txt") into SNPs_sort_stats

	"""
	SNPsplit --paired --no_sort --snp_file ${varFile} -o ${name} ${filtered_reads}
	"""

}

/*
 Step 7. Merge and sort BAM files per sample
*/

// Organise files

split_genome1.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { genome1 }

split_genome2.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { genome2 }

split_unassigned.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { unassigned }

// Process G1 files

process sort_G1_bams{
	publishDir "${resDir}/mapping", mode: 'copy'
	label 'process_bams'

	input:
	set val(name), file(genome1_bam_files) from genome1

	output:
	set val(name), file('*sorted.bam') into G1_bam_merge, G1_bam

	"""
	# Merge
	samtools merge -f ${name}_${params.G1}.bam ${genome1_bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_${params.G1}.bam O=${name}_${params.G1}_sorted.bam SORT_ORDER=coordinate
	"""

}

// Process G2 files

process sort_G2_bams{
	publishDir "${resDir}/mapping", mode: 'copy'
	label 'process_bams'

	input:
	set val(name), file(genome2_bam_files) from genome2

	output:
	set val(name), file('*sorted.bam') into G2_bam_merge, G2_bam

	"""
	# Merge
	samtools merge -f ${name}_${params.G2}.bam ${genome2_bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_${params.G2}.bam O=${name}_${params.G2}_sorted.bam SORT_ORDER=coordinate
	"""

}

// Process unassigned files

process sort_unassigned_bams{
	input:
	set val(name), file(unassigned_bam_files) from unassigned

	output:
	set val(name), file('*sorted.bam') into unassigned_bam_merge

	"""
	# Merge
	samtools merge -f ${name}_unassigned.bam ${unassigned_bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_unassigned.bam O=${name}_unassigned_sorted.bam SORT_ORDER=coordinate
	"""

}

// Merge all reads

process merge_mapped{
	publishDir "${resDir}/mapping", mode: 'copy', pattern: '*bam'
	publishDir "${resDir}/qc/merged_bam", mode: 'copy', pattern: '*stats'

	input:
	set val(name), file(genome1_bam), file(genome2_bam), file(unassigned_bam) from G1_bam_merge.join(G2_bam_merge).join(unassigned_bam_merge)

	output:
	set val(name), file('*Gall_sorted.bam') into bam_sizeFactor, all_bam
	file('*stats') into merging_stats

	"""
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} MergeSamFiles MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true \
		I=${genome1_bam} I=${genome2_bam} I=${unassigned_bam} O=${name}_Gall.bam
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}_Gall.bam O=${name}_Gall_sorted.bam \
		SORT_ORDER=coordinate
	samtools flagstat ${name}_Gall_sorted.bam > ${name}_postMerging.stats
	"""

}

// Set combined channels

G1_bam.mix(G2_bam).mix(all_bam).into { bam_counts; bam_track_fwd; bam_track_rv }

/*
 Step 8. Read counting
*/

process read_counts{
	publishDir "${resDir}/quant", mode: 'copy', pattern: '*.counts'
	publishDir "${resDir}/qc/read_counting", mode: 'copy', pattern: '*.summary'
	label 'count_factor'

	input:
	set val(name), file(bam_fC) from bam_counts

	output:
	file("*.counts") into count_files
	file("*.summary") into counting_summary

	"""
	featureCounts -C -p -s 2 -t exon -g gene_id -a ${featuresFile} -o ${bam_fC.baseName}.counts ${bam_fC}
	"""

}

// Normalize read counts


process normalized_counts{
	publishDir "${resDir}/quant", mode: 'copy'
	label 'count_factor'

	input:
	file(raw_counts) from count_files.collect()

	output:
	file("*counts*.txt") into count_tables

	"""
	awk 'NR > 1 {print \$1"\t"\$6}' ${raw_counts[0]} > table_raw_counts.txt
	for i in ${raw_counts}; do cp table_raw_counts.txt tmp_table.txt; \
		awk 'NR > 1 {print\$7}' \$i | paste tmp_table.txt - > table_raw_counts.txt; \
		rm tmp_table.txt; done
	Rscript ${calTPM} table_raw_counts.txt
	"""

}

/*
 Step 9. Generate normalized signal tracks
*/

// Calculate size factors

process size_factors{
	label 'count_factor'

	input:
	set val(name), file(all_bam_SF) from bam_sizeFactor

	output:
	set val(name), stdout into size_factors_fwd, size_factors_rv

	"""
	samtools view -c ${all_bam_SF} | awk '{print 100000000/\$1}'
	"""

}

// Generate tracks

process signal_tracks_fwd{
	publishDir "${resDir}/signal", mode: 'copy'
	label 'signal_tracks'

	input:
	set val(name), file(bam_cov), val(sizeFactor) from bam_track_fwd.combine(size_factors_fwd, by: 0)

	output:
	file("*.bw") into bigWig_files_fwd

	script:
	if(params.TXline == "yes") {
		"""
		samtools index ${bam_cov}
		bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_fwd.bw --binSize 1 --normalizeUsing CPM \
			--filterRNAstrand forward --numberOfProcessors ${params.numCPUs_Dtools} --outFileFormat bigwig \
			 --ignoreForNormalization chr1 chr8 chr16 chrX --scaleFactor ${sizeFactor}
		"""
	} else if(params.TXline == "no") {
		"""
		samtools index ${bam_cov}
		bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_fwd.bw --binSize 1 --normalizeUsing CPM \
			--filterRNAstrand forward --numberOfProcessors ${params.numCPUs_Dtools} --outFileFormat bigwig --scaleFactor ${sizeFactor}
		"""
	}

}

process signal_tracks_rv{
	publishDir "${resDir}/signal", mode: 'copy'
	label 'signal_tracks'

	input:
	set val(name), file(bam_cov), val(sizeFactor) from bam_track_rv.combine(size_factors_rv, by: 0)

	output:
	file("*.bw") into bigWig_files_rv

	script:
	if(params.TXline == "yes") {
		"""
		samtools index ${bam_cov}
		bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_rv.bw --binSize 1 --normalizeUsing CPM \
			--filterRNAstrand reverse --numberOfProcessors ${params.numCPUs_Dtools} --outFileFormat bigwig \
			 --ignoreForNormalization chr1 chr8 chr16 chrX --scaleFactor ${sizeFactor}
		"""
	} else if(params.TXline == "no") {
		"""
		samtools index ${bam_cov}
		bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_rv.bw --binSize 1 --normalizeUsing CPM \
			--filterRNAstrand reverse --numberOfProcessors ${params.numCPUs_Dtools} --outFileFormat bigwig --scaleFactor ${sizeFactor}
		"""
	}

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Main results are saved in ${resDir}\n" : "There was an error during the execution, check log files." )
}


