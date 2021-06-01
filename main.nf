#!/usr/bin/env nextflow

params.reads = "reads/*.fastq"
params.outdir = "results"
params.centrifuge_db = ""
params.cdhit_identity_threshold = 0.8
params.unicycler_mode = 'normal'
params.do_unicycler_assembly = false
params.do_canu_assembly = false
params.barcode_kit = "SQK-RBK004"

outdir = params.outdir

def helpMessage() {
    log.info"""
    =============================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    =============================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run ${workflow.manifest.name} --reads "/path/to/reads/*.fastq" --centrifuge_db /path/to/centrifuge/db/nt --outdir /output/path -profile standard
    
    Mandatory arguments:
      --reads                   Input reads directory and pattern (default: "${params.reads}")
    Other options:
      --outdir                      The output directory where the results will be saved (default: "${params.outdir}")
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      --do_unicycler_assembly       Run Unicycler assembly? (default: ${params.do_unicycler_assembly})
      --do_canu_assembly            Run CANU assembly? (default: ${params.do_canu_assembly})
      --centrifuge_db               Centrifuge taxonomic classification DB to use (e.g. NCBI nt DB) (default '${params.centrifuge_db}')
      --cdhit_identity_threshold    CD-HIT-EST nucleotide identity threshold for clustering (default ${params.cdhit_identity_threshold} (i.e. ${params.cdhit_identity_threshold * 100}% identity)).
      --unicycler_mode              Unicycler assembly mode [normal, conservative, bold] (default: ${params.unicycler_mode})
      --slurm_queue                 Name of SLURM queue to run workflow on (default ${params.slurm_queue})
      -profile                      Configuration profile to use. [standard, slurm] (default 'standard')
    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}


def check_centrifuge_db(centrifuge_db) {
  file_centrifuge_db = file(centrifuge_db)
  prefix = file_centrifuge_db.getName()
  centrifuge_dir = file_centrifuge_db.getParent()
  if ( !centrifuge_dir.isDirectory() || !centrifuge_dir.exists() ) {
    exit 1, "Centrifuge DB does not exist at '$centrifuge_dir'! Please specify a valid Centrifuge DB."
  }
  any_valid = false
  centrifuge_dir.eachFile { f ->
    if ( f.isFile() ) {
      if ( f.getName() =~ /^$prefix/ && f.getExtension() == 'cf') {
        any_valid = true
      }
    }
  }
  if ( !any_valid ) {
    exit 1, "No valid Centrifuge DB files with prefix '$prefix' in '$centrifuge_dir' and extension 'cf'! Please specify a valid Centrifuge classification DB directory and prefix."
  }
}

if (params.centrifuge_db != "") {
  check_centrifuge_db(params.centrifuge_db)
}

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  log.error "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
  exit 1
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info """=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Barcoding Kit'] = params.barcode_kit
summary['Centrifuge DB'] = params.centrifuge_db
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

outdir = params.outdir



// ==============
// PIPELINE START
// ==============
// get channel of reads from path

Channel
  .fromPath( params.reads )
  .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '/path/to/my-reads/'"}
  .into { ch_basecalled_reads }

// Download all publicly available NCBI Influenza sequences from the NCBI FTP 
// site. Genome metadata also downloaded for report generation.
process download_influenza_db {
  output:
  file('influenza.fna') into ch_influenza_db_fasta
  file('genomeset.dat') into ch_influenza_db_metadata

  """
  wget "ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz" 2> /dev/null
  gunzip influenza.fna.gz
  wget "ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat.gz" 2> /dev/null
  gunzip genomeset.dat
  """
}

// Create a nucleotide BLAST database from the Influenza sequences downloaded 
// from NCBI.
process blast_db {
  input:
  file('influenza.fna') from ch_influenza_db_fasta

  output:
  file("*.{nhr,nin,nsq}") into influenza_blastn_db

  """
  makeblastdb -dbtype nucl -in influenza.fna
  """
}

process split_segments {
  input:
  file('influenza.fna') from ch_influenza_db_fasta
  file('genomeset.dat') from ch_influenza_db_metadata

  output:
  file('*.fasta') into ch_influenza_db_segments

  script:
  """
  split_segments.py -s influenza.fna -o ./ -m genomeset.dat
  """
}

ch_influenza_db_segments
  .flatMap()
  .map {
    seg = file(it).getBaseName()
    [seg, it]
  }
  .set { ch_segment }

process cdhit_segments {
  tag "$seg"
  cpus 8

  input:
  set val(seg), file(fasta) from ch_segment
  output:
  set val(seg), file(fasta), file("${seg}-cdhit.fasta") into ch_cdhit_segment

  script:
  """
  cd-hit-est -i $fasta -o ${seg}-cdhit.fasta -c 0.95 -d 0 -M 0 -T ${task.cpus}
  """
}

process cdhit_cluster_reference_genomes {
  conda 'bioconda::cd-hit=4.6.8'
  cpus 56

  input:
    file 'influenza.fna' from ch_influenza_db_fasta
  output:
    file 'influenza-cdhit.fna' into influenza_db_fasta_cdhit

  script:
  """
  cd-hit-est -i influenza.fna -o influenza-cdhit.fna -c ${params.cdhit_identity_threshold} -d 0 -M ${task.memory.toMega()} -T ${task.cpus}
  """
}

process guppy_barcoder {
  tag "$fastq_name"

  input:
    file(fastq) from ch_basecalled_reads

  output:
    file('**/*.fastq') into ch_reads

  script:
  fastq_name = fastq.getBaseName()
  """
  guppy_barcoder -i ./ -s ./ --barcode_kits  -q 0
  ls -lht **/*.fastq
  """
}

ch_reads
  .flatMap()
  .map { 
    p = file(it).getParent()
    [p.getName(), it]
  }
  .groupTuple()
  .set { ch_barcoded_reads }

process cat_barcoded_fastq {
  tag "$barcode"

  input:
    set val(barcode), file('reads*.fastq') from ch_barcoded_reads
  output:
    set val(barcode), file("${barcode}.fastq") into ch_concat_barcoded_reads

  """
  cat reads*.fastq > ${barcode}.fastq
  """
}

process porechop {
  tag "$barcode"
  conda 'bioconda::porechop'
  publishDir "$outdir/reads", mode: 'copy', saveAs: { "${barcode}.fastq" }
  cpus 6

  input:
    set val(barcode), file(fastq) from ch_concat_barcoded_reads
  output:
    set val(barcode), file('porechop.fastq') into ch_porechop_reads, ch_reads_for_unicycler, ch_reads_for_centrifuge, ch_reads_for_canu, ch_reads_to_combine_with_segments

  """
  porechop -i $fastq -o porechop.fastq -t ${task.cpus}
  """
}

process map_against_ref {
  tag "$barcode"
  conda 'bioconda::minimap2 bioconda::samtools'
  cpus 6

  input:
    set val(barcode), file(fastq) from ch_porechop_reads
    file(refs) from influenza_db_fasta_cdhit
  output:
    set val(barcode), file("${barcode}.bam") into ch_bam_alignments

  script:
  """
  minimap2 -ax map-ont -t${task.cpus} $refs $fastq | samtools sort -@${task.cpus} > ${barcode}.bam
  """
}

ch_reads_to_combine_with_segments
  .combine(ch_cdhit_segment)
  .set { ch_reads_segment }

process map_against_ref_segment {
  tag "$barcode VS segment $seg"

  input:
    set val(barcode), file(fastq), val(seg), file(original_fasta), file(seg_fasta) from ch_reads_segment
  output:
    set val(barcode), val(seg), file(fastq), file(seg_fasta), file("${barcode}-${seg}.bam") into ch_segment_alignments

  script:
  """
  minimap2 -ax map-ont -t${task.cpus} $seg_fasta $fastq | samtools sort -@${task.cpus} > ${barcode}-${seg}.bam
  """
}

process samtools_stats {
  tag "$barcode VS segment $seg"
  publishDir "$outdir/mapping/$barcode/$seg", mode: 'copy', pattern: "*.{tsv,flagstat,idxstats}"

  input:
    set val(barcode), val(seg), file(fastq), file(seg_fasta), file(bam) from ch_segment_alignments

  output:
    set val(barcode), val(seg), file(fastq), file(seg_fasta), file(depths), file(flagstat), file(idxstats) into ch_segment_alignment_stats
  script:
  depths = "${barcode}-${seg}-depths.tsv"
  flagstat = "${barcode}-${seg}.flagstat"
  idxstats = "${barcode}-${seg}.idxstats"
  """
  echo -e "genome\tposition\tdepth" > $depths
  samtools depth -a -d 0 $bam >> $depths
  samtools flagstat $bam > $flagstat
  samtools idxstats $bam > $idxstats
  """
}

// Filter for alignments that did have some reads mapping to each segment
ch_segment_alignment_stats
  .filter {
    depth_linecount = file(it[4]).readLines().size()
    if (depth_linecount == 1) {
      println "No reads from \"${it[0]}\" mapped to segment ${it[1]}"
    }
    depth_linecount > 2
  }
  .set { ch_filtered_segment_alignment_stats }


process process_samtools_stats {
  tag "$barcode VS segment $seg"
  publishDir "$outdir/mapping/$barcode/$seg", mode: 'copy', pattern: "*.tsv"
  publishDir "$outdir/mapping/$barcode/$seg/top_reference", mode: 'copy', pattern: "*.fasta"

  input:
    set val(barcode), val(seg), file(fastq), file(seg_fasta), file(depths), file(flagstat), file(idxstats) from ch_filtered_segment_alignment_stats
    file('genomeset.dat') from ch_influenza_db_metadata

  output:
    set val(barcode), val(seg), file(fastq), file(output_ref_fasta), file(output_cov_stats) into ch_processed_mapping
  script:
  output_ref_fasta = "${barcode}-${seg}-top_reference.fasta"
  output_cov_stats = "${barcode}-${seg}-mapping_summary_stats.tsv"
  """
  process_mapping_stats_get_top_ref.py \
    -s $seg_fasta \
    -m genomeset.dat \
    -d $depths \
    -i $idxstats \
    -o $output_ref_fasta \
    -c $output_cov_stats
  """
}

process map_against_top_ref_segment {
  tag "$barcode VS segment $seg"
  conda 'bioconda::minimap2 bioconda::samtools'
  cpus 4

  input:
    set val(barcode), val(seg), file(fastq), file(output_ref_fasta), file(output_cov_stats) from ch_processed_mapping
  output:
    set val(barcode), val(seg), file(fastq), file(output_ref_fasta), file(output_cov_stats), file("${barcode}-${seg}.bam") into ch_top_segment_alignment, top_segment_alignment_ch2

  script:
  """
  minimap2 -ax map-ont -t${task.cpus} $output_ref_fasta $fastq | samtools sort -@${task.cpus} > ${barcode}-${seg}.bam
  """
}

process coverage_top_ref_segment {
  tag "$barcode VS segment $seg"
  conda 'bioconda::samtools'
  publishDir "$outdir/mapping/$barcode/$seg/top_reference", mode: 'copy', pattern: "*-depths.tsv"

  input:
    set val(barcode), val(seg), file(fastq), file(ref_fasta), file(cov_stats), file(bam) from ch_top_segment_alignment
  output:
    set val(barcode), val(seg), file(fastq), file(ref_fasta), file(cov_stats), file(bam), file(depths) into ch_top_segment_alignment_coverage

  script:
  depths = "${barcode}-${seg}-depths.tsv"
  """
  echo -e "genome\tposition\tdepth" > $depths
  samtools depth -a -d 0 $bam >> $depths
  """
}

process bcftools_mpileup {
  tag "$barcode VS segment $seg"
  conda 'bioconda::bcftools'

  input:
    set val(barcode), val(seg), file(fastq), file(ref_fasta), file(cov_stats), file(bam), file(depths) from ch_top_segment_alignment_coverage
  output:
    set val(barcode), val(seg), file(ref_fasta), file(depths), file(mpileup) into ch_mpileup
  script:
  mpileup = "${barcode}-${seg}.mpileup"
  """
  bcftools mpileup -f $ref_fasta $bam -Q 3 > $mpileup
  """
}

process bcftools_call {
  tag "$barcode VS segment $seg"
  conda 'bioconda::bcftools'
  publishDir "$outdir/mapping/$barcode/$seg/top_reference", mode: 'copy', pattern: '*.vcf.gz'
  cpus 4

  input:
    set val(barcode), val(seg), file(ref_fasta), file(depths), file(mpileup) from ch_mpileup
  output:
    set val(barcode), val(seg), file(ref_fasta), file(depths), file(vcf) into ch_vcf

  script:
  vcf = "${barcode}-${seg}.vcf.gz"
  """
  bcftools call -Oz --threads ${task.cpus} -mv $mpileup > $vcf
  """
}

process bcftools_consensus {
  tag "$barcode VS segment $seg"
  conda 'bioconda::bcftools'
  publishDir "$outdir/mapping/$barcode/$seg/top_reference", mode: 'copy', pattern: '*-consensus.fasta'

  input:
    set val(barcode), val(seg), file(ref_fasta), file(depths), file(vcf) from ch_vcf
  output:
    set val(barcode), val(seg), file(ref_fasta), file(depths), file(vcf), file(consensus) into ch_consensus

  script:
  consensus = "${barcode}-${seg}-consensus.fasta"
  """
  bcftools index $vcf
  bcftools consensus -f $ref_fasta $vcf > $consensus
  """
}

process fixed_consensus {
  tag "$barcode VS segment $seg"
  publishDir "$outdir/mapping/$barcode/$seg/top_reference", mode: 'copy', pattern: "*-fixed-consensus.fasta"

  input:
    set val(barcode), val(seg), file(ref_fasta), file(depths), file(vcf), file(consensus) from ch_consensus
  output:
    set val(barcode), val(seg), file(fixed_consensus) into ch_fixed_consensus

  script:
  fixed_consensus = "${barcode}-${seg}-fixed-consensus.fasta"
  output_cov_stats = "${barcode}-${seg}-mapping_summary_stats.tsv"
  low_cov_char = "N"
  cov_threshold = 1
  sample_name = "${barcode}-${seg}"
  """
  fix_consensus.py \
    -s $ref_fasta \
    -d $depths \
    -o $fixed_consensus \
    --cov-threshold $cov_threshold \
    --low-cov-char $low_cov_char \
    --sample-name "$sample_name"
  """
}

ch_fixed_consensus
  .groupTuple()
  .dump(tag: 'ch_fixed_consensus')
  .set {ch_grouped_consensus}

// Nucleotide BLAST all segment consensus sequences for each sample against the 
// NCBI Influenza sequences.
process blastn_consensus_seqs {
  tag "$barcode"
  cpus 8

  input:
  file(blastn_db) from influenza_blastn_db
  set val(barcode), val(seg), file(fasta) from ch_grouped_consensus

  output:
  set val(barcode), val(seg), file(blast_out) into ch_blastn

  script:
  db_name = blastn_db[0].baseName
  blast_out = "${barcode}-vs-ncbi-influenza.tsv"
  """
  cat $fasta > segs.fa
  blastn \
    -db $db_name \
    -query segs.fa \
    -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle" \
    -out $blast_out \
    -num_threads ${task.cpus} \
    -num_alignments 1000000 \
    -evalue 1e-6
  """
}

ch_blastn
  .collect { it[2] }
  .dump(tag: 'ch_all_blastn_results')
  .set { ch_all_blastn_results }

process subtyping_report {
  publishDir "$outdir/", mode: 'copy'

  input:
  file('genomeset.dat') from ch_influenza_db_metadata
  file(blastn_results) from ch_all_blastn_results

  output:
  file('subtyping_report.xlsx') into ch_subtyping_report_xlsx

  script:
  """
  parse_influenza_blast_results.py \
    --threads ${task.cpus} \
    --flu-metadata genomeset.dat \
    --excel-report subtyping_report.xlsx \
    $blastn_results
  """
}

if (params.do_unicycler_assembly) {
  process unicycler_assembly {
    tag "$barcode"
    publishDir "$outdir/assemblies/unicycler/$barcode", mode: 'copy'
    errorStrategy 'ignore'

    input:
      set val(barcode), file(fastq) from ch_reads_for_unicycler
    output:
      set val(barcode), val('unicycler'), file(output_contigs), file(output_gfa), file(output_unicycler_log) into ch_unicycler_assembly

    script:
    output_contigs = "${barcode}-assembly.fasta"
    output_gfa = "${barcode}-assembly.gfa"
    output_unicycler_log = "${barcode}-unicycler.log"
    """
    unicycler -t ${task.cpus} \
      --mode ${params.unicycler_mode} \
      -o $barcode \
      -l $fastq \
      --min_component_size 100 \
      --min_dead_end_size 100 \
      --linear_seqs 8
    ln -s ${barcode}/assembly.fasta $output_contigs
    ln -s ${barcode}/assembly.gfa $output_gfa
    ln -s ${barcode}/unicycler.log $output_unicycler_log
    """
  }
}

if (params.do_canu_assembly){
  process canu_assembly {
    tag "$barcode"
    publishDir "$outdir/assemblies/canu/$barcode", mode: 'copy'
    errorStrategy 'ignore'
    conda 'bioconda::canu'

    input:
      set val(barcode), file(fastq) from ch_reads_for_canu
    output:
      set val(barcode), val('canu'), file(contigs), file(gfa), file(report) optional true into ch_canu_assembly
    when:
      fastq.size() >= 10000

    script:
    contigs = "${barcode}.contigs.fasta"
    gfa = "${barcode}.contigs.gfa"
    report = "${barcode}.report"
    """
    canu -p $barcode -d $barcode minReadLength=100 minOverlapLength=50 stopOnReadQuality=false genomeSize=2k -nanopore-raw $fastq
    ln -s $barcode/${barcode}.contigs.fasta $contigs
    ln -s $barcode/${barcode}.contigs.gfa $gfa
    ln -s $barcode/${barcode}.report $report
    """
  }
}

if (params.centrifuge_db != "") {
  Channel.value( 
    [ 
      file(params.centrifuge_db).getName(), 
      file(params.centrifuge_db).getParent() 
    ] )
    .set { ch_centrifuge_db }

  process centrifuge {
    tag "$barcode"
    publishDir "$outdir/centrifuge/$barcode", pattern: "*.tsv", mode: 'copy'

    input:
      set db_name, file(centrifuge_db_dir) from ch_centrifuge_db
      set val(barcode), file(reads) from ch_reads_for_centrifuge
    output:
      set val(barcode), file(reads), file(results), file(report) into ch_centrifuge_results

    script:
    results = "${barcode}-centrifuge_results.tsv"
    report = "${barcode}-report.tsv"
    """
    centrifuge -x ${centrifuge_db_dir}/${db_name} -U $reads -S $results --report-file $report --mm
    """
  }

  process centrifuge_kraken_report {
    tag "$barcode"
    publishDir "$outdir/centrifuge", pattern: "*-kreport.tsv", mode: 'copy'

    input:
      set db_name, file(centrifuge_db_dir) from ch_centrifuge_db
      set val(barcode), file(reads), file(results), file(report) from ch_centrifuge_results
    output:
      set val(barcode), file(reads), file(results), file(kreport) into ch_centrifuge_kraken_report_results

    script:
    kreport = "${barcode}-kreport.tsv"
    """
    centrifuge-kreport -x ${centrifuge_db_dir}/${db_name} $results > $kreport
    """
  }
}


/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
