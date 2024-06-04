version 1.0

workflow AlignFASTQ {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    File tar_star_reference
    File white_list
    Int chemistry
    String star_strand_mode
    String output_bam_basename
    # exercise 1: add in new input field for the workflow:
    String output_bam_experimentname
  }
  call STARsoloFastq {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      tar_star_reference = tar_star_reference,
      white_list = white_list,
      chemistry = chemistry,
      star_strand_mode = star_strand_mode,
      output_bam_basename = output_bam_basename,
      # exercise 1: add in new input for the task call:
      output_bam_experimentname = output_bam_experimentname
  }
  
  # exercise 3: add in a call to the new task to calculate UMIs per barcode
  call CalculateCounts { 
    input: 
        mtx_barcodes = STARsoloFastq.barcodes,     
        mtx_features = STARsoloFastq.features, 
        mtx_matrix = STARsoloFastq.matrix 
  }


  output {
    File barcodes = STARsoloFastq.barcodes
    File features = STARsoloFastq.features
    File matrix = STARsoloFastq.matrix
    # exercise 2: add the bam file as a new output
    File bam = STARsoloFastq.bam
    # exercise 3: add the UMI counts violin plot as a new output
    File violin_plot = CalculateCounts.violin_plot

  }
}

task STARsoloFastq {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    File tar_star_reference
    File white_list
    Int chemistry
    String star_strand_mode
    String output_bam_basename
    # exercise 1: add in new input for the task
    String output_bam_experimentname

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.1-2.7.11a-1692706072"
    Int machine_mem_gb = 40
    Int cpu = 8
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 3)) + ceil(size(r1_fastq, "Gi") * 20) +  ceil(size(r2_fastq, "Gi") * 20)
    # by default request non preemptible machine to make sure the slow star alignment step completes

  }

  meta {
    description: "Aligns reads in bam_input to the reference genome in tar_star_reference"
  }

  parameter_meta {
    r1_fastq: "input FASTQ file array"
    r2_fastq: "array of forward read FASTQ files"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    star_strand_mode: "STAR mode for handling stranded reads. Options are 'Forward', 'Reverse, or 'Unstranded'"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_gb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
  }

  command <<<
       set -e

    # Set 10x chemistry parameters
    UMILen=10
    CBLen=16
    if [ "~{chemistry}" == 2 ]
    then
        ## V2
        UMILen=10
        CBLen=16
    elif [ "~{chemistry}" == 3 ]
    then
        ## V3
        UMILen=12
        CBLen=16
    else
        echo Error: unknown chemistry value: "$chemistry". Should be one of "tenX_v2" or "texX_v3".
        exit 1;
    fi

    # Check that the star strand mode matches STARsolo aligner options
    if [[ "~{star_strand_mode}" == "Forward" ]] || [[ "~{star_strand_mode}" == "Reverse" ]] || [[ "~{star_strand_mode}" == "Unstranded" ]]
    then
        ## single cell or whole cell
        echo STAR mode is assigned
    else
        echo Error: unknown STAR strand mode: "~{star_strand_mode}". Should be Forward, Reverse, or Unstranded.
        exit 1;
    fi

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"

    # run STARsolo software
    STAR \
        --soloType Droplet \
        --soloStrand ~{star_strand_mode} \
        --runThreadN ~{cpu} \
        --genomeDir genome_reference \
        --readFilesIn "~{sep=',' r2_fastq}" "~{sep=',' r1_fastq}" \
        --readFilesCommand "gunzip -c" \
        --soloCBwhitelist ~{white_list} \
        --soloUMIlen $UMILen --soloCBlen $CBLen \
        --soloFeatures GeneFull_Ex50pAS \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30  \
        --soloCBmatchWLtype 1MM_multi \
        --soloUMIdedup 1MM_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes UB UR UY CR CB CY NH GX GN sF \
        --soloBarcodeReadLength 0 \
        --soloCellReadStats Standard \
        --soloUMIfiltering MultiGeneUMI_CR
      
    SoloDirectory="Solo.out/GeneFull_Ex50pAS/raw"
    echo "SoloDirectory is $SoloDirectory"
    find "$SoloDirectory" -maxdepth 1 -type f -name "*.mtx" -print0 | xargs -0 -I{}  echo mv {} /cromwell_root/
    find "$SoloDirectory" -maxdepth 1 -type f -name "*.mtx" -print0 | xargs -0 -I{} mv {} /cromwell_root/
    # Move files from the STARSolo output directory to the cromwell_root directory
    mv "Solo.out/GeneFull_Ex50pAS/raw/barcodes.tsv" barcodes.tsv
    mv "Solo.out/GeneFull_Ex50pAS/raw/features.tsv" features.tsv
    mv "Solo.out/GeneFull_Ex50pAS/CellReads.stats" CellReads.stats
    mv "Solo.out/GeneFull_Ex50pAS/Features.stats" Features.stats
    mv "Solo.out/GeneFull_Ex50pAS/Summary.csv" Summary.csv
    mv "Solo.out/GeneFull_Ex50pAS/UMIperCellSorted.txt" UMIperCellSorted.txt
    
    # exercise 2: update the bam's named to combine output_bam_basename and output_bam_experimentname
    mv Aligned.sortedByCoord.out.bam ~{output_bam_basename + output_bam_experimentname}.bam

  >>>

  runtime {
    docker: docker
    memory: "~{machine_mem_gb} GiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
  }

  output {
    File barcodes = "barcodes.tsv"
    File features = "features.tsv"
    File matrix = "matrix.mtx"
    # exercise 2: include the aligned bam file in the outputs
    File bam = "~{output_bam_basename + output_bam_experimentname}.bam"
  }
}

# exercise 3: add a new task to calculate UMIs per barcode
task CalculateCounts {
  input {
    File mtx_features
    File mtx_barcodes
    File mtx_matrix
  }
  command <<<
    echo printing directory
    pwd
    # make directory for the matrix files
    mkdir sample
    # copy the matrix files into the new directory
    cp ~{mtx_features} sample
    cp ~{mtx_barcodes} sample
    cp ~{mtx_matrix} sample
    gzip /cromwell_root/sample/*
    echo "Checking that gzipped files are in the sample directory"
    ls sample
    # Run python script and specify the directory path where matrix files are
    echo "Running script to calculate counts from matrix and print violin plot"
    python3 /data/umis_per_barcode.py /cromwell_root/sample/
    
    # Finished output file will be in the figures folder in cromwell_root. 
  >>>
  
  output {
    File violin_plot= "figures/violinviolin_plot.png"
  }
  runtime {
    docker: "ekiernan/umis_python:v2" # custom pre-built docker with the Python script used to count & plot UMIs per barcode
    memory: "5 GiB"
  }
}
