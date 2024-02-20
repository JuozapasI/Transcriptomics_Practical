configfile: "config.yaml"

rule all:
    input:
        "multiqc_report.html",
        "multiqc_report_trimmed.html"

rule fastqc:
    input:
        "data/samples/{sample}.fastq.gz"
    output:
        "fastqc/{sample}_fastqc.html",
        "fastqc/{sample}_fastqc.zip"
    shell:
        "~/Kita/FastQC/fastqc -o fastqc/ {input}"

rule multiqc:
    input:
        expand("fastqc/{sample}_fastqc.html", sample=config['samples'])
    output:
        "multiqc_report.html"
    shell:
        "multiqc -f fastqc/"

rule trim:
    input:
        "data/samples/{sample}.fastq.gz"
    output:
        "samples_trimmed/{sample}_trimmed.fastq.gz"
    shell:
        "~/Kita/bbmap/bbduk.sh ref=data/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 in={input} out={output}"

rule trimmed_fastqc:
    input:
        "samples_trimmed/{sample}.fastq.gz"
    output:
        "trimmed_fastqc/{sample}_fastqc.html",
        "trimmed_fastqc/{sample}_fastqc.zip"
    shell:
        "~/Kita/FastQC/fastqc -o trimmed_fastqc/ {input}"

rule trimmed_multiqc:
    input:
        expand("trimmed_fastqc/{sample}_trimmed_fastqc.html", sample=config['samples'])
    output:
        "multiqc_report_trimmed.html"
    shell:
        "multiqc -f -n multiqc_report_trimmed.html trimmed_fastqc/"

rule star_genome_index:
    input:
        "data/chr19_20Mb.gtf",
        "data/chr19_20Mb.fa"
    output:
        directory("index/")
    shell:
        "STAR "
        "--runMode genomeGenerate "
        "--genomeSAindexNbases 11 "
        "--runThreadN 4 "
        "--genomeDir index/ "
        "--genomeFastaFiles data/chr19_20Mb.fa "
        "--sjdbGTFfile data/chr19_20Mb.gtf"
