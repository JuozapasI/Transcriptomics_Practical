configfile: "config.yaml"

samples = config['samples']
samples_pairs = [x[:-7] for x in samples]
samples_pairs = list(set(samples_pairs))

rule all:
    input:
        "multiqc_report.html",
        "multiqc_report_trimmed.html",
        "counts/counts_incorrect_strand.txt",
        expand("counts/{sample}.txt", sample=samples_pairs),
        "practical.pdf"

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
        expand("fastqc/{sample}_fastqc.html", sample=samples)
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
        expand("trimmed_fastqc/{sample}_trimmed_fastqc.html", sample=samples)
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

rule star_map:
    input:
        "samples_trimmed/{sample}_R1_001_trimmed.fastq.gz",
        "samples_trimmed/{sample}_R2_001_trimmed.fastq.gz",
        "index/"
    output:
        "bam/{sample}.bam"
    shell:
        "STAR "
        "--runThreadN 4 "
        "--genomeDir index/ "
        "--readFilesIn samples_trimmed/{wildcards.sample}_R1_001_trimmed.fastq.gz samples_trimmed/{wildcards.sample}_R2_001_trimmed.fastq.gz "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix bam/ "
        "--readFilesCommand zcat && "
        "mv bam/Aligned.sortedByCoord.out.bam bam/{wildcards.sample}.bam"

rule index_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule count:
    input:
        "bam/{sample}.bam"
    output:
        "counts/{sample}.txt"
    shell:
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id -a data/chr19_20Mb.gtf -o "
        "counts/tmp1.txt {input} -s 1 && "
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id -a data/chr19_20Mb.gtf -o "
        "counts/tmp2.txt {input} -s 2 && " 
        "entry_file1=$(awk 'NR==2 {{print $2}}' counts/tmp1.txt.summary) && "
        "entry_file2=$(awk 'NR==2 {{print $2}}' counts/tmp2.txt.summary) && "
        "if ((entry_file1 > entry_file2)); then mv counts/tmp1.txt {output} && mv counts/tmp1.txt.summary {output}.summary; "
        "else mv counts/tmp2.txt {output} && mv counts/tmp2.txt.summary {output}.summary; fi"

rule count_collibri:
    input:
        expand("bam/{sample}.bam", sample=[x for x in samples_pairs if x.startswith("Collibri")])
    output:
        "counts/collibri.txt"
    shell:
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id "
        "-a data/chr19_20Mb.gtf -o {output} {input} -s 1"

rule count_kapa:
    input:
        expand("bam/{sample}.bam", sample=[x for x in samples_pairs if x.startswith("KAPA")])
    output:
        "counts/kapa.txt"
    shell:
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id "
        "-a data/chr19_20Mb.gtf -o {output} {input} -s 2"

rule count_incorrect_strand:
    output:
        "counts/counts_incorrect_strand.txt"
    shell:
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id "
        "-a data/chr19_20Mb.gtf -o counts/counts_incorrect_strand.txt "
        "bam/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001.bam -s 2"

rule DESeq:
    input:
        "counts/{sample}.txt"
    output:
        "deseq/{sample}.csv",
        "deseq/{sample}.png"
    shell:
        "Rscript scripts/deseq2.R {wildcards.sample}"

rule PCA:
    input:
        "counts/collibri.txt",
        "counts/kapa.txt",
        "deseq/collibri.csv",
        "deseq/kapa.csv"
    output:
        "pca/collibri.png",
        "pca/kapa.png"
    shell:
        "Rscript scripts/pca.R {input}"

rule fgsea:
    input:
         "deseq/{sample}.csv"
    output:
         "fgsea/{sample}_pathways.tsv"
    shell:
         "Rscript scripts/fgsea.R {wildcards.sample}"

rule reportPDF:
    input:
         "deseq/collibri.png",
         "deseq/kapa.png",
         "practical.tex"
    output:
         "practical.pdf"
    shell:
         "pdflatex practical.tex && "
         "pdflatex practical.tex"
