configfile: "config.yaml"

samples = config['samples']
samples_pairs = [x[:-7] for x in samples]
samples_pairs = list(set(samples_pairs))

rule all:
    input:
        "multiqc_report.html",
        "multiqc_report_trimmed.html",
        "fgsea/KAPA_pathways.tsv",
        "fgsea/Collibri_pathways.tsv",
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

def get_bam_samples(wildcards):
    return ["bam/" + string + ".bam" for string in samples_pairs if string.startswith(wildcards.sample)]

rule count_multi:
    input:
        get_bam_samples
    output:
        "counts/{sample}.txt.tmp1",
        "counts/{sample}.txt.tmp2",
        "counts/{sample}.txt.tmp1.summary",
        "counts/{sample}.txt.tmp2.summary"
    shell:
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id -a data/chr19_20Mb.gtf -o "
        "{output[0]} {input} -s 1 && "
        "~/Kita/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id -a data/chr19_20Mb.gtf -o "
        "{output[1]} {input} -s 2"

rule check_strand:
    input:
        "counts/{sample}.txt.tmp1",
        "counts/{sample}.txt.tmp2",
        "counts/{sample}.txt.tmp1.summary",
        "counts/{sample}.txt.tmp2.summary"
    output:
        "counts/{sample}.txt",
        "counts/{sample}.txt.summary"
    run:
        import pandas as pd
        import os
        df1 = pd.read_csv(input[2], sep="\t", index_col = 0)
        df2 = pd.read_csv(input[3], sep="\t", index_col = 0)
        row_comparison = (df1.iloc[0] > df2.iloc[0])
        if row_comparison.all():
            os.rename(input[0], output[0])
            os.rename(input[2], output[1])
        elif not row_comparison.any():
            os.rename(input[1], output[0])
            os.rename(input[3], output[1])
        else:
            print("Different strands were detected, check input files!")

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
        "counts/Collibri.txt",
        "counts/KAPA.txt",
        "deseq/Collibri.csv",
        "deseq/KAPA.csv"
    output:
        "pca/pca.png"
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
         "deseq/Collibri.png",
         "deseq/KAPA.png",
         "pca/pca.png",
         "practical.tex"
    output:
         "practical.pdf"
    shell:
         "pdflatex practical.tex && "
         "pdflatex practical.tex"
