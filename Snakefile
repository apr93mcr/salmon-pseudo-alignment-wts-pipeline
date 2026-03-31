
#salmon pseudo aligment gene expression quantification wts pipeline

out_dir=config["output_dir"]
sample_sheet = config["sample_sheet"]
salmon_index = config["salmon_index"]
gene_quant_script = config["gene_quant_script"]


import pandas as pd


data_list=pd.read_csv(sample_sheet)
print("Columns in the CSV:", data_list.columns)

sample_dict = {
    row["Samples"]: { # Keys are sample names
        "folder": row["Folder"],
        "fastq_r1": row["Folder"] + row["R1"],
        "fastq_r2": row["Folder"] + row["R2"],
    }
    for _, row in data_list.iterrows()
}

# Define all samples
SAMPLES = list(sample_dict.keys())
print(path_to_data)
print(len(SAMPLES))
print(out_dir)

rule all:
    input:
        expand(out_dir+"trimmed/{sample}_R1.trimmed.fastq.gz",sample=SAMPLES),
        expand(out_dir+"trimmed/{sample}_R2.trimmed.fastq.gz",sample=SAMPLES),
        expand(out_dir+"/quants/{sample}_ISR/quant.sf",sample=SAMPLES),
        expand(out_dir + "/quants/{sample}_ISR/{sample}.txt", sample=SAMPLES),
        expand(out_dir+"/quants/{sample}_ISR/quant.gene.sf",sample=SAMPLES)


rule fastp_pe:
    input:
        r1=lambda wildcards: sample_dict[wildcards.sample]["fastq_r1"],
        r2=lambda wildcards: sample_dict[wildcards.sample]["fastq_r2"]
    output:
        r1=out_dir+"trimmed/{sample}_R1.trimmed.fastq.gz",
        r2=out_dir+"trimmed/{sample}_R2.trimmed.fastq.gz",
        html=out_dir+"qc/fastp/{sample}.html",
        json=out_dir+"qc/fastp/{sample}.json"


    conda:
	"envs/RNAseq_norm"

    threads: 8
    
    priority:
        5

    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --qualified_quality_phred 20 \
            --length_required 35 \
            --json {output.json} \
            --html {output.html}
        """


rule salmon_stranded:
    input:
        r1=out_dir+"trimmed/{sample}_R1.trimmed.fastq.gz",
        r2=out_dir+"trimmed/{sample}_R2.trimmed.fastq.gz",
    
    output: 
        quant=out_dir+"/quants/{sample}_ISR/quant.sf",

    log:
        ".snakemake/log/salmon_stranded_log_{sample}"
    
    conda:
        "envs/salmon2"

    threads:
        20

    params:
        folder_output=out_dir+"/quants/{sample}_ISR"
	index=salmon_index
    priority:
        4

    shell:
        """
        salmon quant -i {params.index} -l ISR -1 {input.r1} -2 {input.r2} -p 20 -o {params.folder_output} --validateMappings --seqBias --gcBias >> {log} 2>> {log}
        """


rule quant_gene_stranded:
    input:
        quant=out_dir+"/quants/{sample}_ISR/quant.sf",
    
    output: 
        quant_genes=out_dir+"/quants/{sample}_ISR/quant.gene.sf"

    log:
        ".snakemake/log/gene_quant_{sample}"

    conda:
        "salmon2"

    params:
        folder_output=out_dir+"/quants/{sample}_ISR"
	gene_quant_script=gene_quant_script    
    priority:
        3

    shell:
        """
        Rscript {params.gene_quant_script} {params.folder_output} >> {log} 2>> {log}
        """


#snamake call comand:
#snakemake -s Scripts/Snakefile --configfile Scripts/stranded_config.yaml --jobs 10 --cores 80 --use-conda --conda-frontend conda --until salmon_stranded
