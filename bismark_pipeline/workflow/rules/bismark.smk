########################## BISMARK RULES

#copy reference file for BISMARK, cause BISMARK is stupid and cant put indexes anywhere else
rule copy_reference:
    input:
        os.path.join(config["reference_folder"],config["genome_fasta"])
    output:
        os.path.join("results/references",config["genome_fasta"])
    shell:
        "cp {input} {output}"

# likely will need - sudo apt-get install libtbb2
rule bismark_genome_preparation_fa_gz:
    input:
        os.path.join("results/references",config["genome_fasta"])
    output:
        directory("results/references/Bisulfite_Genome")   
    log:
        "logs/bismark/bismark_genome_preparation.log"
    params:
        "--verbose --bowtie2"  # optional params string
    threads: 
        6
    conda:
        os.path.join(workflow.basedir, "envs/bismark.yaml")
    script: 
        os.path.join(workflow.basedir, "scripts/bismark_genome_preparation_fa_gz.py")

rule bam2nuc_for_genome:
    input:
        genome_fa=os.path.join("results/references",config["genome_fasta"])
    output:
        report="results/references/genomic_nucleotide_frequencies.txt"
    log:
        "logs/bismark/bismark2nuc.log"
    threads: 
        6
    conda:
        os.path.join(workflow.basedir, "envs/bismark.yaml")
    script: 
        os.path.join(workflow.basedir, "scripts/bam2nuc_for_genome.py")


rule bismark_pe:
    input:
        fq_1="results/trimmed_reads/{SAMPLE}_R1_val_1.fq.gz",
        fq_2="results/trimmed_reads/{SAMPLE}_R2_val_2.fq.gz",
        genome=os.path.join("results/references",config["genome_fasta"]),
        bismark_indexes_dir="results/references/Bisulfite_Genome",
        genomic_freq="results/references/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark/{SAMPLE}/{SAMPLE}.bam"),
        report="results/reports/{SAMPLE}/{SAMPLE}_PE_report.txt",
        nucleotide_stats="results/reports/{SAMPLE}/{SAMPLE}.nucleotide_stats.txt"
        # bam_unmapped_1="bams/{SAMPLE}_{genome}_unmapped_reads_1.fq.gz",
        # bam_unmapped_2="bams/{SAMPLE}_{genome}_unmapped_reads_2.fq.gz",
        # ambiguous_1="bams/{SAMPLE}_{genome}_ambiguous_reads_1.fq.gz",
        # ambiguous_2="bams/{SAMPLE}_{genome}_ambiguous_reads_2.fq.gz"
    log:
        "logs/bismark/{SAMPLE}_bismark_aligning.log"
    threads: 
        3
    params:
        extra="--nucleotide_coverage",
        basename="{SAMPLE}"
    conda:
        os.path.join(workflow.basedir, "envs/bismark.yaml")
    script: 
        os.path.join(workflow.basedir, "scripts/bismark_pe.py")



rule bismark_methylation_extractor:
    input: 
        "results/bismark/{SAMPLE}/{SAMPLE}.bam"
    output:
        mbias_r1="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.M-bias_R1.png",
        mbias_r2="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.M-bias_R2.png",
        mbias_report="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.M-bias.txt",
        splitting_report="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}_splitting_report.txt",

            # coverafe file - 1-based chromosome coordinates ('inclusive') methylation info: % and counts        
        methylome_CpG_cov="results/bismark_methylation_not_deduplicated/{SAMPLE}/{SAMPLE}.bismark.cov.gz",
            # BedGraph with methylation percentage: 0-based start, 1-based end (exclusive)
        methylome_CpG_mlevel_bedGraph="results/bismark_methylation_not_deduplicated/{SAMPLE}/{SAMPLE}.bedGraph.gz",
            # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="results/bismark_methylation_not_deduplicated/{SAMPLE}/CpG_context_{SAMPLE}.txt.gz",
            # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="results/bismark_methylation_not_deduplicated/{SAMPLE}/CHG_context_{SAMPLE}.txt.gz",
        read_base_meth_state_chh="results/bismark_methylation_not_deduplicated/{SAMPLE}/CHH_context_{SAMPLE}.txt.gz"
        # genome wide cytosine report; contains strand information, ideal for bsses read.bismark(); 1-based chromosome coordinates 
        #cytosine_report="results/methylation/{SAMPLE}/{SAMPLE}.cytosine_report.txt.gz",
    log:
        "logs/bismark/{SAMPLE}_methylation_extraction_no_deduplication.log"
    threads: 
        3
    params:
        extra="--gzip --comprehensive --bedGraph --multicore 2"  # optional params string
    conda:
        os.path.join(workflow.basedir, "envs/bismark.yaml")
    script: 
        os.path.join(workflow.basedir, "scripts/bismark_methylation_extractor.py")

rule bismark2report_pe:
    input:    
        alignment_report="results/reports/{SAMPLE}/{SAMPLE}_PE_report.txt",
        nucleotide_report="results/reports/{SAMPLE}/{SAMPLE}.nucleotide_stats.txt",
        mbias_report="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.M-bias.txt",
        splitting_report="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}_splitting_report.txt"
    output:
        html="results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.bismark2report.html"
    log:
       "logs/bismark/{SAMPLE}_bismark2report_no_deduplication.log"
    params:
        skip_optional_reports=True
    conda:
        os.path.join(workflow.basedir, "envs/bismark.yaml")
    script: 
        os.path.join(workflow.basedir, "scripts/bismark2report_pe.py")