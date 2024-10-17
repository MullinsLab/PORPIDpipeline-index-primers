import subprocess, sys

configfile: "2023-08_HIVDR_Run9_Pool_01.yaml"

DATASETS = [d for d in config for s in config[d]]
SAMPLES = [s for d in config for s in config[d]]
VERSION = "1.1"
COMMIT = subprocess.check_output(['git', 'rev-parse', '--verify', 'HEAD']).strip().decode()

sys.stderr.write("Running PorpidPostproc\n")
sys.stderr.write("Version: {0}\n".format(VERSION))
sys.stderr.write("Commit ID: {0}\n".format(COMMIT))

def contam_input(wildcards):
    SAMPLES = [s for s in config[wildcards.dataset]]
    return expand("porpid/{dataset}/consensus/{sample}.fasta",
        dataset = wildcards.dataset,
        sample = SAMPLES
    )

# PorpidPostproc parameters
# demux
chunk_size = 100000  		# default 100000
index_type = "Index_primer" # default "Index_Primer", also accepts "Nextera_primer"
error_rate = 0.01    		# default 0.01
min_length = 2000    		# default 2100
max_length = 4000    		# default 4000
#porpid
fs_thresh = 5        		# default 5
lda_thresh = 0.995   		# default 0.995
#contam
cluster_thresh = 0.015   	# default 0.015
proportion_thresh = 0.2  	# default 0.2
dist_thresh = 0.015      	# default 0.015
#postproc
agreement_thresh = 0.7   	# default 0.7
panel_thresh = 50        	# default 50
degap = "true"           	# default "true" # degapped postproc sequences


rule all:
    input:
        expand("porpid/{dataset}-porpid.tar.gz",
            dataset = DATASETS),
        expand("postproc/{dataset}-postproc.tar.gz",
            dataset = DATASETS)

rule demux:
    input:
        "raw-reads/{dataset}.fastq.gz"
    output:
        directory("porpid/{dataset}/demux"),
        "porpid/{dataset}/{dataset}_quality_report.csv",
        "porpid/{dataset}/{dataset}_demux_report.csv",
        "porpid/{dataset}/{dataset}_reject_report.csv"
    params:
        chunk_size = chunk_size,
        error_rate = error_rate,
        min_length = min_length,
        max_length = max_length,
        index_type = index_type,
        config = lambda wc: config[wc.dataset]
    script:
        "scripts/demux.jl"

rule porpid:
    input:
        "porpid/{dataset}/demux"
    output:
        directory("porpid/{dataset}/porpid/{sample}.fastq"),
        "porpid/{dataset}/tags/{sample}.csv"
    params:
        config = lambda wc: config[wc.dataset][wc.sample],
	    fs_thresh = fs_thresh,
        lda_thresh= lda_thresh
    script:
        "scripts/porpid.jl"

rule consensus:
    input:
        "porpid/{dataset}/porpid/{sample}.fastq"
    output:
        "porpid/{dataset}/consensus/{sample}.fasta"
    params:
        config = lambda wc: config[wc.dataset][wc.sample]
    script:
        "scripts/consensus.jl"

rule contam:
    input:
        files = contam_input,
        panel = "panels/contam_panel.fasta"
    output:
        directory("porpid/{dataset}/contam_passed"),
        directory("porpid/{dataset}/contam_failed"),
        "porpid/{dataset}/{dataset}_contam_report.csv",
        "porpid/{dataset}/{dataset}_contam_suspect.csv"
    params:
        proportion_thresh =proportion_thresh,
        cluster_thresh = cluster_thresh,
        dist_thresh = dist_thresh
    script:
        "scripts/contam.jl"

rule postproc:
    input:
        "porpid/{dataset}/contam_passed",
        "porpid/{dataset}/tags/{sample}.csv",
    output:
        report("postproc/{dataset}/{sample}/{sample}.fasta.mds.png", category = "postproc", caption = "report-rst/mds.rst"),
        "postproc/{dataset}/{sample}/{sample}.fasta.apobec.csv",
        report("postproc/{dataset}/{sample}/{sample}.fasta.tre.svg", category = "postproc", caption = "report-rst/highlighter.rst"),
        "postproc/{dataset}/{sample}/{sample}.fasta",
        report("postproc/{dataset}/{sample}/{sample}_qc_bins.png", category = "postproc", caption = "report-rst/bins.rst"),
        "postproc/{dataset}/{sample}/{sample}_qc_bins.csv",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.fasta",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.csv",
        "postproc/{dataset}/{dataset}_fasta/{sample}.fasta",
        "postproc/{dataset}/{dataset}_fasta/{sample}_consensus.fasta"
    params:
        panel = lambda wc: config[wc.dataset][wc.sample]["panel"],
        fs_thresh = fs_thresh,
        agreement_thresh = agreement_thresh,
        panel_thresh = panel_thresh,
        degap = degap,
    script:
        "scripts/postproc.jl"

rule report:
    input:
        "postproc/{dataset}/{sample}/{sample}_qc_bins.png",
        "postproc/{dataset}/{sample}/{sample}_qc_bins.csv",
        "postproc/{dataset}/{sample}/{sample}.fasta.mds.png",
        "postproc/{dataset}/{sample}/{sample}.fasta.tre.svg",
        "postproc/{dataset}/{sample}/{sample}.fasta",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.fasta",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.csv"
    params:
        VERSION = VERSION,
        COMMIT = COMMIT
    output:
        "postproc/{dataset}/{sample}/{sample}-report.html",
        "postproc/{dataset}/{sample}/{sample}-blast.html"
    script:
        "scripts/report.jl"
        
rule counts:
    input:
        expand("postproc/{dataset}/{sample}/{sample}-report.html", zip, dataset = DATASETS, sample = SAMPLES),
        expand("postproc/{dataset}/{sample}/{sample}_qc_bins.csv", zip, dataset = DATASETS, sample = SAMPLES)
    params:
        VERSION = VERSION,
        COMMIT = COMMIT,
        SAMPLES = SAMPLES,
        chunk_size = chunk_size,
        error_rate = error_rate,
        min_length = min_length,
        max_length = max_length,
        proportion_thresh =proportion_thresh,
        cluster_thresh = cluster_thresh,
        dist_thresh = dist_thresh,
        fs_thresh = fs_thresh, 
        lda_thresh= lda_thresh
    output:
        "postproc/{dataset}/{dataset}-counts.html",
        "postproc/{dataset}/{dataset}-blast.html",
        "postproc/{dataset}/{dataset}-seq-counts.csv"
    script:
        "scripts/counts.jl"        
               
rule tar:
    input:
        "postproc/{dataset}/{dataset}-counts.html"
    output:
        "porpid/{dataset}-porpid.tar.gz",
        "postproc/{dataset}-postproc.tar.gz"
    params:
        datasets = DATASETS,
        samples = SAMPLES
    script:
        "scripts/tar.jl"
