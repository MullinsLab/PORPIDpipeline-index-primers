ENV["MPLBACKEND"] = "Agg"
using Base64, CSV, DataFrames, NextGenSeqUtils, FASTX, WebBlast


function get_image_str(file)
    fig_str = open(file,"r") do io
        read(io)
    end
    return Base64.base64encode(fig_str)
end

function format_tbl(df)
    # df = CSV.read(file, DataFrame)
    tbl_html = DataFrames.repr("text/html", df)
    tbl_fmt = replace(tbl_html,
        "<table class=\"data-frame\">" => "<table class=\"table table-striped\";>")
    tbl_fmt = replace(tbl_fmt, r"<p>.+columns</p>" => "") #remove display description
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String15\">String15</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Int64\">Int64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Int64\">Int64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Float64\">Float64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<tr><th></th></tr>" => "") #remove type rows
    # tbl_fmt = replace(tbl_fmt, "<tr><th></th><th title=\"String\">String</th><th title=\"Int64\">Int64</th></tr>" => "") #remove types
    # tbl_fmt = replace(tbl_fmt, r"<thead>.+</thead>" => "") #remove heading row
    tbl_fmt = replace(tbl_fmt, r"<th></th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th>.</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th>..</th>" => "") #remove row count column
    return tbl_fmt
end



sample = snakemake.wildcards["sample"]
dataset = snakemake.wildcards["dataset"]

qc_fig = get_image_str(snakemake.input[1]);
qc_tbl = format_tbl(CSV.read(snakemake.input[2], DataFrame));
mds_fig = get_image_str(snakemake.input[3]);
phylo_fig = get_image_str(snakemake.input[4]);
pr_tbl = format_tbl(CSV.read(snakemake.input[7], DataFrame));

br_tbl = "click <a href=$(sample)-blast.html> here </a> for blast report"

html_str_hdr = """
<html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
        <style>
            body{ margin:0 100; background:white; }
        </style>
        <style type="text/css">
            .container{
            display: flex;
            }
            .col{
                flex: 1;
                text-align: center;
                align-self: center;
            }
            img {
                margin: auto;
            }
        </style>
    </head>
"""
html_str = html_str_hdr * """
    <body>
        <div style="max-width: 700px;">
            <h2>PorpidPostproc report</h2>
            <h3><i>Sample: $(sample)</i></h3>
            <h3><i>Dataset: $(dataset)</i></h3>
            <h4><i>Version: $(snakemake.params["VERSION"])</i></h3>
            <h4><i>Commit ID: $(snakemake.params["COMMIT"])</i></h3>
            <h3>
                Results of PORPID processing
            </h3>
            <h4>UMI stripplot</h4>
            <p>
                The stripplot below displays UMI family size (the number of individual
                CCS reads with a given UMI) as a function of the UMI length determined
                by PORPID. Only "likely_real" UMIs are kept for downstream analysis.
                "LDA-rejects" are likely offspring from other UMI bins, "fs<5" indicates
                sequencing depth was under 5 CCS and too low for consensus analysis.
                UMIs were flagged as "heteroduplex" when the set of reads had a signature
                drop in quality in the UMI region indicative of a superimposed signal of
                two different UMI sequences during circular consensus generation.
                Some UMIs of length other than 8 pass all other criteria, but are still
                excluded from analysis ("UMI_len != 8").
            </p>
            <img src="data:image/png;base64,$(qc_fig)
            "alt="" width=100% style="max-width: 700px;">
            <h4>UMI family and CCS totals for each type of PORPID result</h4>
            <p>
               Note that any "BPB-rejects" in the table below are sequences that were
               discarded due to <b> bad primer blocks </b>.
               These sequences were not included in the stripplot above as they
               were discarded before the breakdown into UMI families.
               Sequences rejected by the <tt>porpid</tt> script can be recovered by
               looking in the <tt>porpid</tt> archive for the file: <br>
               <tt> porpid/$(dataset)/porpid/$(sample).fastq/$(sample)_rejects.fastq </tt>
            </p>
            $(qc_tbl)
            <h4>
               Panel file check
            </h4>
            <p>
               Further sequences may have been rejected if they are deemed to be NOT close enough
               to the profile of the collection of sequences stored in the <b>panel</b> file.
               If the table below displays a rejection count greater than zero then the
               blast report should be consulted to provide further insight.
            </p>
            $(pr_tbl)
            <h3>
            Post-processing of single-template consensus sequences
            </h3>
            <h4>Multimensional scaling (MDS)  of template sequences and G>A APOBEC model</h4>
            <p>
                Classical MDS was used to represent all pairwise distances between template
                sequences in 2D space. Individual points are scaled by their family size
                (see inset key). A Bayesian model for APOBEC hypermutation is run on the
                single nucleotide substitution matrix between the global consensus and
                query sequence, and estimates the overall mutation rate and a G>A
                accelerator parameter from the data. Points are colored by probability that
                the G>A accelerator parameter is greater than 1.
            <p>
            <img src="data:image/png;base64,$(mds_fig)
            "alt="" width=100% style="max-width: 700px;">
            <h4>Phylogeny and highlighter plot of collapsed nucleotide variants</h4>
            <p>
                After collapsing by nucleotide sequence identity, variants are numbered in
                descending order. The most frequent variant is used as the master sequence for
                the highlighter plot and the root for the phylogeny. Any variant above 10% of the
                population is colored in red on the tree. If all sequences are identical, 
                a collapsed tree/highlighter figure cannot be created. 
            <p>
            <img src="data:image/svg+xml;base64,$(phylo_fig)
            "alt="" width=100% style="max-width: 700px;">
            <p>
            <h4> BLAST results for clades identified in this sample </h4>
            $(br_tbl)
        </div>
    </body>
 </html>
""";

open(snakemake.output[1],"w") do io
    write(io,html_str)
end

html_str = html_str_hdr * """
<body>
<h3> Blast results ... </h3>

No blast results available yet, to get blast results run

<code>
snakemake -s SnakeBlast -k --rerun-incomplete -j1
</code>
</body>
"""

open(snakemake.output[2],"w") do io
    write(io,html_str)
end
