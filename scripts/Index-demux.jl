"""
This julia script takes an input CCS.fastq file along with a primer sheet specifying
the PCR primers used to amplify and each samples as well as the unique "index" primer combination
used to tag each sample and demultiplexes into individual .fastq files. For usage with
Snakemake pipeline.

This was file was originally written by Alec Pankow for use with an older version of this pipeline. 
It has been added to the updated porpid/postproc pipeline created by Hugh Murrell by Dylan Westfall
to enable the pipeline to accept samples labeled with Index primers rather than using unique
Sample IDs present in the cDNA primer.
"""

#load required packages
ENV["MPLBACKEND"] = "Agg"
using NextGenSeqUtils, DataFrames, DataFramesMeta, CSV,
StatsBase, IterTools, StringDistances, Statistics

input_fastq_path = snakemake.input[1]
filtered_path = snakemake.output["filtered_path"]
dataset = split(basename(input_fastq_path),".")[1]
error_rate = snakemake.params["error_rate"]
min_length = snakemake.params["min_length"]
max_length = snakemake.params["max_length"]
index_type = snakemake.params["index_type"]
#templates = snakemake.params["templates"] #Alec originally used "templates" but this snakemake config and pipeline uses "SAMPLES" instead of "templates" so all instances of "templates" were swapped to SAMPLES
SAMPLES = snakemake.params["config"]

mkdir(snakemake.output[1])

println("Performing Quality Filtering with the following parameters...")
println("Input fastq: $(input_fastq_path)")
println("Minimum library length: $(min_length)")
println("Maximum library length: $(max_length)")
println("Error rate: $(error_rate)")

#defining functions
function unique_not_substr(a)
    out = []
    for i in unique(a)
        res = true
        for j in unique(a)
            if occursin(i, j) & (i != j)
                res = false
            end
        end
        if res
            push!(out, i)
        end
    end
    return out
end

function iterative_primer_match(seqs,full_primers,window::Int,slide_by::Int;tol_one_error=true)
    if(slide_by + window - 1 > minimum(length.(full_primers)))
        @warn ("Matching window extends beyond shortest primer. This is ok, but check that you aren't matching something too short.")
    end
    primers = [p[1:min(window,minimum(length.(full_primers)))] for p in full_primers]
    filter = fast_primer_match(seqs,primers,tol_one_error=tol_one_error);
    for i in 2:slide_by
        unresolved = filter .== 0
        primers = [p[i:min(i + window - 1,minimum(length.(full_primers)))] for p in full_primers]
        filter[unresolved] = fast_primer_match(seqs[unresolved],primers,tol_one_error=tol_one_error);
    end
    return filter
end

function sliding_demux_dict(seqs,fwd_primers,window::Int,slide_by::Int; verbose = true, phreds = nothing, tol_one_error = true)
    fwd_matches = iterative_primer_match(seqs,fwd_primers,window,slide_by,tol_one_error=tol_one_error)
    rev_comp_bool = fwd_matches .< 0
    keepers = abs.(fwd_matches) .> 0
    fwd_matches = abs.(fwd_matches)
    pair_keeps = fwd_matches[keepers]
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            println(s[1], " => ", s[2])
        end
    end
    if phreds == nothing
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]

                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            end
        end
        return seq_dict
    else
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Vector{Int8},Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            end
        end
        return seq_dict
    end
end

#define nextera indexes
N7_dic = Dict(
    "N701" => "TCGCCTTA",
    "N702" => "CTAGTACG",
    "N703" => "TTCTGCCT",
    "N704" => "GCTCAGGA",
    "N705" => "AGGAGTCC",
    "N706" => "CATGCCTA",
    "N707" => "GTAGAGAG",
    "N708" => "CCTCTCTG",
    "N709" => "AGCGTAGC",
    "N710" => "CAGCCTCG",
    "N711" => "TGCCTCTT",
    "N712" => "TCCTCTAC"
);

S5_dic = Dict(
    "S501" => "TAGATCGC",
    "S502" => "CTCTCTAT",
    "S503" => "TATCCTCT",
    "S504" => "AGAGTAGA",
    "S505" => "GTAAGGAG",
    "S506" => "ACTGCATA",
    "S507" => "AAGGAGTA",
    "S508" => "CTAAGCCT",
    "S517" => "GCGTAAGA"
);

#define universal adapter sequences
N7_univ = "CAAGCAGAAGACGGCATACGAGAT";
S5_univ = "AATGATACGGCGACCACCGAGATCTACAC";

N7_suffix = "GTCTCGTGGGCTCGG"
S5_suffix = "TCGTCGGCAGCGTC"

#"Index" primers
Index_primers_f = Dict(
  "Index_F01" => "CTACACTCGCCTTATCGTCGGCAGCGTC",
  "Index_F02" => "CTACACCTAGTACGTCGTCGGCAGCGTC",
  "Index_F03" => "CTACACTTCTGCCTTCGTCGGCAGCGTC",
  "Index_F04" => "CTACACGCTCAGGATCGTCGGCAGCGTC",
  "Index_F05" => "CTACACAGGAGTCCTCGTCGGCAGCGTC",
  "Index_F06" => "CTACACCATGCCTATCGTCGGCAGCGTC",
  "Index_F07" => "CTACACGTAGAGAGTCGTCGGCAGCGTC",
  "Index_F08" => "CTACACCAGCCTCGTCGTCGGCAGCGTC",
  "Index_F09" => "CTACACTGCCTCTTTCGTCGGCAGCGTC",
  "Index_F10" => "CTACACTCCTCTACTCGTCGGCAGCGTC",
  "Index_F11" => "CTACACTCATGAGCTCGTCGGCAGCGTC",
  "Index_F12" => "CTACACCCTGAGATTCGTCGGCAGCGTC",
  "Index_F13" => "CTACACTAGCGAGTTCGTCGGCAGCGTC",
  "Index_F14" => "CTACACGTAGCTCCTCGTCGGCAGCGTC",
  "Index_F15" => "CTACACTACTACGCTCGTCGGCAGCGTC",
  "Index_F16" => "CTACACAGGCTCCGTCGTCGGCAGCGTC",
  "Index_F17" => "CTACACGCAGCGTATCGTCGGCAGCGTC",
  "Index_F18" => "CTACACCTGCGCATTCGTCGGCAGCGTC",
  "Index_F19" => "CTACACGAGCGCTATCGTCGGCAGCGTC",
  "Index_F20" => "CTACACCGCTCAGTTCGTCGGCAGCGTC",
  "Index_F21" => "CTACACGTCTTAGGTCGTCGGCAGCGTC",
  "Index_F22" => "CTACACACTGATCGTCGTCGGCAGCGTC",
  "Index_F23" => "CTACACTAGCTGCATCGTCGGCAGCGTC",
  "Index_F24" => "CTACACGACGTCGATCGTCGGCAGCGTC"
);

Index_primers_r = Dict(
  "Index_R01" => "CGAGATCTCTCTATGTCTCGTGGGCTCGG",
  "Index_R02" => "CGAGATTATCCTCTGTCTCGTGGGCTCGG",
  "Index_R03" => "CGAGATGTAAGGAGGTCTCGTGGGCTCGG",
  "Index_R04" => "CGAGATACTGCATAGTCTCGTGGGCTCGG",
);

Index_F_univ = "CTACACNNNNNNNNTCGTCGGCAGCGTC"
Index_R_univ = "CGAGATNNNNNNNNGTCTCGTGGGCTCGG"

function find_nextera_suffix(query_seq, query_phred, suffix; start_ix = 34, end_ix = 12, try_reverse_comp = true)
    for ix in start_ix:-1:end_ix
        window = query_seq[ix : ix + length(suffix) - 1]
        d = evaluate(Hamming(), window, suffix)
        if d < 2
            return(query_seq[ix - 8:end], query_phred[ix - 8:end])
        end
    end
    #try reverse complement
    if try_reverse_comp
        query_seq = reverse_complement(query_seq)
        query_phred = query_phred[end:-1:1]
        for ix in start_ix:-1:end_ix
            window = query_seq[ix : ix + length(suffix) - 1]
            d = evaluate(Hamming(), window, suffix)
            if d < 2
                return(query_seq[ix - 8:end], query_phred[ix - 8:end])
            end
        end
    end
end

function get_nextera_matches(seqs, phreds)
    #define universal adapter sequences
    N7_univ = "CAAGCAGAAGACGGCATACGAGAT";
    S5_univ = "AATGATACGGCGACCACCGAGATCTACAC";

    N7_suffix = "GTCTCGTGGGCTCGG";
    S5_suffix = "TCGTCGGCAGCGTC";

    #find N7
    N7_matches = find_nextera_suffix.(seqs, phreds, N7_suffix;
        start_ix = length(N7_univ) + 10, end_ix = 12, try_reverse_comp = true)

    N7_coords = [i for (i,m) in enumerate(N7_matches) if !isnothing(m)]
    N7_keeps = N7_matches[N7_coords]

    #find S5 (no revc)
    matches = find_nextera_suffix.([s for (s,p) in N7_keeps],
        [p for (s,p) in N7_keeps], S5_suffix;
        start_ix = length(S5_univ) + 10, end_ix = 12, try_reverse_comp = true) #find a better solution than this
    coords = [i for (i,m) in enumerate(matches) if !isnothing(m)]
    keeps = matches[coords]

    return [s for (s,p) in keeps], [p for (s,p) in keeps], coords
end

"""
demux CCS based on nextera illumina adapter sequences using a sliding primer match.
Writes collections of reads to .fastq named by index.
"""
function demux_nextera(file; verbose = true)
    if verbose println("Demultiplexing $(file)...") end
    seqs, phreds, seqnames = read_fastq(file);

    #proper usage
    @time matched_seqs, matched_phreds, coords = get_nextera_matches(seqs, phreds);
    names_N7S5 = seqnames[coords];
    nextera_demux_dic = demux_dict(matched_seqs,collect(values(N7_dic)),collect(values(S5_dic));
        phreds = matched_phreds,tol_one_error = false,verbose = false);
    return nextera_demux_dic, names_N7S5
end

t1 = time()
#filter .fastq
#filtered_path = snakemake.output[3] #julia temp paths don't work
@time fastq_filter(input_fastq_path,
                   filtered_path, #path here
                   error_rate = error_rate,
                   min_length = min_length,
                   max_length = max_length)
template_counts = Dict()

println("Running demux with the following parameters...")
println("Index type: $(index_type)")

if snakemake.params["index_type"] == "Nextera_primer"
    nextera_demux_dic,seqnames = demux_nextera(filtered_path)
    nex_tuples = collect(keys(nextera_demux_dic))
    N7 = collect(keys(N7_dic))
    S5 = collect(keys(S5_dic))
    index_tuples = [(N7[x],S5[y]) for (x,y) in nex_tuples]
    indexes2tuples = Dict(zip(index_tuples,nex_tuples))
    for template in collect(keys(SAMPLES))
        indexes = SAMPLES[template]
        if (indexes["fwd_index"],indexes["rev_index"]) in index_tuples
            template_seqs = nextera_demux_dic[indexes2tuples[(indexes["fwd_index"],indexes["rev_index"])]]
            #match template sequences, length here

            if length(template_seqs) < 3 @warn "Less than 3 reads for $(template): $(indexes)" end
            trimmed_seqs = [
                double_primer_trim(s,p,
                N7_suffix*SAMPLES[template]["rev_primer"],S5_suffix*SAMPLES[template]["sec_str_primer"];
                buffer = 8)
            for (s,p) in template_seqs
            ]
            write_fastq(snakemake.output[1]*"/$(template).fastq",
                        [i[1] for i in trimmed_seqs],
                        [i[2] for i in trimmed_seqs];
                        names = seqnames[[i[3] for i in template_seqs]])
        else
            @warn "No reads found for $(template): $(indexes)"
        end
    end
elseif snakemake.params["index_type"] == "Index_primer"
    seqs, phreds, seqnames = read_fastq(filtered_path)
    demux_dic = demux_dict(seqs,
        [i[1:16] for i in collect(values(Index_primers_f))], #setting length to run demux
        [i[1:16] for i in collect(values(Index_primers_r))];
        phreds = phreds,
        tol_one_error = true,
        verbose = false);
    nex_tuples = collect(keys(demux_dic))
    Index_F = collect(keys(Index_primers_f))
    Index_R = collect(keys(Index_primers_r))
    index_tuples = [(Index_F[x],Index_R[y]) for (x,y) in nex_tuples]
    indexes2tuples = Dict(zip(index_tuples,nex_tuples))
    for template in collect(keys(SAMPLES))
    	println("Demultiplexing Sample $(template)")
        indexes = SAMPLES[template]
        if (indexes["fwd_index"],indexes["rev_index"]) in index_tuples
            template_seqs = demux_dic[indexes2tuples[(indexes["fwd_index"],indexes["rev_index"])]]
			index_trimmed = [double_primer_trim(s,p,Index_F_univ,Index_R_univ*SAMPLES[template]["rev_primer"]) for (s,p,n) in template_seqs];
            #match template with forward primer
            keeps = iterative_primer_match([s for (s,p) in index_trimmed], [SAMPLES[template]["sec_str_primer"]],12,25; tol_one_error=true) .> 0 #primer matching, primers needs to be array
            seqs_keeping = index_trimmed[keeps]
            if length(seqs_keeping) == 0
                @warn "No reads found for $(template): $(indexes)";
                template_counts[template] = 0
                continue
            end

            filtered_seqs = [s for (s,p) in seqs_keeping]
            filtered_phreds = [p for (s,p) in seqs_keeping]
            filtered_names = seqnames[[i[3] for i in template_seqs]][keeps]
            template_counts[template] = length(filtered_seqs)
                                  
            if length(filtered_seqs) < 3
                @warn "Less than 3 reads for $(template): $(indexes)"
                template_counts[template] = length(filtered_seqs)
                continue
            end
            
            write_fastq(snakemake.output[1]*"/$(template).fastq",
                        filtered_seqs,
                        filtered_phreds;
                        names = filtered_names)
        else
            @warn "No reads found for $(template): $(indexes)"
            template_counts[template] = 0
        end
    end
else
    @warn "index_type $(snakemake.params["index_type"]) not recognized."
end
t2 = time()
println("Demultiplex took $(t2 - t1) seconds.")
template_df = DataFrame([[(k) for (k,v) in template_counts],[(v) for (k,v) in template_counts]], [:template,:count])
template_df.fwd_index = [SAMPLES[k]["fwd_index"] for (k,v) in template_counts]
template_df.rev_index = [SAMPLES[k]["rev_index"] for (k,v) in template_counts]

template_df.fwd_primer = [SAMPLES[k]["sec_str_primer"] for (k,v) in template_counts]
template_df.rev_primer = [SAMPLES[k]["rev_primer"] for (k,v) in template_counts]
CSV.write(snakemake.output["counts"], template_df)