# load required packages
using NextGenSeqUtils, StatsBase
using BioSequences, FASTX
using CodecZlib: GzipDecompressorStream

function pp_demux_dict(seqs,fwd_primers,rev_primers; verbose = true, phreds = nothing, tol_one_error = true, demux_dir = "demux")
    if rev_primers == nothing
        fwd_matches = fast_primer_match(seqs,fwd_primers,tol_one_error=tol_one_error)
        rev_comp_bool = fwd_matches .< 0
        keepers = abs.(fwd_matches) .> 0
        fwd_matches = abs.(fwd_matches)
        pair_keeps = fwd_matches[keepers]
    else
         keepers,fwd_matches,rev_matches,rev_comp_bool = fast_primer_pair_match(seqs,fwd_primers,rev_primers,tol_one_error=tol_one_error)
        f_keeps = fwd_matches[keepers]
        r_keeps = rev_matches[keepers]
        pair_keeps = [(f_keeps[i],r_keeps[i]) for i in 1:length(f_keeps)]
    end
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
                if rev_primers == nothing
                    d_key = fwd_matches[i]
                else
                    d_key = (fwd_matches[i],rev_matches[i])
                end
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            else
                write_fasta(demux_dir*"/REJECTS_PRIMER_TRIM.fasta",
                  [seqs[i]],names=["reject-$(i)"],append=true)
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
                if rev_primers == nothing
                    d_key = fwd_matches[i]
                else
                    d_key = (fwd_matches[i],rev_matches[i])
                end
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            else
                write_fastq(demux_dir*"/REJECTS_PRIMER_TRIM.fastq",
                  [seqs[i]],[phreds[i]],names=["reject-$(i)"],append=true)
            end
        end
        return seq_dict
    end
end

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

function iterative_primer_match(seqs,full_primers,
  window::Int,slide_by::Int;tol_one_error=true)
    if(slide_by + window - 1 > minimum(length.(full_primers)))
        @warn("Matching window extends beyond shortest primer. This is ok, but check that you aren't matching something too short.",maxlog=1)
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


function sliding_demux_dict(seqs,fwd_primers,window::Int,slide_by::Int; verbose = true, phreds = nothing, tol_one_error = true, demux_dir = "demux")
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
            else
                write_fasta(demux_dir*"/REJECTS_PRIMER_FWD.fasta",
                  [seqs[i]],names=["reject-$(i)"],append=true)
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
            else
                write_fastq(demux_dir*"/REJECTS_PRIMER_FWD.fastq",
                  [seqs[i]],[phreds[i]],names=["reject-$(i)"],append=true)
            end
        end
        return seq_dict
    end
end


function longest_conserved_5p(seqs)
    for i in 1:length(seqs[1])
        if length(unique(getindex.(seqs,i))) != 1
            return seqs[1][1:i-1]
        end
    end
    return seqs[1]
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
  "Index_R05" => "CGAGATAAGGAGTAGTCTCGTGGGCTCGG",
  "Index_R06" => "CGAGATCTAAGCCTGTCTCGTGGGCTCGG",
  "Index_R07" => "CGAGATCGTCTAATGTCTCGTGGGCTCGG",
  "Index_R08" => "CGAGATTCTCTCCGGTCTCGTGGGCTCGG",
  "Index_R09" => "CGAGATTCGACTAGGTCTCGTGGGCTCGG",
  "Index_R10" => "CGAGATTTCTAGCTGTCTCGTGGGCTCGG",
  "Index_R11" => "CGAGATCCTAGAGTGTCTCGTGGGCTCGG",
  "Index_R12" => "CGAGATGCGTAAGAGTCTCGTGGGCTCGG",
  "Index_R13" => "CGAGATCTATTAAGGTCTCGTGGGCTCGG",
  "Index_R14" => "CGAGATAAGGCTATGTCTCGTGGGCTCGG",
  "Index_R15" => "CGAGATGAGCCTTAGTCTCGTGGGCTCGG",
  "Index_R16" => "CGAGATTTATGCGAGTCTCGTGGGCTCGG",
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


function chunked_filter_apply(fpath, func::Function; chunk_size=10000, f_kwargs = [], verbose = false)
    if endswith(fpath, ".gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(fpath)))
    else
        reader = FASTQ.Reader(open(fpath))
    end
    # if !hasmethod(func, Int64, Int64, Tuple{Array{Any,1}, Array{Phred,1}, Array{Any,1}})
    #    @error "Function argument must accept func(seqs::Array{Any,1}, phreds::Array{Phred,1}, names::Array{Any,1})!"
    # end
    seqs, phreds, names = [], Vector{Phred}[], []
    i = 0
    read_counts = [0, 0, 0]
    chunk = 0
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        push!(seqs, FASTQ.sequence(String, record))
        push!(phreds, FASTQ.quality(record, :sanger))
        push!(names, FASTQ.identifier(record))
        i += 1
        if i == chunk_size
            #apply func...
            chunk += 1
            println("processing chunk $(chunk), of size $(chunk_size)")
            read_counts += func(chunk, chunk_size, seqs, phreds, names; f_kwargs...)
            seqs, phreds, names = [], Vector{Phred}[], []
            i = 0
        end
    end
    close(reader)
    if i > 0
        #apply func...
        chunk += 1
        println("processing chunk $(chunk), of size $(i)")
        read_counts += func(chunk, i, seqs, phreds, names; f_kwargs...)
    end
    return read_counts   # [total_reads, quality_reads, demuxed_reads]
end

#------Chunked quality demux function--------
"""
function chunked_quality_demux(chunk, chunk_size, seqs, phreds, names;
    demux_dir = "demux", samples = Dic(), error_rate = 0.01, min_length = 30,
    max_length = 1000000, label_prefix = "seq", error_out = true, verbose = false)

Intended to be passed to `chunked_filter_apply()` for quality filtering and demultiplexing of
FASTQ files in one pass.

The code to demultiplex samples using indexes (lines 410-504) was originally written by 
Alec Pankow for use with an older version of porpidpostproc. It has been adapted and added
here by Dylan Westfall to replace the demultiplexing code that was originally part of this
script that was created by Hugh Murrell as part of the h705mod1 branch of porpidpostproc. 

"""
function chunked_quality_demux(chunk, chunk_size, seqs, phreds, names;
    demux_dir = "demux",
    samples = Dic(),
    error_rate = 0.01,
    index_type = "Index_primer",
    min_length = 30,
    max_length = 1000000,
    label_prefix = "seq",
    error_out = true,
    verbose = false)
    
    t1 = time()
    
    total_reads, quality_reads, demuxed_reads = 0, 0, 0
    
    total_reads = length(seqs)
    println("total reads =$(total_reads)")  

    # quality filter
    if verbose
        println("filtering chunk on mean phred scores ...")
    end
    lengths = length.(seqs)
    mean_errors = [mean(phred_to_p.(phred)) for phred in phreds]
    inds = [1:length(seqs);][(lengths .< max_length) .& (lengths .> min_length) .& (mean_errors .< error_rate)]
    
    # save the rejected sequences as a fastq file
    rejects = setdiff(1:length(seqs),inds)
    write_fastq(demux_dir*"/REJECTS_DEMUX_QUAL.fastq",
        seqs[rejects],phreds[rejects],names=names[rejects],append=true)
    
    # process the passed sequences
    if error_out == true
        names = ["$label_prefix$((chunk-1)*chunk_size+i)|ee=$(mean_errors[i])" for i in inds]
    else
        names = ["$label_prefix$((chunk-1)*chunk_size+i)" for i in 1:length(inds)]
    end
    seqs, phreds = seqs[inds], phreds[inds]
    
    quality_reads = length(seqs)
    println("quality_reads =$(quality_reads)")  

	template_counts = Dict()

	if index_type == "Nextera_primer"
    	nextera_demux_dic,seqnames = demux_nextera(seqs)
    	nex_tuples = collect(keys(nextera_demux_dic))
    	N7 = collect(keys(N7_dic))
    	S5 = collect(keys(S5_dic))
    	index_tuples = [(N7[x],S5[y]) for (x,y) in nex_tuples]
    	indexes2tuples = Dict(zip(index_tuples,nex_tuples))
    	for template in collect(keys(samples))
        	indexes = samples[template]
        	if (indexes["fwd_index"],indexes["rev_index"]) in index_tuples
            	template_seqs = nextera_demux_dic[indexes2tuples[(indexes["fwd_index"],indexes["rev_index"])]]
            	#match template sequences, length here

            	if length(template_seqs) < 3 @warn "Less than 3 reads for $(template): $(indexes)" end
            	trimmed_seqs = [
                	double_primer_trim(s,p,
                	N7_suffix*samples[template]["rev_primer"],S5_suffix*samples[template]["sec_str_primer"];
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
	elseif index_type == "Index_primer"
    	#seqs, phreds, seqnames = read_fastq(seqs)
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
    	for template in collect(keys(samples))
    		#println("Demultiplexing Sample $(template)")
        	indexes = samples[template]
        	if (indexes["fwd_index"],indexes["rev_index"]) in index_tuples
           		template_seqs = demux_dic[indexes2tuples[(indexes["fwd_index"],indexes["rev_index"])]]
				index_trimmed = [double_primer_trim(s,p,Index_F_univ,Index_R_univ*samples[template]["rev_primer"]) for (s,p,n) in template_seqs];
            	#match template with forward primer
            	keeps = iterative_primer_match([s for (s,p) in index_trimmed], [samples[template]["sec_str_primer"]],12,25; tol_one_error=true) .> 0 #primer matching, primers needs to be array
            	seqs_keeping = index_trimmed[keeps]
            	final_trim = [primer_trim(s,p,samples[template]["sec_str_primer"]) for (s,p) in seqs_keeping]; #now trim off forward primer
            	
            	if length(final_trim) == 0
                	@warn "No reads found for $(template)";
                	template_counts[template] = 0
                	continue
            	end

            	filtered_seqs = [s for (s,p) in final_trim]
            	filtered_phreds = [p for (s,p) in final_trim]
            	filtered_names = names[[i[3] for i in template_seqs]][keeps]
            	template_counts[template] = length(filtered_seqs)
                                  
            	if length(filtered_seqs) < 3
                	@warn "Less than 3 reads for $(template)"
                	template_counts[template] = length(filtered_seqs)
                	continue
            	end
                
                #because reading and writing the fastq files takes time
                #recommend using larger chunk sizes like 100000  
            	if isfile(demux_dir*"/$(template).fastq")
    				seqs2,phreds2,names2 = read_fastq(demux_dir*"/$(template).fastq")
    				cat_seqs = vcat(filtered_seqs, seqs2)
    				cat_phreds = vcat(filtered_phreds, phreds2)
    				cat_names = vcat(filtered_names, names2)
    				write_fastq(demux_dir*"/$(template).fastq",
                        	cat_seqs,
                        	cat_phreds;
                        	names = cat_names)
				else 
					write_fastq(demux_dir*"/$(template).fastq",
                        	filtered_seqs,
                        	filtered_phreds;
                        	names = filtered_names)
				end
            	
        	else
            	@warn "No reads found for $(template)"
            	template_counts[template] = 0
        	end
    	end
	else
    	@warn "index_type $(snakemake.params["index_type"]) not recognized."
	end
	    
    t2 = time()
	println("Demultiplexing chunk took $(t2 - t1) seconds.")  
    
    return [total_reads, quality_reads, demuxed_reads]

end
