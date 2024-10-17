ENV["MPLBACKEND"] = "Agg"
using PorpidPostproc, NextGenSeqUtils, StatsBase
using BioSequences, FASTX
using DataFrames, CSV
# using CodecZlib: GzipDecompressorStream
# include("../../src/postproc_functions.jl")

function NextGenSeqUtils.write_fastq(filename, seqs, phreds::Vector{Vector{Phred}};
                     names=String[], LongSequence = false,
                     append = false)
    if !LongSequence
        seqs = [LongCharSeq(s) for s in seqs]
    end
    stream = open(FASTQ.Writer, filename, append=append)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, q, n) in zip(seqs, phreds, names)
        i += 1
        write(stream, FASTQ.Record(n, s, q))
    end
    close(stream)
end

println("using Julia version: $(VERSION)")

t1 = time()

SAMPLE_CONFIGS = snakemake.params["config"]
index_type = snakemake.params["index_type"]
mkdir(snakemake.output[1])


f_kwargs = [
    :demux_dir => snakemake.output[1],
    :samples => SAMPLE_CONFIGS,
    :verbose => false,
    :index_type => snakemake.params["index_type"],
    :error_rate => snakemake.params["error_rate"],
    :min_length => snakemake.params["min_length"],
    :max_length => snakemake.params["max_length"],
    :label_prefix => "seq",
    :error_out => true
    ]

println("performing chunked quality filtering and demux on $(snakemake.input[1])")
chunk_size = snakemake.params["chunk_size"]
println("Index type: $(index_type)")

reads = chunked_filter_apply(snakemake.input[1], chunked_quality_demux; chunk_size=chunk_size, f_kwargs, verbose = false)

# create empty fastq files for those samples with no reads
# this does not work, too much trouble with downstream scripts
# for sample in snakemake.params["SAMPLES"]
#     println(sample)
#     run(`touch "$(snakemake.output[1])/$(sample).fastq"`)
# end

# now do some accounting
total_reads = reads[1]
quality_reads = reads[2]
demuxed_reads = reads[3]

# report on total sequences for each sample
println()
println("total reads => $(total_reads)")
println("quality reads => $(quality_reads)")
println("-------------------------------------")
println("demuxed reads:")
no_assigned = 0
no_rejected = 0
filepaths = readdir(snakemake.output[1],join=true)
df_demux = DataFrame(Sample = [], Read_Count = [])
df_reject = DataFrame(Sample = [], Read_Count = [])
#println("demuxed reads => $(demuxed_reads)")
for path in filepaths
    stream = open(FASTQ.Reader, path)
    count = 0
    for record in stream
        count += 1
        if occursin("REJECT",path)
            global no_rejected += 1
        else
            global no_assigned += 1
        end
    end
    sample_name = replace(basename(path), ".fastq" => "")
    println(sample_name," => ",count)
    if occursin("REJECT",path)
        push!(df_reject,[sample_name,count])
    else
        push!(df_demux,[sample_name,count])
    end
end
println("-------------------------------------")
println("total demuxed => $(no_assigned)")
println("total rejected => $(no_rejected)")
println("")
println("If read counts are lower than expected check that min and mix read lengths match amplicon size")
println("")
CSV.write("$(snakemake.output[3])", df_demux)
CSV.write("$(snakemake.output[4])", df_reject)

# save a quality_filter report
# no_of_fails = no_of_reads - no_assigned # fix this, get chunked_filter_apply to return no_of_fails
df_qual = DataFrame(Description = [], Count = [])
push!(df_qual,["Initial number of raw reads",total_reads])
push!(df_qual,["Number of quality reads",quality_reads])
push!(df_qual,["Number assigned by demux",no_assigned])
CSV.write("$(snakemake.output[2])", df_qual)


t2 = time()
println("Quality filtering and demultiplexing took $(t2 - t1) seconds.")
