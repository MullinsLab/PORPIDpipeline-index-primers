ENV["MPLBACKEND"] = "Agg"
using PorpidPostproc, CSV, NextGenSeqUtils, RobustAmpliconDenoising, BioSequences, DataFrames, DataFramesMeta, FASTX

# include("../../src/functions.jl")
# include("../../src/apobec_model.jl")
# include("../../src/molev_functions.jl")
# include("../../src/postproc_functions.jl")

fasta_collection = snakemake.input[1]*"/"*snakemake.wildcards["sample"]*".fasta"
tag_df = CSV.read(snakemake.input[2], DataFrame)
sample = snakemake.wildcards["sample"]
dataset = snakemake.wildcards["dataset"]
fs_thresh = snakemake.params["fs_thresh"]
agreement_thresh = snakemake.params["agreement_thresh"]
panel_thresh = snakemake.params["panel_thresh"]
degap_flag = snakemake.params["degap"]


#env_seqs = read_fasta("panels/env_column_stripped_panel.fasta")
#env_profile = seqs2profile(uppercase.(env_seqs))

#re_seqs = read_fasta("panels/HIV1_COM_2017_5970-8795_DNA_stripped.fasta")
#re_profile = seqs2profile(uppercase.(re_seqs))

# check for fasta collection
if !isfile(fasta_collection)
    @error "No input FASTA file at $(fasta_collection)! Check if no output sequences from previous step."
    exit()
end

# check for panel file
if !isfile(snakemake.params["panel"])
    @error "Panel $(snakemake.params["panel"]) not found!"
    exit()
end
panel_file = snakemake.params["panel"]


ali_seqs,seqnames = H704_init_template_proc(fasta_collection, panel_file, snakemake.output[1], snakemake.output[2],  snakemake.output[3], snakemake.output[4],  agreement_thresh=agreement_thresh, panel_thresh=panel_thresh)


#copy final seq set to fasta directory
#degap if flag is set to true
println( "degap flag = $(degap_flag)")
if degap_flag == "true"
	write_fasta(snakemake.output[9],
    	degap.(ali_seqs),
    	names = seqnames)
else write_fasta(snakemake.output[9],
    	ali_seqs,
    	names = seqnames)
end    	     	


   
   
#create consensus sequence for each sample from final sequence set
#println("creating consensus from final read set...")  

seq_names,seqs = read_fasta(snakemake.output[9])

#consensus_seqs function can't handle lowercase characters so must be changed to uppercase
seqs = replace.(seqs, r"a" => "A")
seqs = replace.(seqs, r"g" => "G")
seqs = replace.(seqs, r"c" => "C")
seqs = replace.(seqs, r"t" => "T")

draft = consensus_seq(seqs)
draft2 = refine_ref(draft, seqs)
final_cons = [refine_ref(draft2,seqs)] #must be array
cons_name = [sample*"_consensus"] #must be array
out = snakemake.output[10]
write_fasta(out, final_cons, names = cons_name)
      

sp_selected = @linq tag_df |> where(:Sample .== sample)
sp_selected = @linq sp_selected |> where(:tags .!= "BPB-rejects")
fig = family_size_umi_len_stripplot(sp_selected,fs_thresh=fs_thresh)
fig.savefig(snakemake.output[5];
    transparent = true,
    dpi = 200,
    bbox_inches = "tight")
selected = @linq tag_df |> where(:Sample .== sample)
gdf = DataFramesMeta.groupby(selected, :tags)
summary = @combine gdf cols(AsTable) = ( porpid_result=first(:tags), n_UMI_families=length(:fs), n_CCS=sum(:fs) )
# println( summary[!, [:porpid_result,:n_UMI_families,:n_CCS]] )
CSV.write(snakemake.output[6],summary[!, [:porpid_result,:n_UMI_families,:n_CCS]])
