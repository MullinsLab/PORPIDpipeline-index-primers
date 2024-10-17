# import MolecularEvolution
# const ME = MolecularEvolution
using Compose, Colors, StatsBase

############## stolen from MolecularEvolution typed_felsenstein.jl ###################

include("molev_abstracttree.jl")

mutable struct GeneralFelNode <: AbstractTreeNode
  parent::Union{GeneralFelNode,Nothing}
  children::Array{GeneralFelNode,1}
  branchlength::Float64 #Should I just load this into branch_params? In a way, it is special, because it gets stored in the newick file
  name::AbstractString
  nodeindex::Int
  seqindex::Int
  # state_path::Vector{StatePath}
  # branch_params::Vector{Float64} #Parameters that apply to each branch of the tree.
  # message::Vector{Partition} #Used to be called "probs". Represents P(value,all_obs|model), but re-scaled.
  # parent_message::Vector{Partition} #This is meant to be the downward message, but at the top of the branch, before it has been propped over the branch.
  # child_messages::Vector{Vector{Partition}} #One of the above for each child, after prop. For memory reasons, try avoid using this.
  # scaling::Vector{Vector{Float64}}  #A vector of scaling constants for each partition. One scaling constant per site.
  # parent_scaling::Vector{Vector{Float64}}
  #Note: We do not need a child scaling. This is because children are set during the up-pass, and so they get scaled through the message scaling.
  GeneralFelNode(branchlength::Float64, name::AbstractString) = new(nothing,GeneralFelNode[],branchlength,name)
  GeneralFelNode() = GeneralFelNode(0.0, "")
end

broadcastable(x::GeneralFelNode) = Ref(x) #???

function Base.show(io::IO, z::GeneralFelNode)
    println(io, "GeneralFelNode")
    #println(io, "Type: ", typeof(z.branchlength))
    println(io,"nodeindex: $(z.nodeindex)\nRoot: $(isroot(z))\nLeaf: $(isleafnode(z))")
    println(io, "Defined fields:")
    for f in fieldnames(typeof(z))
        println(io, isdefined(z, f), "\t ", f)
    end
end

function print_traversal(node::GeneralFelNode)
    print(node.nodeindex," ", node.branchlength, " ", node.name)
    if isdefined(node,:branch_params)
        print(" ", node.branch_params)
    end
    if !isleafnode(node)
        println(" ",[n.nodeindex for n in node.children])
        for nod in node.children
            print_traversal(nod)
        end
    else
        println()
    end
end

include("molev_tree_draw.jl")

############### end of theft ###################################

function leaf_rename(str)
    #s1 = split(str,"H704_")
    return str #s1[1]*join(split(s1[2],"_")[5:6],"_")
end

function get_size(s)
    return parse(Int64,split(s,"_")[end])
end

function draw_variant_tree(treestring, rename_func, filepath)
    newt = gettreefromnewick(treestring, GeneralFelNode);
    max_dist = maximum(root2tip_distances(newt)[1])
    newt.name = string(round(max_dist,sigdigits = 4))

    for n in getleaflist(newt)
        n.name = rename_func(n.name)
    end

    dot_size_dict = Dict()
    all_sizes = []
    for n in getleaflist(newt)
        dot_size_dict[n] = sqrt(get_size(n.name))
        push!(all_sizes,get_size(n.name))
    end
    tot_size = sum(all_sizes)

    label_color_dict = Dict()
    for n in getleaflist(newt)
        if get_size(n.name) > tot_size/10
            label_color_dict[n] = "red"
        else
            label_color_dict[n] = "black"
        end
    end

    img = tree_draw(newt,dot_color_default = "black",draw_labels = true,label_color_dict = label_color_dict,
        line_width = 0.3, dot_size_dict = dot_size_dict,
        max_dot_size = 0.03, dot_color_dict = label_color_dict,
        dot_opacity = 0.85,
    canvas_height = length(getleaflist(newt))*0.2cm, canvas_width = 10cm
    )
    img |> SVG(filepath*".svg", 10cm, length(getleaflist(newt))*0.2cm)
end

const AA_colors = [
    'Q' => RGB(1.0,0.0,0.8),
    'E' => RGB(1.0,0.0,0.4),
    'D' => RGB(1.0,0.0,0.0),
    'S' => RGB(1.0,0.2,0.0),
    'T' => RGB(1.0,0.4,0.0),
    'G' => RGB(1.0,0.6,0.0),
    'P' => RGB(1.0,0.8,0.0),
    'C' => RGB(1.0,1.0,0.0),
    'A' => RGB(0.8,1.0,0.0),
    'V' => RGB(0.6,1.0,0.0),
    'I' => RGB(0.4,1.0,0.0),
    'L' => RGB(0.2,1.0,0.0),
    'M' => RGB(0.0,1.0,0.0),
    'F' => RGB(0.0,1.0,0.4),
    'Y' => RGB(0.0,1.0,0.8),
    'W' => RGB(0.0,0.8,1.0),
    'H' => RGB(0.0,0.4,1.0),
    'R' => RGB(0.0,0.0,1.0),
    'K' => RGB(0.4,0.0,1.0),
    'N' => RGB(0.8,0.0,1.0),
    '-' => "grey60"
];

const NT_colors = [
      'A' => "salmon",
      'G' => "gold",
      'T' => "lightgreen",
      'C' => "lightblue",
      '-' => "grey60",
]

function highlighter_figure(fasta_collection; out_path = "figure.png")
    #read in and collapse
    seqnames, ali_seqs = read_fasta_with_names(fasta_collection);
    collapsed_seqs, collapsed_sizes, collapsed_names = variant_collapse(ali_seqs; prefix = "v")

	if length(collapsed_seqs) == 1
		println("------------------------------------------------------------------------------")
		println("")
    	@error "All sequences in $(split(basename(fasta_collection),".fast")[1]) are identical, so no collapsed tree can be created"
    	println("")
    	println("------------------------------------------------------------------------------")
    	SVG(out_path, 20cm, 20cm) #blank SVG
	else
   		treestring = fasttree_nuc(collapsed_seqs, collapsed_names; quiet = true)[1];
    	newt = gettreefromnewick(treestring, GeneralFelNode);
    	root = [n for n in getleaflist(newt) if split(n.name,'_')[1] == "v1"][1]
    	rerooted = ladderize(reroot(root))

    	dot_size_dict = Dict()
    	all_sizes = []
    	for n in getleaflist(rerooted)
       		dot_size_dict[n] = sqrt(get_size(n.name))
        	push!(all_sizes,get_size(n.name))
    	end
    	tot_size = sum(all_sizes)

    	color_dict = Dict()
    	for n in getleaflist(rerooted)
        	if get_size(n.name) > tot_size/10
            	color_dict[n] = "red"
        	else
            	color_dict[n] = "black"
        	end
    	end

    	img = highlighter_tree_draw(rerooted, uppercase.(collapsed_seqs), collapsed_names, uppercase(collapsed_seqs[1]);
        	legend_colors = NT_colors, legend_padding = 1cm,
        	tree_args = [:line_width => 0.5mm, :font_size => 1, :dot_size_dict => dot_size_dict, :label_color_dict => color_dict,
        	:dot_color_dict => color_dict, :max_dot_size => 0.1, :dot_opacity => 0.6, :name_opacity => 0.8])
    	compose(context(0.05,0.01,0.95,0.99), img) |> SVG(out_path, 20cm, 20cm)
    end	
end
