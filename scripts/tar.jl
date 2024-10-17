using PorpidPostproc, NextGenSeqUtils

# zip porpid and postproc directories for easy download

dataset = snakemake.wildcards["dataset"]
porpid_dir = "porpid/$(dataset)"
postproc_dir = "postproc/$(dataset)"


# and now tar and zip both porpid and postproc directories
# first copy some reports from porpid to postproc
println("archiving porpid ...")
run(`cp $(porpid_dir)/$(dataset)_contam_report.csv $(postproc_dir)/$(dataset)_contam_report.csv`)

# now rename porpid_dir, zip and rename back
run(`mv $(porpid_dir) $(porpid_dir)-porpid`)
run(`tar -C porpid -czf porpid/$(dataset)-porpid.tar.gz $(dataset)-porpid`)
run(`mv $(porpid_dir)-porpid $(porpid_dir)`)
println("porpid directory archived and zipped ...")

# now do the same with postproc_dir
println("archiving postproc ...")
run(`mv $(postproc_dir) $(postproc_dir)-postproc`)
run(`tar -C postproc -czf postproc/$(dataset)-postproc.tar.gz $(dataset)-postproc`)
run(`mv $(postproc_dir)-postproc $(postproc_dir)`)
println("postproc directory archived and zipped ...")
