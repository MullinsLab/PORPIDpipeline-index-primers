# PorpidPostproc using Indexed Reads

This pipeline has been adapted from the PORPID-postproc repository (https://github.com/MurrellGroup/PORPIDpipeline/)
by Dylan Westfall to enable the pipeline to accept samples labeled with Index primers rather 
than unique donor ID barcodes present in the cDNA primer. This primarily revolved around 
replacing the demultiplexing code in demux_functions.jl with code which Alec Pankow had 
written previously for an older version of PORPID-postproc that could accept Indexed samples. 

The new code block does not report counts for reads that were rejected so the following lines
are now not included in the {dataset}-count.html report. 
REJECTS_PRIMER_REV
REJECTS_PRIMER_TRIM

To avoid confusion the index rule, script index.jl, and report {dataset}-index.html in the
h705mod1 branch were renamed to counts, counts.jl, and {dataset}-counts.html to avoid confusion
now that the demultiplexing step uses Index primers.

## Index Version 1.1 changes
1) All csv reports (contam, demux, quality, etc) now include dataset name
2) Added a directory in postproc **{dataset}_fasta**, which contains the final postproc sequence 
   set for each sample. If degap is TRUE then these files are degapped. The fasta file in 
   the sample's postproc folder remains an alignment of the final sequences regardless of degap flag.
   This moved the degap flag parameter code for creation of the final sequence file from tar to postproc.
3) Error 1 - Highlighter/tree plots could not be created in the report.html for samples with 
   fewer than 3 final sequences and caused the pipeline to abort. Now the report.html file 
   will contain a blank highlighter plot and tree image for samples with fewer than 3 final sequences.
4) Error 2 - Contamination report table in counts.html could not be created for samples with only 1
   discarded sequence. Code has been edited to fix this issue and table displays correctly
   for any number of discarded sequences.
5) Error 3 - Samples with a "." in the name were inadvertently truncated during contamination
   check step by splitting fasta file names using '.' This has been resolved by splitting using
   ".fast"


## Quick start

### Dependencies (on an ubuntu machine)

- first update all apps
   - `apt update`
   - `apt upgrade`
- Snakemake
   - `apt-get install -y snakemake`
- mafft
   - `apt-get install -y mafft`
- fasttree
   - `apt-get install -y fasttree`
- python3 packages
  - `apt-get install python3-pandas`
  - `apt-get install python3-seaborn`


### Julia version 1.7

Download and unpack the latest Julia (we recommend version 1.7.1) from: 

[https://julialang.org/downloads/](https://julialang.org/downloads/)

Make sure you can enter the julia REPL from the command line, on an ubuntu machine you would do:

```bash
# move the julia system to a lib directory
mv julia-1.7.1 /usr/lib/julia-1.7.1
# make julia v1.7.1 executable from the command line
ln -s /usr/lib/julia-1.7.1/bin/julia /usr/local/bin/julia
# check that you can enter the julia REPL
julia --version
```

### cloning the PorpidPostproc repository

Now that the dependencies are setup we clone the PorpidPostproc repository

```bash
cd ~
git clone -b h705mod1 https://gitlab.com/hugh.murrell/porpidpostproc.git
```

### setting up the Julia package environment

then navigate to the `porpidpostproc` project folder and start the Julia REPL. 
Enter the package manager using `]` and then enter

```julia
activate .
instantiate
precompile
```

This will activate, install, and precompile the `julia` environment specified by the 
`Project.toml` and `Manifest.toml` files. The `precompile` command
above is not strictly needed but is useful if there are issues with installing
the `julia` packages listed in `Project.toml`

Next, add the following text to your Julia startup file (typically at `~/.julia/config/startup.jl`; 
you may need to create the directory if not present, `mkdir -p ~/.julia/config`).

```julia
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
```

This will activate the local environment at Julia startup.


### Configuration

To configure the PorpidPostproc workflow, first edit the demo `config.yaml` file to reflect
your library construction. 
It should follow the same format as shown below in **example.yaml**

```yaml
Dataset1:
  Sample1:
    cDNA_primer: CCGCTCCGTCCGACGACTCACTATAcactcaNNNNNNNNCAATCAKCACCTGCCATCTGTTTTCCAT
    fwd_index: Index_F06
    index_type: Index_primer
    panel: panels/POL-FmodMozFl-INTR2_panel.fasta
    rev_index: Index_R01
    rev_primer: CCGCTCCGTCCGACGACTCACTATA
    sec_str_primer: CTGTACCAGTAAAATTAAAGCCAGGAATGGATGG
  Sample2:
    cDNA_primer: CCGCTCCGTCCGACGACTCACTATAcactcaNNNNNNNNCAATCAKCACCTGCCATCTGTTTTCCAT
    fwd_index: Index_F04
    index_type: Index_primer
    panel: panels/POL-FmodMozFl-INTR2_panel.fasta
    rev_index: Index_R01
    rev_primer: CCGCTCCGTCCGACGACTCACTATA
    sec_str_primer: CTGTACCAGTAAAATTAAAGCCAGGAATGGATGG
  Sample3:
    cDNA_primer: CCGCTCCGTCCGACGACTCACTATAcactcaNNNNNNNNCAATCAKCACCTGCCATCTGTTTTCCAT
    fwd_index: Index_F05
    index_type: Index_primer
    panel: panels/POL-FmodMozFl-INTR2_panel.fasta
    rev_index: Index_R01
    rev_primer: CCGCTCCGTCCGACGACTCACTATA
    sec_str_primer: CTGTACCAGTAAAATTAAAGCCAGGAATGGATGG
```

The primer sequences provided will be used for demultiplexing and will be trimmed
from the final sequences. **sec_str_primer** and **rev_primer** are the PCR primers 
used in the 2nd rd of PCR. **fwd_index** and **rev_index** are the index primers used to 
tag each sample. Illumina Nextera indexes can be used as well as the Index primer set 
developed in Jim Mullins lab. Information on the Index primer set can be found in the file
Mullins_Index_Primers.csv in the docs folder. 

Note that the **donor ID barcode** is in lowercase
and the **Unique Molecular Identifier (UMI) barcode** is indicated with N's.

The **panel** arg should be a path to a `.fasta` alignment spanning your amplicon, 
with all gaps stripped. This will be used only in the postproccessing step to remove 
off-target seqs and trim to the correct coordinates.

Raw CCS .fastq files should be placed in the `raw-reads/` subdirectory and named 
according to the the dataset name used in the `config.yaml` file, ie, `Dataset1.fastq`
for the example dataset.

An introduction to PacBio sequencing and an explanation for each *PorpidPostproc* rule 
is given in the set of introductory slides packaged with this repository, although the info 
for the demux rule has not been updated to reflect usage with Index tags.
[docs/slides/PorpidPostprocPipeline.pdf](docs/slides/PorpidPostprocPipeline.pdf)


### Preview and Execution

Preview jobs with Snakemake and run with {n} cores.

```bash
#preview jobs
snakemake -np

#run
snakemake -j{n}
```

For more info on Snakemake, see https://snakemake.readthedocs.io/en/stable/

## Conda setup

Some (without root access) may prefer to setup PorpidPostproc in a **conda** environment.

To accomplish this, first install `anaconda` locally. (the install script allows you to choose
the location for anaconda, by default `/home/user` but choose something else if
you want something accessable to a group of users)

```bash
curl –O https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh > Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
```

then log out and log in again and check that you are in the `base` environment.

`conda` is very slow, so we suggest installing `mamba` in the conda `base` environment:

```bash
conda install -n base -c conda-forge mamba
```
clone the PorpidPostproc repository

```bash
cd ~  # or some other directory used for your anaconda installation
git clone -b h705mod1 https://gitlab.com/hugh.murrell/porpidpostproc.git
```

and then all the PorpidPostproc dependencies including `julia` version `1.7.1`
( as listed in the PorpidPostproc conda environment spec in `environment.yaml`),
can be installed in a `conda` environment via `mamba` using the commands:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
mamba env create --file environment.yaml
```

Note that if you did use *some other directory* than your home directory for
installing the PorpidPostproc repository then you have to inform Julia where
your packages are stored by placing the following command in your `.bashrc`
file:

```bash
# set path to .julia files
export JULIA_DEPOT_PATH="/some/other/directory/.julia"
```

to complete the setup, activate the new PorpidPostproc conda environment, 

```bash
conda activate PorpidPostproc
```

and continue with the `julia` package environment setup as outlined above in the *quick start* section.

## Cluster setup

Seting up a `snakemake` pipeline on a cluster is a *dark* art. Here we describe an attempt
at installing PorpidPostproc on a two node cluster, (one node a *controller* node with 16 cores
and the other node a *compute* node with 64 cores).

**Firstly**, since the cluster administrator is hardly likely to give you root access we
suggest you follow the `conda` installation for PorpidPostproc. If you expect more
than one user of your PorpidPostproc pipeline then install in a directory that all
your users can see and that is visible from both the *contoller* and *compute* nodes. 
ie use `some other directory` rather than the standard home directory and make 
sure to inform `julia` about this choice of directory as
outlined in the `conda` section above.

**Secondly**, cluster administrators usually insist that large data sets are stored
in an appropriate volume and **not** in the usual user's space. On our cluster the
administrator required the PorpidPostproc code to be installed in a `\tools\porpid\`
directory and the large data sets (input, output and temporary) to be stored in
a `\data\porpid\` directory so we installed PorpidPostproc into `\tools\porpid\porpidpostproc`
and then replaced some of the directories in the `porpidpostproc`
directory with symbolic links to an appropriate directory in the `\data\porpid\` directory
as shown below

```
config.yaml -> /raw/porpid/config/demo.yaml
panels -> /raw/porpid/panels/
porpid -> /raw/porpid/porpid/
postproc -> /raw/porpid/postproc/
raw-reads -> /raw/porpid/raw-reads/
```

Naturally, one must copy contents of the installation to the `/raw/porpid/` directory
before deleting the installation directory and replacing it with a symbolic link to the
appropriate place on the `raw` volume.

**Job submission**, after setting up like this we are ready to run the `demo` study through PorpidPostproc
by submitting the `snakemake` command to the cluster managemant system.
On our cluster that management system is `slurm` and the following shell script
stored in `porpid_job.sh` facilitated that submission:

```bash
#!/bin/bash
#SBATCH --job-name==porpid
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --partition=main

if [ "$#" -lt 1 ]; then
    echo "please supply a config file name as first parameter"
    exit
fi
echo "config file is $1"

echo "${SLURM_JOB_NAME} job submited using ${SLURM_NTASKS} cores"

# create a symbolic link for the snakemake config file to point to the config for the current study
rm -f /tools/PorpidPostproc/porpidpostproc/config.yaml
ln -s /RAW/PORPID/CONFIG/$1.yaml /tools/PorpidPostproc/porpidpostproc/config.yaml

# tell slurm where anaconda is and conda activate the PorpidPostproc environment
source /tools/PorpidPostproc/anaconda3/etc/profile.d/conda.sh
conda activate PorpidPostproc

# navigate to the porpidpostproc directory and run snakemake
# add -F to to the snakemake command to force re-run of all rules
cd /tools/PorpidPostproc/porpidpostproc
snakemake --rerun-incomplete -j${SLURM_NTASKS}  
```

To submit the `demo` to run as a `slurm` batch job one just uses

```bash
sbatch porpid_job.sh demo
```
The script above sets some environment variables for `slurm` and then resets
the symbolic link to the appropriate config file for the `demo` study.
It then activates the conda environment switches to the installation
directory and runs the snakemake pipeline.

With this structure it is easy to run a new study through PorpidPostproc.
One copies the new config file into the `/raw/porpid/config/` directory,
transfers the `fastq` data to the `/raw/porpid/raw-reads/` directory
and then issues the `sbatch` command using the appropriate study name
instead of `demo`

Note that with this method you must predetermine the number of cores
you intend to use on your cluster's node. In the `demo` study this is set
to 7 ( 6 cores for the samples to run in parallel plus 1 core for snakemake )

Each study will be different. To see how many samples can be run in parallel
you can do a `snakemake` dry run using the `porpid_dry_run.sh` script below:

```bash
#!/bin/bash
if [ "$#" -lt 1 ]; then
        echo "please supply a config file name as first parameter"
        exit
fi
echo "config file is $1"
# create a symbolic link for the snakemake config file to
# point to the config for the current study
rm -f /tools/PorpidPostproc/porpidpostproc/config.yaml
ln -s /RAW/PORPID/CONFIG/$1.yaml /tools/PorpidPostproc/porpidpostproc/config.yaml
# activate the conda environment
source /tools/PorpidPostproc/anaconda3/etc/profile.d/conda.sh
conda activate PorpidPostproc
# perform a snakemake dry run
# remove the -f for a partial dry run of what's left to do
cd /tools/PorpidPostproc/porpidpostproc
snakemake -F --rerun-incomplete -np
```

Note that this dry run is not compute intensive and can ve executed on the
*controller* machine without using the `sbatch` command as follows:

```bash
./porpid_dry_run.sh demo
```

### Caveat

The above suggestion for running a `snakemake` pipeline under `slurm`
is rudamentary. Maximum cores must be requested at the start of execution
and they are probably held throughout the run.

However, it is alledged that `snakemake` can play nicely with `slurm` and
it should be possible to have `snakemake` invoke `slurm` for each rule in
the pipeline. In this case `snakemake` would request the optimal number
of cores needed for each step in the pipeline.

We have not attempted this yet, and it would probably require writing a
`slurm` efficient version of the `snakefile`. 

Watch this space for further developments.

## Documentation

### Workflow

The graph below summarizes the overall organization of the workflow. 
Each node in the graph is a *rule* in the The [Snakefile](Snakefile).

![rulegraph](rulegraph.png)





