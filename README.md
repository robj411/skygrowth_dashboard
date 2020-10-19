# Skygrowth dashboard

## 1. Prerequisites

### System requirements

To run in an environment:

```
module load anaconda3/personal
anaconda-setup
conda create -n snakemake
source activate snakemake
conda install -c bioconda snakemake
conda install iqtree
conda install R
```

To use a profile for cluster parallelisation:

```
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake
cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git profile_name=profile
```

To the job script, one can add module loads and paths as required:

```
module load anaconda3/personal
```

### R requirements

```
install.packages('ape')
install.packages('phangorn')
install.packages('lubridate')
install.packages('limSolve')
install.packages('treedater')
install.packages('dplyr')
install.packages('DT')
install.packages('stringr')
install.packages('shiny')
install.packages('shinyjs')
install.packages('grid')
install.packages('scales')
install.packages('ggplot2')
install.packages('plotly')
install.packages('RColorBrewer')
install.packages('plotrix')
install.packages('svDialogs')
install.packages('devtools')
require('devtools')
install_github('emvolz-phylodynamics/sarscov2Rutils')
install_github('YuLab-SMU/ggtree')
install_github('mrc-ide/skygrowth')
```

### Other requirements

An alignment is required, such as can be downloaded and extracted from [www.gisaid.org](https://www.gisaid.org), e.g. following [JorgensenD/sarscov2_phylo_pipeline](https://github.com/JorgensenD/sarscov2_phylo_pipeline), steps 1 to 3. It should be saved as `algn3.fasta` in the target directory, e.g. [Peru/algn3.fasta](Peru/algn3.fasta). Sequence names should include a GISAID ID and a date in the following format (pipe separated, in positions 2 and 3 respectively): `<name> | <GISDAID ID> | <date> | ...`, e.g. `hCoV-19/Germany/BAV-V2010492/2020|EPI_ISL_420899|2020-03-11|2020.1912568306|_Il`. Sequences whose names include the string "exog" will be identified as exogenous.

## 2. Running the pipeline

```
module load anaconda3/personal
source activate snakemake
snakemake --config region="Norway" --profile profile
```

or, to run locally using one core,

```
snakemake --config region="Norway" -j1
```

## 3. Running the app

In R, run the app locally, or remotely:

```
library(shiny)
runGitHub('robj411/skygrowth_dashboard')
```
