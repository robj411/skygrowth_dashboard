# Skygrowth dashboard

## Pre-requisites

### System requirements

To run in an environment:

```
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

An alignment is required, such as can be downloaded and extracted from [www.gisaid.org](www.gisaid.org). It should be saved as `algn3.fasta` in the target directory, e.g. [Peru/algn3.fasta](Peru/algn3.fasta). 

## Running the pipeline

```
module load anaconda3/personal
source activate snakemake
snakemake --config population=5400000 region="Norway" --profile profile
```

or, to run locally using one core,

```
snakemake --config population=5400000 region="Norway" -j1
```

## Running the app

In R, run the app locally, or remotely:

```
library(shiny)
runGitHub('robj411/skygrowth_dashboard')
```
