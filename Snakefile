wildcard_constraints:
  region='[^/]+'

configfile: "config.yaml"

rule all:
  input:
    "input_files.tar",
    expand("input_files/{rg}.Rds",rg=config['region'])

rule zip:
  input:
    expand("input_files/{rg}.Rds",rg=config['region'])
  output:
    "input_files.tar"
  shell:
    "tar cvf {output} input_files" 

rule process_for_shiny:
  input:
    "{region}/trees.Rds",
    "{region}/skygrowth3/skygrowth3-sgs.rds"
  output:
    "input_files/{region}.Rds"
  params:
    rg="{region}"
  shell:
    "R --vanilla --args '{params.rg}' < process_for_shiny.R > '{params.rg}.log' 2> '{params.rg}.log'"

rule run_skygrowth:
  input:
    "{region}/trees.Rds"
  output:
    "{region}/skygrowth3/skygrowth3-sgs.rds"
  params:
    rg="{region}"
  shell:
    "R --vanilla --args '{params.rg}' < run_skygrowth.R  > '{params.rg}.log' 2> '{params.rg}.log'"


rule make_trees:
  input:
    "{region}/algn3.fasta"
  output:
    "{region}/trees.Rds",
    "{region}/startTrees.nwk"
  params:
    rg="{region}"
  shell:
    "R --vanilla --args '{params.rg}' < make_trees.R  > '{params.rg}.log' 2> '{params.rg}.log'"

rule make_alignment:
  output:
    "{region}/algn3.fasta"
  params:
    rg="{region}",
    nrg=config['n_region'],
    nres=config['n_reservoir'],
    miss=config['missingness_threshold']
  shell:
    "R --vanilla --args '{params.rg}' '{params.nrg}' '{params.nres}' '{params.miss}' < make_alignment.R  > '{params.rg}.log' 2> '{params.rg}.log'"




