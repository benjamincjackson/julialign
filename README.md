# Julialign

A few functions for working with alignments in FASTA format

Requires Julia version `>=` `1.3.1`

Julialign uses a slightly modified version of the bit-level coding scheme for nucleotides developed by Emmanuel Paradis (described [here](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html), and implemented in the R package [ape](https://doi.org/10.1093/bioinformatics/btg412)).

### Commands



| run              | description                                                                        |
|------------------|------------------------------------------------------------------------------------|
| src/bootstrap.jl | Bootstrap an alignment by sampling sites with replacement                          |
| src/closest.jl   | Find the closest sequence to a query by SNP-distance                             |
| src/collapse.jl  | Heuristic for stripping out the redundancy from a set of similar sequences |
| src/del_typer.jl | Type alignments for pre-specified deletions                                        |


For example, issue

`julia src/collapse.jl -i input.fasta` at the command line to run the collapse function.

Run any function with the `-h` flag to get a full list of options, e.g. `julia src/bootsrap.jl -h`




