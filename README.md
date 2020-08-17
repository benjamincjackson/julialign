# Julialign

A few functions for working with alignments in FASTA format

Requires Julia version `>=` `1.3.1` to be available.

### Commands



| run              | description                                                                        |
|------------------|------------------------------------------------------------------------------------|
| src/bootstrap.jl | Bootstrap an alignment by sampling sites with replacement                          |
| src/collapse.jl  | Heuristic for stripping out the redundancy from a set of similar sequences |
| src/del_typer.jl | Type alignments for pre-specified deletions                                        |


For example, issue

`julia src/collapse.jl -i input.fasta` at the command line to run the collapse function.

Run any function with the `-h` flag to get a full list of options, e.g. `julia src/bootsrap.jl -h`




