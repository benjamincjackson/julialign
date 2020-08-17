using ArgParse

include("functions_align.jl")
include("functions_deletions.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile", "-i"
            help = "the alignment to type, in fasta format"
            required = true
        "--csv-out"
            help = "the name of csv file to write with deletion genotypes"
            default = "out.csv"
        "--dels-file"
            help = "file of deletions to type and include as variants when collapsing"
            required = true
        "--fasta-out"
            help = "the name of fasta file to write if --append-as-snp is invoked"
            default = "out.fasta"
        "--append-as-snp"
            action = "store_true"
            help = "if invoked, append the deletion genotypes as SNPs in the output alignment"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    dels = read_del_file(parsed_args["dels-file"])

    fasta_IDs = get_fasta_descriptions(parsed_args["infile"])
    nucleotide_array = populate_byte_array(parsed_args["infile"])

    if parsed_args["append-as-snp"]
        nucleotide_array_dels, genotypes_table = type_deletions_and_append_as_SNP(nucleotide_array, dels)

        open(parsed_args["csv-out"], "w") do io
            println(io, "sequence," * join([join(x, "_") for x in dels], ","))
            for i in 1:size(genotypes_table, 2)
                println(io, fasta_IDs[i] * "," * join(genotypes_table[:,i], ","))
            end
        close(io)
        end

        open(parsed_args["fasta-out"], "w") do io
            for i in 1:size(nucleotide_array_dels, 2)
                println(io, ">" * fasta_IDs[i])
                println(io, get_seq_from_1D_byte_array(nucleotide_array_dels[:,i]))
            end
        close(io)
        end

    else
        genotypes_table = type_deletions(nucleotide_array, dels)

        open(parsed_args["csv-out"], "w") do io
            println(io, "sequence," * join([join(x, "_") for x in dels], ","))
            for i in 1:size(genotypes_table, 2)
                println(io, fasta_IDs[i] * "," * join(genotypes_table[:,i], ","))
            end
        close(io)
        end
    end
end

#=
NB this is slow at the moment because everything goes into an array before
typing it (side effect of writing collapse first)
=#

main()
