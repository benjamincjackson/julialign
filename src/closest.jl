using ArgParse

include("functions_closest.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--target", "-t"
            help = "the alignment to search for closest matches, in fasta format"
            required = true
        "--query", "-q"
            help = "the sequence(s) to find matches for, in fasta format"
            required = true
        "--csv-out"
            help = "the name of csv file containing closest matches to write"
            default = "out.csv"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    closest(parsed_args["target"], parsed_args["query"], parsed_args["csv-out"])
end

main()
