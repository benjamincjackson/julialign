using ArgParse
include("functions_bootstrap.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile", "-i"
            help = "the alignment to bootstrap, in fasta format"
            required = true
        "--out-prefix", "-p"
            help = "prefix for output filenames"
            default = "bootstrap_"
        "--n-bootstraps", "-n"
            help = "how many bootstraps to carry out"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    bootstrap(parsed_args["infile"], parsed_args["n-bootstraps"], parsed_args["out-prefix"])
end

main()
