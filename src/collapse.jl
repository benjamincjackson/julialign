using ArgParse

include("functions_align.jl")
include("functions_collapse.jl")
include("functions_deletions.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--infile", "-i"
            help = "the alignment to collapse, in fasta format"
            required = true
        "--outfile", "-o"
            help = "the name of the reduced alignment to write"
            default = "out.fasta"
        "--retain", "-r"
            help = "force the retention of this sequence in the output"
            action = "append_arg"
        "--dels-file"
            help = "file of deletions to type and include as variants when collapsing"
        "--append-as-snp"
            action = "store_true"
            help = "if invoked, and a dels file is provided, append the deletion genotypes as SNPs in the output alignment"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    fasta_IDs = get_fasta_descriptions(parsed_args["infile"])
    nucleotide_array = populate_byte_array(parsed_args["infile"])

    if parsed_args["dels-file"] != nothing
        dels = read_del_file(parsed_args["dels-file"])
        nucleotide_array_dels, A_genotypes_table = type_deletions_and_append_as_SNP(nucleotide_array, dels)

        pairwise_identity = is_same(nucleotide_array_dels)
        # for i in 1:size(pairwise_identity, 1)
        #     println(pairwise_identity[i,:])
        # end
        println("sequence comparisons complete")

        final_sets = get_sets_in_one_go(pairwise_identity)
        # println(final_sets)

        levels = TEST_get_levels(final_sets)
        println("number of sequences in input dataset: ", length(fasta_IDs))
        println("number of sequences in output dataset: ", length(levels))

        safe_sets_test = TEST_compare_within_sets(collect(final_sets), nucleotide_array_dels)

        if parsed_args["append-as-snp"]
            tip_to_seq_relationships = write_highest_scoring_unused_seq(final_sets, parsed_args["retain"], nucleotide_array_dels, fasta_IDs, parsed_args["outfile"], false)
        else
            tip_to_seq_relationships = write_highest_scoring_unused_seq(final_sets, parsed_args["retain"], nucleotide_array, fasta_IDs, parsed_args["outfile"], false)
        end

    else
        pairwise_identity = is_same(nucleotide_array)
        # for i in 1:size(pairwise_identity, 1)
        #     println(pairwise_identity[i,:])
        # end
        println("sequence comparisons complete")

        final_sets = get_sets_in_one_go(pairwise_identity)
        # println(final_sets)

        levels = TEST_get_levels(final_sets)
        println("number of sequences in input dataset: ", length(fasta_IDs))
        println("number of sequences in output dataset: ", length(levels))

        safe_sets_test = TEST_compare_within_sets(collect(final_sets), nucleotide_array)

        tip_to_seq_relationships = write_highest_scoring_unused_seq(final_sets, parsed_args["retain"], nucleotide_array, fasta_IDs, parsed_args["outfile"], false)
    end

    println("safe sets test (true/false): ", safe_sets_test)
    println("number of sequences in output alignment: ", length(keys(tip_to_seq_relationships)))

    write_tip_to_redundant_seq_relationships(tip_to_seq_relationships, "tip_to_redundants.csv", parsed_args["retain"])

    seq_to_tip_relationships = get_redundant_seq_to_tip_relationships(tip_to_seq_relationships)
    write_redundant_seq_to_tip_relationships(seq_to_tip_relationships, "redundant_to_tips.csv")
end

main()
