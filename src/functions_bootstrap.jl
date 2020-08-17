using StatsBase

include("functions_align.jl")

function bootstrap(path_to_fasta, n_rep, out_prefix)
    # get the dimensions of the alignment
    array_dim = get_alignment_dimensions(path_to_fasta)
    l = array_dim[1] # number of sites in each sequence
    n = array_dim[2] # number of sequences in the fasta file

    for bootstrap in 1:n_rep
        open(out_prefix * string(bootstrap) * ".fasta", "w") do io

            # V is the sample of sites:
            V = sample(collect(1:l), (l,1), replace = true)

            records = Channel{fasta_record}((channel_arg) -> read_fasta_alignment(path_to_fasta, channel_arg))

            for record in records
                id = record.description
                sequence = record.seq

                # A is an array of one fasta entry's nucleotides
                A = Array{Char,1}(undef, l)
                # BS is the same array, but sampled with replacement (i.e., one bootstrap)
                BS = Array{Char,1}(undef, l)

                counter = 1
                for nuc in sequence
                    A[counter] = nuc
                    counter+=1
                end

                for i in 1:length(V)
                    BS[i] = A[V[i]]
                end

                new_seq = join(BS)

                println(io, ">" * id)
                println(io, new_seq)
            end

        close(io)
        end

    end
end
