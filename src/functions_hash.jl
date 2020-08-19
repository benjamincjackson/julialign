include("functions_align.jl")

function seq_hash(target_file, query_file, outfile)

    target_dim = get_alignment_dimensions(target_file)
    query_dim = get_alignment_dimensions(query_file)

    target_hash_Dict = Dict{UInt64, Array{String, 1}}()

    target_channel = Channel{fasta_record}((channel_arg) -> read_fasta_alignment(target_file, channel_arg))

    for record in target_channel
        h = hash(record.seq, hash(record.seq))
        if in(h, keys(target_hash_Dict))
            push!(target_hash_Dict, h => vcat(target_hash_Dict[h], record.description))
        else
            push!(target_hash_Dict, h => [record.description])
        end
    end

    query_channel = Channel{fasta_record}((channel_arg) -> read_fasta_alignment(query_file, channel_arg))

    for record in query_channel
        h = hash(record.seq, hash(record.seq))
        if in(h, keys(target_hash_Dict))
            println(record.description * " is exactly the same as: " * join(target_hash_Dict[h], ", "))
        end
    end

end
















#
