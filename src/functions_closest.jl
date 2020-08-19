using StatsBase

include("functions_align.jl")


function get_difference_matrix(target_A, query_A)

    difference_array = Array{Float64,2}(undef, size(target_A, 2), size(query_A, 2))

    for query_index in 1:size(query_A, 2)
        for target_index in 1:size(target_A, 2)
            differences = 0.0
            denominator = 0.0
            for r in 1:size(target_A, 1)

                @inbounds x = target_A[r,target_index]
                @inbounds y = query_A[r,query_index]

                different = (x & y) < 16
                same = (x & 8 == 8) && x == y

                if different
                    differences += 1.0
                end

                if different || same
                    denominator += 1.0
                end
            end

            difference_array[target_index, query_index] = (differences / denominator)

        end
    end

    return difference_array
end

function get_min_float_indices(array_of_scores)

    min::Float64 = array_of_scores[1]
    indices::Array{Int64,1} = [1]

    for (i, score) in enumerate(array_of_scores)
        if i == 1
            min = score
            indices = [i]
        elseif score == min
            indices = vcat(indices, i)
        elseif score < min
            min = score
            indices = [i]
        end
    end

    return min, indices
end

function get_max_int_indices(array_of_scores)

    max::Int64 = array_of_scores[1]
    indices::Array{Int64,1} = [1]

    for (i, score) in enumerate(array_of_scores)
        if i == 1
            max = score
            indices = [i]
        elseif score == max
            indices = vcat(indices, i)
        elseif score > max
            max = score
            indices = [i]
        end
    end

    return max, indices
end

function get_best_target_index(differences_V, completeness_V)
    distance_min, distance_min_indices = get_min_float_indices(differences_V)

    if length(distance_min_indices) > 1
        completenesses = completeness_V[distance_min_indices]
        completeness_max, completeness_max_indices = get_max_int_indices(completenesses)
        indx = distance_min_indices[completeness_max_indices[1]]
    else
        indx = distance_min_indices[1]
    end

    return indx
end

function get_SNPs(target_V, query_V)
    nuc_dict = Dict{UInt8, Char}(136 => 'A',
                                 72 => 'G',
                                 40 => 'C',
                                 24 => 'T',
                                 192 => 'R',
                                 160 => 'M',
                                 144 => 'W',
                                 96 => 'S',
                                 80 => 'K',
                                 48 => 'Y',
                                 224 => 'V',
                                 176 => 'H',
                                 208 => 'D',
                                 112 => 'B',
                                 240 => 'N',
                                 244 => '-',
                                 242 => '?')
    SNPs = Array{String, 1}()
    for r in 1:length(target_V)
        @inbounds x = target_V[r]
        @inbounds y = query_V[r]
        different = (x & y) < 16
        if different
            nuc_T = nuc_dict[x]
            nuc_Q = nuc_dict[y]
            snp = string(r) * nuc_T * nuc_Q
            push!(SNPs, snp)
        end
    end

    return join(SNPs, ";")
end

function closest(target_file, query_file, outfile)
    T_names = get_fasta_descriptions(target_file)
    Q_names = get_fasta_descriptions(query_file)

    T_A = populate_byte_array(target_file)
    Q_A = populate_byte_array(query_file)

    differences = get_difference_matrix(T_A, Q_A)

    target_completeness = score_alignment(T_A)

    open(outfile, "w") do io

        println(io, "query,target,distance,SNPs")

        for query_index in 1:size(differences, 2)
            best_indx = get_best_target_index(differences[:,query_index], target_completeness)
            distance = differences[best_indx, query_index]
            SNPs = get_SNPs(Q_A[:,query_index], T_A[:,best_indx])
            Q_name = Q_names[query_index]
            T_name = T_names[best_indx]

            println(io, Q_name * "," * T_name * "," * string(round(distance, digits = 9)) * "," * SNPs)
        end
    end

end
