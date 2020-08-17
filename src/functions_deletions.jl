function read_del_file(file_path)
    # input file is in the format:
    # start (1-based), reference allele
    # e.g.:
    # 1605,TTT
    A = Array{Tuple{Int64, String}}(undef,0)
    for line in eachline(file_path)
        l = split(line, ",")
        start = parse(Int64, l[1])
        ref = l[2]
        t = (start, ref)
        push!(A, t)
    end
    return A
end

function type_deletions(nuc_bit_array, del_tuple_array)
    # type specific deletions in an alignment
    byte_dict = Dict{Char, UInt8}('A' => 136,
                                    'C' => 72,
                                    'G' => 40,
                                    'T' => 24,
                                    'R' => 192,
                                    'M' => 160,
                                    'W' => 144,
                                    'S' => 96,
                                    'K' => 80,
                                    'Y' => 48,
                                    'V' => 224,
                                    'H' => 176,
                                    'D' => 208,
                                    'B' => 112,
                                    'N' => 240,
                                    '-' => 244,
                                    '?' => 242)

    # A will be returned: it is an array of strings, one column for every sequence
    # in the alignment, one row for every deletion to be typed. Entries in cells can
    # be "ref" (Wuhan-Hu-1 allele), "del" ('-' at all sites) or "X" (any other allele)
    A_genotypes_table = Array{String,2}(undef, length(del_tuple_array), size(nuc_bit_array, 2))

    for sequence_column in 1:size(nuc_bit_array, 2)
        for (del_row,t) in enumerate(del_tuple_array)
            start = t[1]
            ref = t[2]

            del_array = fill(byte_dict['-'], length(ref))
            ref_array = Array{UInt8,1}(undef, length(ref))

            for (j,c) in enumerate(ref)
                ref_array[j] = byte_dict[c]
            end

            QUERY = nuc_bit_array[start:start + length(ref) - 1, sequence_column]

            if QUERY == del_array
                genotype = "del"

            elseif QUERY == ref_array
                genotype = "ref"

            else
                genotype = "X"
            end

            A_genotypes_table[del_row, sequence_column] = genotype
        end
    end

    return A_genotypes_table
end

function type_deletions_and_append_as_SNP(nuc_bit_array, del_tuple_array)
    # type specific deletions in an alignment
    # and append(/prepend?) them as genotypes
    byte_dict = Dict{Char, UInt8}('A' => 136,
                                    'C' => 72,
                                    'G' => 40,
                                    'T' => 24,
                                    'R' => 192,
                                    'M' => 160,
                                    'W' => 144,
                                    'S' => 96,
                                    'K' => 80,
                                    'Y' => 48,
                                    'V' => 224,
                                    'H' => 176,
                                    'D' => 208,
                                    'B' => 112,
                                    'N' => 240,
                                    '-' => 244,
                                    '?' => 242)

    # A will be returned: it is an array of strings, one column for every sequence
    # in the alignment, one row for eery deletion to be typed. Entries in cells can
    # be "ref" (Wuhan-Hu-1 allele), "del" ('-' at all sites) or "X" (any other allele)
    A_genotypes_table = Array{String,2}(undef, length(del_tuple_array), size(nuc_bit_array, 2))

    # this is for the appending of a SNP:
    A_dels_as_SNPs_array = Array{UInt8,2}(undef, length(del_tuple_array), size(nuc_bit_array, 2))

    for sequence in 1:size(nuc_bit_array, 2)
        for (del,t) in enumerate(del_tuple_array)
            start = t[1]
            ref = t[2]

            del_array = fill(byte_dict['-'], length(ref))
            ref_array = Array{UInt8,1}(undef, length(ref))

            for (j,c) in enumerate(ref)
                ref_array[j] = byte_dict[c]
            end

            QUERY = nuc_bit_array[start:start + length(ref) - 1, sequence]

            if QUERY == del_array
                nuc = byte_dict['C']
                genotype = "del"

            elseif QUERY == ref_array
                nuc = byte_dict['A']
                genotype = "ref"

            else
                nuc = byte_dict['N']
                genotype = "X"
            end

            A_genotypes_table[del, sequence] = genotype
            A_dels_as_SNPs_array[del, sequence] = nuc
        end
    end

    # vcat A_dels_as_SNPs_array on the bottom of nuc_bit_array to get an alignment
    # with deletions added as SNPs
    appended_array = vcat(nuc_bit_array, A_dels_as_SNPs_array)

    return appended_array, A_genotypes_table
end
