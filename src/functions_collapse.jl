using StatsBase

include("functions_align.jl")

# compare sites (rows) until a difference is found, for all pairs of sequences (columns)
# - break as soon as you know they're different - don't need any more information from this pair
# NB - if I declare the type (UInt8) of the input array, will things be faster?
# NB - Consider using views for slices (https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-tips-1)
function is_same(nuc_bit_array)
    same_array_size = size(nuc_bit_array, 2)
    # println(same_array_size)
    same_array = trues(same_array_size, same_array_size)
    # comb_test=0
    for a in 1:size(nuc_bit_array,2) - 1
        for b in (a + 1):size(nuc_bit_array,2)
            # comb_test+=1
            for r in Iterators.reverse(1:size(nuc_bit_array,1))
                @inbounds x = nuc_bit_array[r,a]
                @inbounds y = nuc_bit_array[r,b]
                different = (x & y) < 16

                # store the different True/False result in an array:
                if different
                    same_array[a,b] = !different
                    same_array[b,a] = !different
                    break
                end
            end
        end
    end
    # println(comb_test)
    return same_array
end

# # attempt at parallel version of the above
# using SharedArrays
# using Distributed
#
# function parallel_is_same(nuc_bit_array)
#     same_array_size = size(nuc_bit_array, 2)
#
#     # initialises as falses so have to change this:
#     same_array = SharedArray{Bool,2}((same_array_size,same_array_size))
#
#     for i in 1:size(same_array, 2)
#         same_array[i,i] = true
#     end
#
#     # comb_test=0
#     @sync @distributed for a in 1:size(nuc_bit_array,2) - 1
#         for b in (a + 1):size(nuc_bit_array,2)
#             # comb_test+=1
#             all_same = true
#             for r in Iterators.reverse(1:size(nuc_bit_array,1))
#                 @inbounds x = nuc_bit_array[r,a]
#                 @inbounds y = nuc_bit_array[r,b]
#                 different = (x & y) < 16
#
#                 # store the different True/False result in an array:
#                 if different
#                     all_same = false
#                     break
#                 end
#             end
#             same_array[a,b] = all_same
#             same_array[b,a] = all_same
#         end
#     end
#     # println(comb_test)
#     return sdata(same_array)
# end

function get_one_set_from_view(bool_view)
    # one subset of rows of the view
    A = Array{Int64,1}()
    original_row_numbers = bool_view.indices[1]

    for (index, value) in enumerate(bool_view)
        if value
            push!(A, original_row_numbers[index])
        end
    end
    return Set(A)
end

function get_sets_from_view(bool_view)
    S = Set{Set}()
    for i in 1:size(bool_view, 2)
        push!(S, get_one_set_from_view(view(bool_view, :, i)))
    end
    return collect(S)
end

function check_the_view_is_good(a_bool_array)
    # because these logical arrays have symmetry only need
    # to check one triangle, but this is harder to code, so haven't
    for i in a_bool_array
        if !i
            return false
        end
    end
    return true
end

# 1) test all sets for being supersets, and throw them out if they are
# NB - this is 1/10 of the run time and allocates a lot of memory in total (GBs - not all at once)
# NB could probably be better optimised
function get_subsets(array_of_sets)

    A = Array{BitSet,1}()
    sizehint!(A, length(array_of_sets))

    for i in 1:size(array_of_sets, 1)
        valid = true
        for j in 1:size(array_of_sets, 1)
            if i == j
                continue
            end
            # if array_of_sets[i] is a superset of comparison set, then break, else push to A
            # (note the reverse logic here, because there is no superset function)
            if issubset(array_of_sets[j], array_of_sets[i])
                valid = false
                break
            end
        end

        if valid
            push!(A, array_of_sets[i])
        end
    end

    return A
end

function get_one_truefalse_set_from_vector(bool_vector)
    # get the indices of Trues in bool_vector as a BitSet
    # bool vector must be 1-dimensional
    T = Array{Int64,1}()
    F = Array{Int64,1}()

    for i in 1:length(bool_vector)
        if bool_vector[i]
            push!(T, i)
        else
            push!(F, i)
        end
    end

    BT = BitSet(T)
    BF = BitSet(F)

    return BT, BF
end

function get_truefalse_row_indices(bool_array)
    # will return two arrays of BitSets which are:
    # 1) the rows in bool_array that are true, for every column in bool_array - matches between pairs of sequences
    # 2) the rows in bool_array that are false, for every column in bool_array - mismatches between pairs of sequences

    AT = Array{BitSet,1}(undef, size(bool_array, 2))
    AF = Array{BitSet,1}(undef, size(bool_array, 2))

    for i in 1:size(bool_array, 2)
        AT[i], AF[i] = get_one_truefalse_set_from_vector(view(bool_array, :, i))
    end

    return AT, AF
end

function get_sets_in_one_go(the_whole_bool_array)
    # returns an array of bitsets which are the sets

    good_views = []
    bad_views = []

    the_whole_thing_is_fine = check_the_view_is_good(the_whole_bool_array)

    if the_whole_thing_is_fine
        return [BitSet(collect(1:size(the_whole_bool_array, 2)))]
    else
        push!(bad_views, the_whole_bool_array)
    end

    while length(bad_views) > 0
        new_bad_views = []
        for bv in bad_views
            new_bitsets = get_subsets(get_sets_from_view(bv))
            for bitset in new_bitsets
                set = collect(bitset)
                new_view = @view the_whole_bool_array[set, set]
                this_view_is_good = check_the_view_is_good(new_view)
                if this_view_is_good
                    push!(good_views, new_view)
                else
                    push!(new_bad_views, new_view)
                end
            end
        end
        bad_views = new_bad_views
    end

    A = Array{BitSet,1}()

    for gv in good_views
        inds = gv.indices[2]
        push!(A, BitSet(inds))
    end

    return Set(A)
end

function write_highest_scoring_unused_seq(final_sets, retained, whole_nuc_bit_array, whole_ID_array, alignment_out, append_IDs)

    all_scores = score_alignment(whole_nuc_bit_array)

    # convert the set to an Array
    final_sets_as_array = collect(final_sets)

    # split up singletons and others to relieve the burden on sorting?
    singletons = Array{BitSet,1}()
    non_single = Array{BitSet,1}()

    for set in final_sets_as_array
        if length(set) == 1
            push!(singletons, set)
        else
            push!(non_single, set)
        end
    end

    # sort non-singletons by length in increasing order
    sort!(non_single, by = length)

    # and here are singletons + sorted longer sets, concatenated together
    sorted_final_sets = vcat(singletons, non_single)

    # initiate a set that we will use as a check for whether a particular sequence
    # has been used (to represent a set) yet
    # used_names = Set{String}()
    # NB doing this with a dictionary instead now
    used_names = Dict{String, Array{String, 1}}()

    open(alignment_out, "w") do io
        for (i, set) in enumerate(sorted_final_sets)

            # cs is the columns in the alignment that this set represents
            cs = collect(set)

            # If it's a singleton it can't have been so we can just get on with things.
            if length(cs) == 1

                id = ">" * whole_ID_array[cs[1]]
                seq = get_seq_from_1D_byte_array(whole_nuc_bit_array[:,cs[1]])

                # write fasta header and sequence:
                println(io, id)
                println(io, seq)

                used_names[whole_ID_array[cs[1]]] = []

            # otherwise get the highest scoring sequence
            else

                set_scores = all_scores[cs]
                set_nuc_array = whole_nuc_bit_array[:,cs]
                set_IDs = whole_ID_array[cs]

                # max_score, max_indices = get_max_indices(set_scores)
                max_indices = sortperm(set_scores, rev = true)

                indx = 1
                best_column = max_indices[indx]
                best_ID = set_IDs[best_column]

                # A while loop to check that the highest scoring sequence hasn't
                # already been used to represent a set - if it has, go through
                # each next best sequence in turn.
                while in(best_ID, keys(used_names)) && indx < length(max_indices)
                    indx += 1
                    best_column = max_indices[indx]
                    best_ID = set_IDs[best_column]
                end

                # if, after the iteration above, all of the IDs in this set are already
                # representing sets (are in used_names), then we can just skip this set (because
                # all it's members are assigned to a set anyway).
                if in(best_ID, keys(used_names))
                    continue
                end

                best_seq = set_nuc_array[:,best_column]
                other_IDs = setdiff(set_IDs, [best_ID])

                used_names[best_ID] = collect(other_IDs)

                if !append_IDs
                    id = ">" * string(best_ID)
                else
                    id = ">" * join(vcat(best_ID, other_IDs), "|")
                end

                seq = get_seq_from_1D_byte_array(best_seq)

                # write fasta header and sequence:
                println(io, id)
                println(io, seq)

            end

        end

    # add retained sequences to the alignment if they aren't included already,
    # and add them to used names so that don't get duplicated in the mapping
    for seqname in retained
        if !in(seqname, keys(used_names))

            indx = findfirst(isequal(seqname), whole_ID_array)

            if indx == nothing
                println(stderr, "warning: " * seqname * " can't be retained because it isn't in the input alignment")
                continue
            end

            if seqname != whole_ID_array[indx]
                throw("bad indexing when forced to retain " * seqname)
            end

            id = ">" * whole_ID_array[indx]
            seq = get_seq_from_1D_byte_array(whole_nuc_bit_array[:,indx])

            println(io, id)
            println(io, seq)

            used_names[seqname] = []
        end
    end

    close(io)
    end

    return used_names
end

function get_redundant_seq_to_tip_relationships(representative_to_set)
    #= representative_to_set is a dict with representative
       haplotypes to be written to file as keys, and the other
       members of the set as values.

       Want to get a map from each possible member => representatives
    =#

    D = Dict{String, Array{String, 1}}()

    for (key, value) in representative_to_set
        # key is the representative sequence for this set
        # value is an array of other sequences in the set
        for v in value
            # if this guy is already a representative of a set,
            # then skip it
            if in(v, keys(representative_to_set))
                continue
            end

            if in(v, keys(D))
                D[v] = vcat(D[v], [key])
            else
                D[v] = [key]
            end
        end
    end

    return D
end

function write_redundant_seq_to_tip_relationships(D, filepath)
    # write a table of sequences that are not included in the
    # collapsed alignment, with information about their place(s)
    # on the tree i.e., for each excluded seq, what tip it can
    # be placed at
    open(filepath, "w") do io
        println(io, "sequence,tips")
        for (key, value) in D
            println(io, key * "," * join(value, "|"))
        end
    close(io)
    end
end

function write_tip_to_redundant_seq_relationships(D, filepath, retained)
    # write the converse file to write_redundant_seq_to_tip_relationships()
    # i.e., tips of the tree and what other seqs they represent
    open(filepath, "w") do io
        println(io, "tip,sequences")
        for (key, value) in D
            # subtract retained sequences from the set that is represented
            value = setdiff(value, retained)
            # ignore sequences that represent no other sequences:
            if length(value) == 0
                continue
            end
            println(io, key * "," * join(value, "|"))
        end
    close(io)
    end
end

function TEST_get_levels(collection_of_bitsets)
    #=
    return the set of all rows
    that make it into the output - if things
    have worked, then the length of this should
    be the same as the number of sequences - all
    sequences are included in the dataset (at least once)
    =#
    levels = []
    for x in collection_of_bitsets
        for y in x
            push!(levels, y)
        end
    end
    return Set(levels)
end

function TEST_compare_within_sets(array_of_bitsets, nuc_bit_array)
    #=
    compare sequences within sets using the original data -
    =#

    function TEST_all_same(my_small_array)
        for a in 1:size(my_small_array,2) - 1
            for b in (a + 1):size(my_small_array,2)
                # comb_test+=1
                for r in 1:size(my_small_array,1)
                    @inbounds x = my_small_array[r,a]
                    @inbounds y = my_small_array[r,b]
                    different = (x & y) < 16

                    # store the different True/False result in an array:
                    if different
                        return false
                    end
                end
            end
        end

        return true
    end

    samples = sample(collect(1:length(array_of_bitsets)), length(array_of_bitsets), replace = false)
    fails = []
    for smp in samples
        bitset = array_of_bitsets[smp]
        columns = collect(bitset)
        if length(columns) > 1
            # println(columns)
            # A = @view nuc_bit_array[:,columns]
            sub_nuc_bit_array = nuc_bit_array[:,columns]
            # println(size(sub_nuc_bit_array))
            test = TEST_all_same(sub_nuc_bit_array)
            # println(test)
            # println()
            if !test
                push!(fails, columns)
            end
        end
    end

    if length(fails) > 0
        return false
        # return false, fails
    end

    return true
end
