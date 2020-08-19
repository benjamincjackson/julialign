#=
See - http://ape-package.ird.fr/misc/BitLevelCodingScheme.html for the bit-level
coding scheme :

Nucleotide	    IUPAC code	 Bit-level code	  Value
---------------------------------------------------
A	                  A	      10001000	       136
G	                  G	      01001000	        72
C	                  C	      00101000	        40
T	                  T	      00011000	        24
A or G	              R	      11000000	       192
A or C	              M	      10100000	       160
A or T	              W	      10010000	       144
G or C	              S	      01100000	        96
G or T	              K	      01010000	        80
C or T	              Y	      00110000	        48
A or G or C	          V	      11100000	       224
A or C or T	          H	      10110000	       176
A or G or T	          D	      11010000	       208
G or C or T	          B	      01110000	       112
A or G or C or T      N	      11110000	       240
Alignment gap        (–)	  00000100	         4
Unknown character    (?)	  00000010	         2
---------------------------------------------------


This has been adapated (last two rows have changed) to
allow matches between alignment gaps / ? and any nucleotide:

Nucleotide	    IUPAC code	 Bit-level code	  Value
---------------------------------------------------
A	                  A	      10001000	       136
G	                  G	      01001000	        72
C	                  C	      00101000	        40
T	                  T	      00011000	        24
A or G	              R	      11000000	       192
A or C	              M	      10100000	       160
A or T	              W	      10010000	       144
G or C	              S	      01100000	        96
G or T	              K	      01010000	        80
C or T	              Y	      00110000	        48
A or G or C	          V	      11100000	       224
A or C or T	          H	      10110000	       176
A or G or T	          D	      11010000	       208
G or C or T	          B	      01110000	       112
A or G or C or T      N	      11110000	       240
Alignment gap        (–)	  11110100	       244  ###### CHANGED
Unknown character    (?)	  11110010	       242  ###### CHANGED
---------------------------------------------------


Function	               C code	                      Value returned
------------------------------------------------------------------------------------------
KnownBase(a)	           a & 8	                      8 if a is known surely
IsAdenine(a)	           a == 136	                      1 if a is adenine
IsGuanine(a)	           a == 72	                      1 if a is guanine
IsCytosine(a)       	   a == 40	                      1 if a is cytosine
IsThymine(a)	           a == 24	                      1 if a is thymine
IsPurine(a)	               a & 55	                      0 if a is a purine
IsPyrimidine(a)	           a & 199	                      0 if a is a pyrimidine
DifferentBase(a, b) 	   (a & b) < 16	                  1 if a and b are different surely
SameBase(a, b)	           KnownBase(a) && a == b	      1 if a and b are the same surely
-------------------------------------------------------------------------------------------
=#

struct fasta_record
    id::String
    description::String
    seq::String
end

function get_alignment_dimensions(filepath::AbstractString)

    n = 0 # number of sequences
    l = 0 # length of sequence
    first = true

    open(filepath, "r") do io
        while !eof(io)

            if n > 1
                first = false
            end

            line = readline(io)

            if string(line[1]) == ">"
                n+=1
            end

            if first && string(line[1]) != ">"
                l+=length(line)
            end
        end

    end

    return l, n
end

function get_fasta_descriptions(filepath::AbstractString)

    A = Array{String,1}()

    open(filepath, "r") do io
        while !eof(io)
            l = readline(io)
            if string(l[1]) == ">"
                description = lstrip(l, '>')
                id = split(description)[1]
                push!(A, description)
            end
        end
    close(io)
    end

    return A
end

function read_fasta_alignment(filepath::AbstractString, channel::Channel)

    open(filepath, "r") do io

        first = true

        seq_buffer = String("")
        id = String("")
        description = String("")

        while !eof(io)

            l = readline(io)

            if first
                if string(l[1]) != ">"
                    throw("badly formed fasta file")
                end

                description = lstrip(l, '>')
                id = split(description)[1]

                first = false

            elseif eof(io)

                seq_buffer = seq_buffer * uppercase(l)

                fr = fasta_record(id, description, seq_buffer)
                put!(channel, fr)

                break

            elseif string(l[1]) == ">"

                fr = fasta_record(id, description, seq_buffer)
                put!(channel, fr)

                description = lstrip(l, '>')
                id = split(description)[1]
                seq_buffer = String("")

            else

                seq_buffer = seq_buffer * uppercase(l)

            end
        end

    close(io)
    end
end

function populate_byte_array(filepath::AbstractString)
    byte_dict = Dict{Char, UInt8}('A' => 136,
                                    'G' => 72,
                                    'C' => 40,
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

    alignment_dim = get_alignment_dimensions(filepath)
    height = alignment_dim[1]
    width = alignment_dim[2]

    ch = Channel{fasta_record}((channel_arg) -> read_fasta_alignment(filepath, channel_arg))

    A = Array{UInt8,2}(undef, height, width)

    j = 1
    for record in ch
        sequence = record.seq
        i = 1
        for nuc in sequence
            byte = byte_dict[nuc]
            A[i, j] = byte
            i+=1
        end
        j+=1
    end

    return A
end

function get_seq_from_1D_byte_array(byte_V)
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

    char_V = Array{Char, 1}(undef, length(byte_V))

    for i in 1:size(byte_V, 1)
        char_V[i] = nuc_dict[byte_V[i]]
    end

    return join(char_V)
end

function score_alignment_column(nuc_bit_sequence)
    score_dict = Dict{UInt8, Int64}(136 => 12,
                                     72 => 12,
                                     40 => 12,
                                     24 => 12,
                                     192 => 6,
                                     160 => 6,
                                     144 => 6,
                                     96 => 6,
                                     80 => 6,
                                     48 => 6,
                                     224 => 4,
                                     176 => 4,
                                     208 => 4,
                                     112 => 4,
                                     240 => 3,
                                     244 => 3,
                                     242 => 3)

    score = 0::Int64
    for nuc in nuc_bit_sequence
        score+=score_dict[nuc]
    end
    return score
end


function score_alignment(nuc_bit_array)

    A = fill(0::Int64, size(nuc_bit_array, 2))

    for i in 1:size(nuc_bit_array, 2)
        A[i] = score_alignment_column(view(nuc_bit_array, :, i))
    end

    return A
end
