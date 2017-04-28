#reference contig

type ReferenceContigs
    count::Int64
    names::Array{String}
    sizes::Array{Int64}
    offsets::Array{Int64}
    
    function ReferenceContigs(count, names, sizes)
        new(count, names, sizes, [sum(sizes[1:i-1]) for i in 1:length(sizes)])
    end
end

using GZip

import Base: eof, close, position
export BamReader, close, value, eof, advance!, eachposition

type BamReader
    bamStream
    readOrientation #useReverseReads::Bool
    done::Bool
    position::Int64
    isforward::Bool
    chrom::Int32
    contigs::ReferenceContigs
end

function BamReader(bamFileName::String, readOrientation)
    f = GZip.open(bamFileName)

    # this is a BAM file
    code = read(f, UInt8, 4)
    @assert code == b"BAM\1"

    # get through the header data
    l_text = read(f, Int32)
    skip(f, l_text)

    # build a reference contig with the bam header. 
    n_ref = read(f, Int32)
    sizes = Int32[]
    names = String[]
    for j in 1:n_ref
        l_name = read(f, Int32)
        refName = convert(String, read(f, UInt8, l_name)[1:end-1]) # ignore the null terminator
        l_ref = read(f, Int32)
        push!(sizes, l_ref)
        push!(names, refName)
    end
    
    contigs = ReferenceContigs(n_ref, names, sizes)

    r = BamReader(f, readOrientation, false, 1, -1, contigs)
    advance!(r)
    r
end

close(reader::BamReader) = GZip.close(reader.bamStream)
value(reader::BamReader) = 1
position(reader::BamReader) = reader.position
isforward(reader::BamReader) = reader.isforward
chrom(reader::BamReader) = get(reader.contigs.names, reader.chrom,"n/a")
eof(reader::BamReader) = reader.position == -1

function advance!(r::BamReader)
    f = r.bamStream
    while !r.done
        if peek(f) == -1 # eof does not work on the BAM files either in C++ or here (BAM vs. gzip issue?)
            r.done = true
            r.position = -1
            return
        end

        buf = Array(Int32, 6) # [block_size, refID, pos, bin_mq_nl, flag_nc, l_seq]
        gzread(f, pointer(buf), 24) # moving the pointer by 24 bytes 4 bytesx6=24
        block_size = buf[1]
        refID = buf[2] + 1 # the reference contig this read maps to
        lseq = buf[6]

        forward = (buf[5] & 1048576) == 0 # see if we are reverse complemented
        skip(f, block_size-20) # skip the rest of the entry
        
        #writing out some informations
        r.chrom = refID
        #get the position and convert to 1 based indexing
        if forward
            r.position = buf[3] + 1
        else
            r.position = buf[3] + 1 + buf[6]
        end
        r.isforward = forward
        
        #work-around here
        if get(r.contigs.sizes, refID, 1) < r.position
            r.position = buf[3] + 1
        end

        # break if we found a read in the right direction
        if refID != 0 && (r.readOrientation == :any || (forward && r.readOrientation == :forward) || (!forward && r.readOrientation == :reverse))
            return
        end
    end
end