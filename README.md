# jlBAM
Very simple BAM file parser module for Julia language (v5.0).

This currently parses through a file with the BAM format and streams DNA reads information.
As this is a snippet of code used for my research, it only does one speicific thing, which is 
to stream **5'-end positions** of both forward and reverse reads. This is little different from 
how other parsers do because they tend to output the leftmost poistions of reads (which are 5'-end
for foward reads and 3'-end for reverse reads).

One can easily change the code to get other types of information about reads if necessary.

## Sample usage:
``` Julia
b = BamReader("test.bam", :any) #say :forward or :reverse to focus on a specific direction.
while !eof(b)
    print(position(b), chrom(b), isforward(b))
end
```

