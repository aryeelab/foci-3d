# Input:    Pairtools pairs file
# Output:   One line per fragment (i.e. each pair becomes two lines)
#           Each fragment is represented by its chromosome, midpoint and length


import sys

input_file = sys.argv[1]  # Get the input file name from command line argument
output_file = sys.argv[2]  # Get the output file name from command line argument

# Add a header to the output file
with open(output_file, "w") as outfile:
    outfile.write("chrom\tmidpoint\tlength\n")

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        # Skip header lines starting with '#'
        if line.startswith("#"):
            continue
        
        # Split the line into columns
        columns = line.strip().split("\t")
        
        # Extract first fragment (chrom1, pos51, pos31)
        chrom1, pos51, pos31 = columns[1], int(columns[8]), int(columns[10])
        start1, end1 = min(pos51, pos31), max(pos51, pos31)
        # Convert to midpoint, length
        midpoint1 = (start1 + end1) // 2
        length1 = end1 - start1 + 1
        outfile.write(f"{chrom1}\t{midpoint1}\t{length1}\n")
        
        # Extract second fragment (chrom2, pos52, pos32)
        chrom2, pos52, pos32 = columns[3], int(columns[9]), int(columns[11])
        start2, end2 = min(pos52, pos32), max(pos52, pos32)
        # Convert to midpoint, length
        midpoint2 = (start2 + end2) // 2
        length2 = end2 - start2 + 1
        outfile.write(f"{chrom2}\t{midpoint2}\t{length2}\n")
