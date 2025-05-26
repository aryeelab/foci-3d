# Input:    Pairtools pairs file
# Output:   One line per fragment (i.e. each pair becomes two lines)
#           Each fragment is represented by its chromosome, midpoint and length
#
# This script can dynamically determine column indices from the header line
# or fall back to default indices for backward compatibility.
#
# The script looks for a header line starting with "#columns:" to determine
# the positions of the required columns (chrom1, chrom2, pos51, pos52, pos31, pos32).
# If the header line is not found or does not contain all required columns,
# the script falls back to default indices.


import sys

def get_column_indices(header_line):
    """
    Parse the header line to get the indices of the required columns.

    Args:
        header_line (str): The header line starting with '#columns:'

    Returns:
        dict: A dictionary mapping column names to their indices

    Raises:
        ValueError: If required columns are missing from the header
    """
    # Remove the '#columns:' prefix and split by whitespace
    columns = header_line.replace('#columns:', '').strip().split()

    # Create a dictionary mapping column names to their indices
    column_indices = {col: idx for idx, col in enumerate(columns)}

    # Check if all required columns are present
    required_columns = ['chrom1', 'chrom2', 'pos51', 'pos52', 'pos31', 'pos32']
    missing_columns = [col for col in required_columns if col not in column_indices]

    if missing_columns:
        raise ValueError(f"Required columns missing from header: {', '.join(missing_columns)}")

    return column_indices

input_file = sys.argv[1]  # Get the input file name from command line argument
output_file = sys.argv[2]  # Get the output file name from command line argument

# Initialize column indices with default values for backward compatibility
default_indices = {
    'chrom1': 1,
    'chrom2': 3,
    'pos51': 8,
    'pos52': 9,
    'pos31': 10,
    'pos32': 11
}

# First pass: look for the header line to get column indices
column_indices = default_indices.copy()
header_found = False

with open(input_file, "r") as infile:
    for line in infile:
        if line.startswith("#columns:"):
            print(f"Found header line: {line.strip()}", file=sys.stderr)
            try:
                column_indices = get_column_indices(line)
                header_found = True
                print(f"Using column indices: {column_indices}", file=sys.stderr)
            except ValueError as e:
                print(f"Warning: {e}. Using default column indices.", file=sys.stderr)
            break

if not header_found:
    print(f"No header line found. Using default column indices: {column_indices}", file=sys.stderr)

# Second pass: process the data
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # We don't write a header to the output file because the test script expects raw data
    # The test script will add the appropriate header

    # Process the input file
    line_count = 0
    data_line_count = 0

    for line in infile:
        line_count += 1

        # Skip header lines starting with '#'
        if line.startswith("#"):
            continue

        data_line_count += 1

        # Split the line into columns
        columns = line.strip().split("\t")

        # Check if the line has enough columns
        if len(columns) <= max(column_indices.values()):
            print(f"Warning: Line {line_count} has fewer columns than expected. Skipping: {line.strip()}", file=sys.stderr)
            continue

        try:
            # Extract both fragments
            chrom1 = columns[column_indices['chrom1']]
            pos51 = int(columns[column_indices['pos51']])
            pos31 = int(columns[column_indices['pos31']])
            start1, end1 = min(pos51, pos31), max(pos51, pos31)
            midpoint1 = (start1 + end1) // 2
            length1 = end1 - start1 + 1

            chrom2 = columns[column_indices['chrom2']]
            pos52 = int(columns[column_indices['pos52']])
            pos32 = int(columns[column_indices['pos32']])
            start2, end2 = min(pos52, pos32), max(pos52, pos32)
            midpoint2 = (start2 + end2) // 2
            length2 = end2 - start2 + 1

            # Always write the fragments in the same order (chrom1 first, then chrom2)
            # This ensures consistent output regardless of column order in the input file
            outfile.write(f"{chrom1}\t{midpoint1}\t{length1}\n")
            outfile.write(f"{chrom2}\t{midpoint2}\t{length2}\n")
        except (IndexError, ValueError) as e:
            print(f"Warning: Error processing line {line_count}: {e}. Skipping: {line.strip()}", file=sys.stderr)

    print(f"Processed {line_count} lines, {data_line_count} data lines", file=sys.stderr)
