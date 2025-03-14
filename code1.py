import pysam

def process_bam(file_path, target_region, quality_threshold, mismatch_threshold):
    """
    Processes a BAM file to:
    1. Count total on-target reads.
    2. Count reads with mapping quality > quality_threshold.
    3. Filter reads with mismatch rate > mismatch_threshold.
    """
    on_target_reads = 0
    high_quality_reads = 0
    filtered_reads = []

    # Open BAM file
    bamfile = pysam.AlignmentFile(file_path, "rb")

    for read in bamfile.fetch(region=target_region):
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        # Skip reads with only "NNNN"
        if read.query_sequence and set(read.query_sequence) == {"N"}:
            continue

        # Check if the read is on-target
        on_target_reads += 1

        # Check mapping quality
        if read.mapping_quality > quality_threshold:
            high_quality_reads += 1

        # Calculate mismatch rate
        if read.has_tag("NM"):  # NM tag gives the number of mismatches
            mismatch_rate = read.get_tag("NM") / read.query_length
            if mismatch_rate > mismatch_threshold:
                filtered_reads.append(read.query_name)

    bamfile.close()

    return on_target_reads, high_quality_reads, filtered_reads

# Parameters
bam_file_path = "input.bam"
target_region = "chr1:10000-20000"  # Example target region (e.g., a genomic region)
mapping_quality_threshold = 20
mismatch_rate_threshold = 0.1

# Process BAM file
on_target, high_quality, filtered = process_bam(bam_file_path, target_region, mapping_quality_threshold, mismatch_rate_threshold)

# Output results
print(f"Total on-target reads: {on_target}")
print(f"Reads with mapping quality > {mapping_quality_threshold}: {high_quality}")
print(f"Filtered reads with mismatch rate > {mismatch_rate_threshold}: {len(filtered)}")