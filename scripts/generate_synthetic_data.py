from pathlib import Path
import random

# Parameters
reference_length = 360
num_reads = 100
read_length = 30

# Output paths
ref_path = Path("synthetic_reads/test_ref.fasta")
fastq_path = Path("synthetic_reads/test_reads.fastq")

# Create synthetic reference
reference_seq = "ATGC" * (reference_length // 4)
ref_path.parent.mkdir(exist_ok=True)
with open(ref_path, "w") as f:
    f.write(f">synthetic_reference\n{reference_seq}\n")

# Generate synthetic reads evenly spaced across reference
with open(fastq_path, "w") as f:
    for i in range(num_reads):
        start = (i * reference_length) // num_reads
        read_seq = reference_seq[start:start + read_length]
        if len(read_seq) < read_length:
            read_seq += reference_seq[:read_length - len(read_seq)]
        f.write(f"@read{i}\n{read_seq}\n+\n{'I'*read_length}\n")

print("✅ Synthetic FASTA and FASTQ files generated.")
