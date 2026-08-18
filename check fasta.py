from pathlib import Path

folder = Path(__file__).parent

files = list(folder.glob("*.txt"))

print("Files found:", files)

if not files:
    print("❌ FASTA file nahi mili")
else:
    file = files[0]
    print("✅ File found:", file.name)

    from Bio import SeqIO

file = "sequence.fasta.txt"

sequences = list(SeqIO.parse(file, "fasta"))

print("Total sequences:", len(sequences))
print("First sequence ID:", sequences[0].id)
print("First sequence length:", len(sequences[0].seq))
print("First sequence:", sequences[0].seq)
gc_values = []

for seq in sequences:
    dna = str(seq.seq)

    g = dna.count("G")
    c = dna.count("C")

    gc = ((g + c) / len(dna)) * 100

    gc_values.append(gc)

print("Average GC Content:", sum(gc_values) / len(gc_values))
print("Minimum GC:", min(gc_values))
print("Maximum GC:", max(gc_values))
import matplotlib.pyplot as plt

plt.hist(gc_values, bins=20)

plt.title("GC Content Distribution of Human Promoters")
plt.xlabel("GC Content (%)")
plt.ylabel("Number of Promoters")

plt.tight_layout()
plt.show()
total_bases = len(sequences) * 300

A = 0
T = 0
G = 0
C = 0

for seq in sequences:
    dna = str(seq.seq).upper()

    A += dna.count("A")
    T += dna.count("T")
    G += dna.count("G")
    C += dna.count("C")

print("A:", (A / total_bases) * 100, "%")
print("T:", (T / total_bases) * 100, "%")
print("G:", (G / total_bases) * 100, "%")
print("C:", (C / total_bases) * 100, "%")

tata_count = 0

for seq in sequences:
    dna = str(seq.seq).upper()

    if "TATA" in dna:
        tata_count += 1

print("Sequences containing TATA motif:", tata_count)
print("Percentage:", (tata_count / len(sequences)) * 100, "%")
ccaat_count = 0

for seq in sequences:
    dna = str(seq.seq).upper()

    if "CCAAT" in dna:
        ccaat_count += 1

print("Sequences containing CCAAT motif:", ccaat_count)
print("Percentage:", (ccaat_count / len(sequences)) * 100, "%")

gc_box_count = 0

for seq in sequences:
    dna = str(seq.seq).upper()

    if "GGGCGG" in dna:
        gc_box_count += 1

print("Sequences containing GC box motif:", gc_box_count)
print("Percentage:", (gc_box_count / len(sequences)) * 100, "%")

import matplotlib.pyplot as plt

motifs = ["TATA", "CCAAT", "GC box"]
counts = [tata_count, ccaat_count, gc_box_count]

plt.bar(motifs, counts)

plt.title("Promoter Motif Distribution")
plt.xlabel("Motif")
plt.ylabel("Number of Promoters")

plt.tight_layout()
plt.show()
tata_positions = []

for seq in sequences:
    dna = str(seq.seq).upper()

    position = dna.find("TATA")

    if position != -1:
        tata_positions.append(position)

print("Number of TATA motifs:", len(tata_positions))
print("First 10 TATA positions:", tata_positions[:10])
print("Average TATA position:", sum(tata_positions) / len(tata_positions))

ccaat_positions = []

for seq in sequences:
    dna = str(seq.seq).upper()
    position = dna.find("CCAAT")

    if position != -1:
        ccaat_positions.append(position)

print("CCAAT-positive sequences:", len(ccaat_positions))
print("First 10 positions:", ccaat_positions[:10])
print("Average CCAAT position:", sum(ccaat_positions) / len(ccaat_positions))

gc_box_positions = []

for seq in sequences:
    dna = str(seq.seq).upper()
    position = dna.find("GGGCGG")

    if position != -1:
        gc_box_positions.append(position)

print("GC-box-positive sequences:", len(gc_box_positions))
print("First 10 positions:", gc_box_positions[:10])
print("Average GC-box position:", sum(gc_box_positions) / len(gc_box_positions))

import matplotlib.pyplot as plt

motif_names = ["TATA", "CCAAT", "GC box"]

average_positions = [
    sum(tata_positions) / len(tata_positions),
    sum(ccaat_positions) / len(ccaat_positions),
    sum(gc_box_positions) / len(gc_box_positions)
]

plt.bar(motif_names, average_positions)

plt.title("Average Position of Promoter Motifs")
plt.xlabel("Motif")
plt.ylabel("Average Position in 300-bp Sequence")

plt.tight_layout()
plt.show()