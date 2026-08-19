from pathlib import Path
import pandas as pd

folder = Path(__file__).parent

# Folder mein .txt file automatically find karega
files = list(folder.glob("*.txt"))

print("TXT files found:", files)

if not files:
    print("❌ TXT file nahi mili")
else:
    file = files[0]

    df = pd.read_csv(
        file,
        sep=r"\s+",
        header=None,
        engine="python"
    )

    print("✅ File successfully loaded!")
    print("File:", file.name)
    print("Total rows:", len(df))
    print(df.head())
    df.columns = [
    "chromosome",
    "start",
    "end",
    "promoter_id",
    "score",
    "strand",
    "promoter_start",
    "promoter_end"
]

print(df.head())
print(df.info())

print("Missing values:")
print(df.isnull().sum())

print("\nDuplicate rows:", df.duplicated().sum())

print("\nData types:")
print(df.dtypes)
print("Total promotor:", len(df))
print("\nPromotors by chromosome:")
print(df["chromosome"]. value_counts().head(10))
print("\nPromotors by strand:")
print(df["strand"].value_counts())
df["promoter_length"] = df["promoter_end"] - df["promoter_start"] + 1

print("Promoter length statistics:")
print(df["promoter_length"].describe())
import matplotlib.pyplot as plt

#chromosome_counts = df["chromosome"].value_counts()

#chromosome_counts.plot(kind="bar")

#plt.title("Human Promoter Distribution by Chromosome")
#plt.xlabel("Chromosome")
#plt.ylabel("Number of Promoters")
#plt.xticks(rotation=45)
#plt.tight_layout()
#plt.show()
strand_counts= df["strand"].value_counts()
strand_counts.plot(kind="bar")
plt.title("Promotor Distributio by DNA Strand")
plt.xlabel("Strand")
plt.ylabel("Number of Promotors")
plt.tight_layout()
plt.show()
