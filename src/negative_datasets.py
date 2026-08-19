from pathlib import Path
import pandas as pd

project_folder = Path(__file__).parent

files = list(project_folder.glob("Hs_EPDnew_006_hg38.*"))

print("Files found:")
for file in files:
    print(file)

promoter_file = files[0]

promoters = pd.read_csv(
    promoter_file,
    sep=r"\s+",
    header=None
)

print("\nTotal promoter regions:", len(promoters))
print("\nFirst 5 promoter regions:")
print(promoters.head())

promoter_regions = promoters.iloc[:, [0, 1, 2]].copy()

promoter_regions.columns = ["chromosome", "start", "end"]

print("\nPromoter regions:")
print(promoter_regions.head())

print("\nNumber of promoter regions:", len(promoter_regions))

# Check the size of promoter regions
promoter_regions["length"] = (
    promoter_regions["end"] - promoter_regions["start"] + 1
)

print("\nPromoter length statistics:")
print(promoter_regions["length"].describe())
import random

# Generate random DNA sequences of 300 bp
def generate_random_sequence(length=300):
    bases = "ACGT"
    return "".join(random.choice(bases) for _ in range(length))


negative_sequences = []

for i in range(len(promoters)):
    seq = generate_random_sequence(300)
    negative_sequences.append(seq)

print("\nTotal negative sequences:", len(negative_sequences))
print("First negative sequence:")
print(negative_sequences[0])
print("Length:", len(negative_sequences[0]))
# Save negative sequences as FASTA
output_file = project_folder / "negative_sequences.fasta"

with open(output_file, "w") as f:
    for i, seq in enumerate(negative_sequences, start=1):
        f.write(f">negative_{i}\n")
        f.write(seq + "\n")

print("\nNegative FASTA file saved at:")
print(output_file)

# Read FASTA sequences
def read_fasta(file_path):
    sequences = []

    with open(file_path, "r") as f:
        sequence = ""

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
            else:
                sequence += line

        if sequence:
            sequences.append(sequence)

    return sequences


# Positive promoter sequences
positive_file = project_folder / "sequence.fasta.txt"
positive_sequences = read_fasta(positive_file)

print("\nPositive sequences:", len(positive_sequences))
print("Negative sequences:", len(negative_sequences))


# Create ML dataset
positive_data = pd.DataFrame({
    "sequence": positive_sequences,
    "label": 1
})

negative_data = pd.DataFrame({
    "sequence": negative_sequences,
    "label": 0
})

dataset = pd.concat(
    [positive_data, negative_data],
    ignore_index=True
)

print("\nTotal dataset:", len(dataset))
print(dataset.head())
print(dataset["label"].value_counts())
# Shuffle the dataset
dataset = dataset.sample(frac=1, random_state=42).reset_index(drop=True)

print("\nShuffled dataset:")
print(dataset.head())

# Save dataset
dataset_file = project_folder / "promoter_dataset.csv"
dataset.to_csv(dataset_file, index=False)

print("\nDataset saved at:")
print(dataset_file)

from sklearn.model_selection import train_test_split

# Split dataset into training and testing
X = dataset["sequence"]
y = dataset["label"]

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)

print("\nTraining sequences:", len(X_train))
print("Testing sequences:", len(X_test))

from sklearn.feature_extraction.text import CountVectorizer

# Convert DNA sequences into 3-mer features
vectorizer = CountVectorizer(
    analyzer="char",
    ngram_range=(3, 3)
)

X_train_kmer = vectorizer.fit_transform(X_train)
X_test_kmer = vectorizer.transform(X_test)

print("\nTraining feature shape:", X_train_kmer.shape)
print("Testing feature shape:", X_test_kmer.shape)

from sklearn.linear_model import LogisticRegression
# train the model
model = LogisticRegression(max_iter=1000)
model.fit(X_train_kmer,y_train)
print("\nModel training completed!")
# make prediction on test data
y_pred = model.predict(X_test_kmer)
print("\nPrediction compeleted!")
print("First 10 preictions:", y_pred[:10])

from sklearn.metrics import accuracy_score
accuracy = accuracy_score(y_test,y_pred)
print("\nModel Accuracy:", accuracy)
print("Accuracy percentage:", accuracy * 100,"%")

from sklearn.metrics import classification_report
print("\nclassification_report:")
print(classification_report(y_test,y_pred))
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

cm = confusion_matrix(y_test, y_pred)

disp = ConfusionMatrixDisplay(confusion_matrix=cm)
disp.plot()

plt.title("Promoter Classification Confusion Matrix")
plt.show()

from sklearn.metrics import roc_curve, roc_auc_score

# Predict probabilities
y_prob = model.predict_proba(X_test_kmer)[:, 1]

# ROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_prob)

# AUC
auc = roc_auc_score(y_test, y_prob)

print("AUC:", auc)

plt.figure()
plt.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve - Promoter Classification")
plt.legend()
plt.show()
import joblib
# save trained model
model_file = project_folder / "promoter_model.pkl"
joblib.dump(model,model_file)
print("\nModel save at:")
print(model_file)
# Test a new DNA sequence
new_sequence = input("\nEnter a DNA sequence: ").upper().strip()

# Convert sequence into the same 3-mer features
new_sequence_vector = vectorizer.transform([new_sequence])

# Predict
prediction = model.predict(new_sequence_vector)[0]

if prediction == 1:
    print("Prediction: PROMOTER")
else:
    print("Prediction: NOT A PROMOTER")
