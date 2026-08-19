# Human Promoter Classifier

## Project Overview

Human Promoter Classifier is a bioinformatics project designed to analyze DNA sequences and identify potential human promoter sequences using computational approaches.

## Objective

The main objective of this project is to process promoter and non-promoter DNA sequences, examine their sequence characteristics, and prepare the data for computational classification.

## Dataset

The project uses human promoter sequence data along with negative/non-promoter sequences.

The datasets include:

- Human promoter sequences
- Negative DNA sequences
- Promoter dataset in CSV format
- FASTA-formatted DNA sequences

## Project Structure

human-promoter-classifier/

├── data/

│   ├── promoter_dataset.csv

│   ├── negative_sequence.fasta

│   └── HSEPD News 006HG38.txt

├── src/

│   ├── check_fasta.py

│   ├── check_promoter.py

│   └── negative_datasets.py

└── README.md

## Methodology

The project workflow includes:

1. Collecting human promoter sequence data.
2. Preparing positive and negative DNA sequence datasets.
3. Checking FASTA sequence format and sequence quality.
4. Processing promoter sequences using Python.
5. Preparing the dataset for computational analysis and classification.

## Technologies Used

- Python
- Bioinformatics
- FASTA
- CSV
- DNA sequence analysis

## Future Improvements

Future work can include feature extraction, machine-learning model training, performance evaluation, and visualization of classification results.

## Author

Molecular Biology Student | Bioinformatics & Computational Biology
