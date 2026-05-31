#!/usr/bin/env python3
import pandas as pd
import argparse

def calculate_roy_score(seq):
    # Roy & Chanfreau 2020: Downstream 19bp. 
    # First 6bp count as 2 points per 'A'. Remaining 13bp count as 1 point per 'A'.
    if pd.isna(seq) or len(seq) == 0:
        return 0
    
    # Ensure we only look at the first 19bp (or less if short)
    seq = seq[:19]
    
    score = 0
    # First 6bp (index 0 to 5)
    score += seq[:6].count('A') * 2
    # Remaining bp (index 6 to 18)
    score += seq[6:].count('A') * 1
    return score

def calculate_max_run_a(seq):
    if pd.isna(seq) or len(seq) == 0:
        return 0
    
    # Check only the immediate downstream region (20bp is standard SQANTI output)
    seq = seq[:20] 
    
    max_run = 0
    current_run = 0
    for char in seq:
        if char == 'A':
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 0
    return max_run

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="SQANTI3 classification.txt file")
    parser.add_argument("output_file", help="Name for the new output file")
    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep="\t")

    # 1. Calculate Roy Score
    print("Calculating Roy & Chanfreau Scores...")
    df['roy_score'] = df['seq_A_downstream_TTS'].apply(calculate_roy_score)

    # 2. Calculate Max Run of A's
    print("Calculating Max Run of A's...")
    df['max_run_A'] = df['seq_A_downstream_TTS'].apply(calculate_max_run_a)

    # Save
    df.to_csv(args.output_file, sep="\t", index=False)
    print(f"Done! Created {args.output_file} with new columns: 'roy_score' and 'max_run_A'")

if __name__ == "__main__":
    main()