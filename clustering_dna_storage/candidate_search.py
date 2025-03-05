import argparse
import json
import os
import glob
import random
import tqdm
from itertools import groupby
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Levenshtein import ratio

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

OLIGO_LENGTH = 200
DEBUG = False

START_ADAPTER = Seq("AATGTACTTCGTTCAGTTACGTATTGCT")
# END_ADAPTER = Seq("GCAATACGTAACTGAACGAAGT")


def filter_junk_reads(records, ids=None, similarity_threshold=0.85):
    """
    Removes all the sequences that are not similar to any others.
    Prints out percentage of sequences removed
    """
    filtered_records = []
    filtered_seqs = set()

    for record in tqdm(records):
        for record_ in records:
            
            # Checking if its the same record
            if record is record_:
                continue
            
            # Checking if the record is already in the filtered pool
            if not record_.seq in filtered_seqs:
                if ratio(record.seq, record_.seq) > similarity_threshold:
                    filtered_records.append(record_)
                    filtered_seqs.add(record_.seq)

    print(f"{100 - (len(filtered_records) * 100 / len(records))} percent sequences filtered out")

    return filtered_records
                    

def create_fasta_file(ids, strands, output_filepath='output.fasta'):
    with open(output_filepath, 'w') as f:
        for i, strand in enumerate(strands):
            f.write(f">{ids[i]}\n")
            f.write(strand + '\n')

def run_length_encode(data):
    """Returns run length encoded Tuples for string"""
    # A memory efficient (lazy) and pythonic solution using generators
    return ((x, sum(1 for _ in y)) for x, y in groupby(data))

def get_recovery_percentage(consensus_strand, original_strand):
    """Gets recovery percentage based on two strands. Chooses the length of the original strand to evaluate identity"""
    min_length = min(len(original_strand), len(consensus_strand))
    return sum([
        1 for i in range(min_length)
        if consensus_strand[i] == original_strand[i]]
        ) / len(original_strand)

def analyseBadread(opt):
    dict = {}
    for seq_record in SeqIO.parse(os.path.join(opt.path, "badread", "reads.fastq"), "fastq"):
        id = seq_record.description.split(",")[0].split(" ")[1]
        dict[id] = dict.get(id, 0) + 1
        
    dict = {k: v for k, v in sorted(dict.items(), key=lambda item: item[1])}
    # print([(k, v) for k, v in dict.items()])

    # drop random_seq & junk_seq
    dict.pop("random_seq", None)
    dict.pop("junk_seq", None)

    depths = list(dict.values())
    print("mean", sum(depths) / len(depths))
    print("std", sum([(x - sum(depths) / len(depths))**2 for x in depths]) / len(depths))
    print("min", min(depths))
    print("max", max(depths))

    plt.figure()
    plt.hist(depths, bins=100, label="depths")
    # plt.show()
    return pd.DataFrame(dict.items(), columns=['UUID', 'Depth'])

def main(opt):
    candidatesPath = glob.glob(os.path.join(opt.path, "candidates", "*.json"))[0]
    oligosPath = glob.glob(os.path.join(opt.path, "fasta", "*.fasta"))[0]

    if DEBUG:
        depth_df = analyseBadread(opt)

    with open(candidatesPath, 'r') as f:
        data = json.load(f)

    rec = []
    count = 0
    unrec = []
    score_df = pd.DataFrame()
    all = []
    candidates = set()

    for candidate in data["candidates"]:
        if len(candidate) == OLIGO_LENGTH:
            candidates.add(candidate)
        elif len(candidate) > OLIGO_LENGTH:
            # check for start adapter
            candidates.add(candidate[:OLIGO_LENGTH])

    # candidates = random.choices(list(candidates), k=1000)
    # print(f"{len(candidates)} / {len(data['candidates'])} candidates have length {OLIGO_LENGTH}")
    # create_fasta_file(["" for _ in range(len(candidates))], [candidate for candidate in candidates], os.path.join(opt.path, "rec", f"recovered_reads2.fasta"))
    # return
    
    print(f"{len(candidates)} / {len(data['candidates'])} candidates have length {OLIGO_LENGTH}")

    total = 0
    original_oligos = []
    for seq_record in SeqIO.parse(oligosPath, "fasta"):
        total += 1
        original_oligos.append(seq_record)
        

    with tqdm.tqdm(total=total) as pbar:
        for seq_record in original_oligos:
            score_max = 0
            oligo = seq_record.seq
            for candidate in candidates:
                if candidate == oligo:
                    rec.append(SeqRecord(
                        Seq(candidate),
                        id=seq_record.id,
                        description=seq_record.description))
                                        
                    count += 1
                    all.append(1)
                    break

                if DEBUG:
                    score_max = max(get_recovery_percentage(candidate, oligo), score_max)

            if DEBUG and score_max > 0:
                tempDF = pd.DataFrame([(seq_record.id, score_max)],columns=['UUID','Score'])
                score_df = pd.concat([score_df, tempDF])
                unrec.append(score_max)
                all.append(0)
            pbar.update(1)

    print(f"{count} / {total} original oligos recovered")
    if DEBUG:
        plt.figure()
        plt.hist(unrec, bins=100, label='unrecovered')
        print("mean", sum(unrec) / len(unrec))
        print("std", sum([(x - sum(unrec) / len(unrec))**2 for x in unrec]) / len(unrec))
        print("min", min(unrec))
        print("max", max(unrec))

        df = depth_df.merge(score_df, how='inner', on='UUID')
        df.dropna(inplace=True)
        print(df.describe())
        plt.figure()
        df.plot(x='Depth', y='Score', style='o')
        plt.show()

    # save recovered oligos to fasta with UUID
    create_fasta_file(["" for seq_record in rec], [str(seq_record.seq) for seq_record in rec], os.path.join(opt.path, "rec", f"recovered_reads.fasta"))

    with open(os.path.join(opt.path, "all", "all.json"), 'w') as f:
        json.dump(all, f)

    return

def match(opt):
    candidatesPath = glob.glob(os.path.join(opt.path, "candidates", "guesses_14_2_25_1.fasta"))[0]
    oligosPath = glob.glob(os.path.join(opt.path, "fasta", "*.fasta"))[0]

    total_can = 0
    count = 0
    score_df = pd.DataFrame()
    all = []

    for candidate in SeqIO.parse(candidatesPath, "fasta"):
        total_can += 1
        for oligo in SeqIO.parse(oligosPath, "fasta"):
            if candidate.seq == oligo.seq:
                all.append(oligo.seq)
                count += 1
        

    print(count)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to the dir")
    opt = parser.parse_args()

    match(opt)
