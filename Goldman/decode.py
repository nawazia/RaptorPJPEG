import argparse
import os
import random
import pickle
import tqdm

from Bio import SeqIO
    
def decimal(ternary_str):
    """
    Converts a ternary (base-3) number to a decimal (base-10) number.

    :param ternary_str: A string representing a ternary number (e.g., "10212")
    :return: The decimal equivalent as an integer.
    """
    decimal_value = 0
    power = 0

    # Process the ternary string from right to left
    for digit in reversed(ternary_str):
        if digit not in {'0', '1', '2'}:
            raise ValueError("Invalid ternary number. Only digits 0, 1, and 2 are allowed.")
        decimal_value += int(digit) * (3 ** power)
        power += 1

    return decimal_value

def trit_huffman_decode(encoded_data, codes):
    """Decodes the encoded trit Huffman data."""
    decoded = []
    rev_codes = {v: k for k, v in tqdm.tqdm(codes.items())}

    key = ""
    for trit in encoded_data:
        key += trit
        symbol = rev_codes.get(str(key), None)
        if symbol:
            decoded.append(rev_codes[key])
            key = ""

    return decoded

def DNAtoTrit(s4, prev='A'):
    # prev: {trit: next}
    hmap = {
        'A': {'C': '0', 'G': '1', 'T': '2'},
        'C': {'G': '0', 'T': '1', 'A': '2'},
        'G': {'T': '0', 'A': '1', 'C': '2'},
        'T': {'A': '0', 'C': '1', 'G': '2'}
    }

    # initial trit
    strand = hmap[prev][s4[0]]
    for idx, trit in enumerate(s4[1:]):
        strand += hmap[s4[idx]][trit]

    return strand

def OligosToDNA(fpath):
    O = []
    errors = 0
    max_idx = 0
    print("Parsing oligo FASTA file:")
    for seq_record in tqdm.tqdm(SeqIO.parse(fpath, "fasta")):
        # print('\n', seq_record.id, sep="")
        # print(len(seq_record))
        oligo = seq_record.seq

        # TODO: simulate PCR reverse (complementing?)
        if random.choice([0, 1]) == 1:
            oligo = oligo[::-1]

        pre, pro = oligo[0], oligo[-1]
        if pre in ['A', 'T']:
            assert pro in ['C', 'G'], pro
        elif pre in ['C', 'G']:
            assert pro in ['A', 'T'], pro
            oligo = oligo[::-1]
            pre, pro = oligo[0], oligo[-1]

        # remove 'orientation' nt
        oligo = oligo[1:-1]

        message = oligo[:100]
        IX = oligo[100:]

        # IX decoding
        IX = DNAtoTrit(IX, prev=message[-1])
        ID = IX[0:2]

        i3 = IX[2:14]

        P = IX[14:]
        if P != str((int(ID[1]) + int(i3[1]) + int(i3[3]) + int(i3[5]) + int(i3[7]) + int(i3[9]) + int(i3[11])) % 3):
            print('P error')
            errors += 1
            continue

        i = i3.lstrip('0')
        if len(i) == 0:
            i = '0'
        idx = decimal(i)
        max_idx = max(idx, max_idx)
        
        # message decoding
        if idx % 2 != 0:
            # odd
            message = message.reverse_complement()

        # print(len(message))
        O.append([idx, message])

    print("# of errors:", errors)

    print("Oligos to DNA:")
    # TODO: FUSE Os - this needs rework with sequencing sims.
    # create str of size max_idx
    s5 = "_" * (max_idx * 25 + 100)
    s5 = bytearray(s5, "ascii")
    for m_idx, message in tqdm.tqdm(O):
        for j_idx, j in enumerate(str(message)):
            s5[(m_idx * 25) + j_idx] = ord(j)
    s5 = s5.decode("ascii")
    return s5

def decode(opt):
    # load ternary Huffman codes
    print("Loading pre-computed ternary Huffman codes")
    codes_path = os.path.join(os.path.dirname(os.path.dirname(opt.path)), "codes", os.path.basename(opt.path).split(".")[0] + ".pickle")
    with open(codes_path, 'rb') as handle:
        codes = pickle.load(handle)

    s5 = OligosToDNA(opt.path)

    print("Decode single DNA-stream into trit-stream")
    s4 = DNAtoTrit(s5)

    print("Calculate ternary length, s2 & s3")
    s2 = s4[-20:]
    i = s2.lstrip('0')
    if len(i) == 0:
        i = '0'
    len_s1 = decimal(i)

    s3 = s4[len_s1:len(s4) - 20]
    assert len(s3) <= 24, len(s3)

    s1 = s4[:len_s1]
    print("Decoding codes -> symbols:")
    file = trit_huffman_decode(s1, codes)
    file = b''.join(file)

    outpath = os.path.join(os.path.dirname(os.path.dirname(opt.path)), os.path.basename(opt.path).split(".")[0] + "_dec.jpg")
    print("Writing to:", outpath)
    with open(outpath, "wb") as binary_file:
        # Write bytes to file
        binary_file.write(file)

    print("Done decoding!\n")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required=True)
    opt = parser.parse_args()
    decode(opt)
