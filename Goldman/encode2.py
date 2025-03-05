import os
import queue
import argparse
import random
import pickle
import tqdm
from uuid import uuid4

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


class HuffmanTreeNode():
    def __init__(self, char, freq):
        self.data = char

        self.freq = freq

        self.left = None
        self.middle = None
        self.right = None

    def __lt__(self, other):
        return self.freq < other.freq
    
def generateTree(pq):
    while pq.qsize() != 1:
        left = pq.get()
        middle = pq.get()
        right = pq.get()

        node = HuffmanTreeNode("$", left.freq + middle.freq + right.freq)
        node.left = left
        node.middle = middle
        node.right = right

        pq.put(node)

    return pq.get()

def generateCodes(root: HuffmanTreeNode, arr, top, codes):

    if root.left:
        arr[top] = 0
        generateCodes(root.left, arr, top+1, codes)
    if root.middle:
        arr[top] = 1
        generateCodes(root.middle, arr, top+1, codes)
    if root.right:
        arr[top] = 2
        generateCodes(root.right, arr, top+1, codes)
    
    # base case
    if not root.left and not root.middle and not root.right:
        code = "".join(str(arr[i]) for i in range(top))

        codes[root.data] = code

def printCodes(root: HuffmanTreeNode, arr, top):
    if root.left:
        arr[top] = 0
        printCodes(root.left, arr, top+1)
    if root.middle:
        arr[top] = 1
        printCodes(root.middle, arr, top+1)
    if root.right:
        arr[top] = 2
        printCodes(root.right, arr, top+1)
    
    # base case
    if not root.left and not root.middle and not root.right:
        print(root.data, end=" ")
        for i in range(top):
            print(arr[i], end="")
        print()

def HuffmanCodes(data):
    pq = queue.PriorityQueue()

    for k, v in data.items():
        newNode = HuffmanTreeNode(k, v)
        pq.put(newNode)

    root = generateTree(pq)

    arr = [0] * 100
    top = 0

    # printCodes(root, arr, top)

    codes = {}
    generateCodes(root, arr, top, codes)
    return codes

markerIdx = {
    "SOI": 0,
    "APPN": 1,
    "DQT": 2,
    "SOF": 3,
    "DHT": 4,
    "SOS": 5,
    "EOI": 6,
    "[EXTRA]": 7
}

def readMarker(current):
    if current == b'\x00':
        return

    hmap = {
        b'\xD8': "SOI",
        b'\xD9': "EOI",
        b'\xDB': "DQT",
        b'\xC0': "SOF",
        b'\xC2': "SOF",
        b'\xC4': "DHT",
        b'\xDA': "SOS",
    }

    # handle APPN markers
    if b'\xE0' <= current <= b'\xEF':
        marker = "APPN"
    else:
        marker = hmap[current]

    return marker

def readFile(fpath, data):
    file = []
    markers = [None]

    prev = False
    with open(fpath, "rb") as f:
        while (byte := f.read(1)):
            marker = None
            # check marker and append
            if prev:
                # return marker string, None if real-FF
                marker = readMarker(byte)
                prev = False

            if byte == b'\xFF':
                prev = True
            
            # increment frequency
            data[byte] = data.get(byte, 0) + 1
            # append to file
            file.append(byte)

            if marker:
                markers.append(marker)
            else:
                markers.append(markers[-1])
    
    markers = markers[2:]
    markers.append("EOI")

    assert len(file) == len(markers), f"{len(file)} != {len(markers)}"

    return file, markers

def trit_huffman_encode(data, codes):
    """Encodes the input data using the Trit Huffman codes."""
    return "".join(codes[char] for char in tqdm.tqdm(data))

def trit_huffman_encode2(data, codes, markers):
    """Encodes the input data using the Trit Huffman codes."""
    enc_code = ""
    enc_marker = ""

    for idx, char in enumerate(tqdm.tqdm(data)):
        marker = markers[idx]
        code = codes[char]
        # print(code)
        enc_code += code
        enc_marker += str(markerIdx[marker]) * len(code)

    assert len(enc_code) == len(enc_marker)
    return enc_code, enc_marker


def trit_huffman_decode(encoded_data, codes):
    """Decodes the encoded trit Huffman data."""
    decoded = []
    rev_codes = {v: k for k, v in codes.items()}

    key = ""
    for trit in encoded_data:
        key += trit
        symbol = rev_codes.get(str(key), None)
        if symbol:
            decoded.append(rev_codes[key])
            key = ""

    return decoded

def ternary(n):
    if n == 0:
        return '0'
    nums = []
    while n:
        n, r = divmod(n, 3)
        nums.append(str(r))
    return ''.join(reversed(nums))

def tritToDNA(s4, prev='A'):
    # prev: {trit: next}
    hmap = {
        'A': {'0': 'C', '1': 'G', '2': 'T'},
        'C': {'0': 'G', '1': 'T', '2': 'A'},
        'G': {'0': 'T', '1': 'A', '2': 'C'},
        'T': {'0': 'A', '1': 'C', '2': 'G'}
    }

    # initial trit
    strand = hmap[prev][s4[0]]
    for trit in s4[1:]:
        strand += hmap[strand[-1]][trit]

    return strand

def DNAtoOligos(s5, N, ID, markers):
    F = []

    for idx, pos in enumerate(tqdm.tqdm(range(0, N-75, 25))):
        marker = markers[pos:pos+100]
        message = Seq(s5[pos:pos+100])
        # print(marker, message)
        assert len(message) == 100, len(message)
        if idx % 2 != 0:
            # odd
            message = message.reverse_complement()
        i3 = '0' * (12 - len(ternary(idx))) + ternary(idx)
        P = (int(ID[1]) + int(i3[1]) + int(i3[3]) + int(i3[5]) + int(i3[7]) + int(i3[9]) + int(i3[11])) % 3

        IX = ID + i3 + str(P) # 2 + 12 + 1
        IX = tritToDNA(IX, prev=message[-1])

        oligo = message + IX

        # add 'orientation' nt
        if oligo[0] == 'A':
            oligo = 'T' + oligo
        elif oligo[0] == 'T':
            oligo = 'A' + oligo
        else:
            oligo = random.choice(['A', 'T']) + oligo

        if oligo[-1] == 'C':
            oligo = oligo + 'G'
        elif oligo[-1] == 'G':
            oligo = oligo + 'C'
        else:
            oligo = oligo + random.choice(['C', 'G'])
        
        # print(oligo)
        # print("GC content:", gc_fraction(oligo))
        for _ in range(1000):
            F.append(SeqRecord(
                oligo,
                id=f"{uuid4()}",                           # id="gi|14150838|gb|AAK54648.1|AF376133_1",
                description=f"{marker[0]}_{pos}_{pos + 100}",                  # description="chalcone synthase [Cucumis sativus]",
            ))

    return F

        

def encode(opt):
    data = {}       # used for Huffman encoding
    file, markers = readFile(opt.path, data)
    if len(data) % 2 == 0:
        print("Adding extra byte to make even")
        data["[EXTRA]"] = 0

    print("Generating ternary Huffman codes")
    codes = HuffmanCodes(data)
    codes_path = os.path.join(os.path.dirname(opt.path), "codes", os.path.basename(opt.path).split(".")[0] + ".pickle")
    with open(codes_path, 'wb') as handle:
        pickle.dump(codes, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print("Encoding symbols -> codes:")
    # s1 = trit_huffman_encode(file, codes)
    s1, s1_markers = trit_huffman_encode2(file, codes, markers)

    print("Adding ternary length, s2 & s3")
    tern = ternary(len(s1))
    s2 = '0' * (20 - len(tern)) + tern
    s3 = (25 - ((len(s1) + len(s2)) % 25)) * '0' # always <=24
    s4 = s1 + s3 + s2   # X*25 = len(data) + (<=24) + 20
    assert len(s4) % 25 == 0

    s4_markers = s1_markers + ('7' * (len(s3+s2)))
    assert len(s4) == len(s4_markers)

    # encode trit-stream into DNA
    print("Encode trit-stream into single DNA-stream")
    s5 = tritToDNA(s4)
    N = len(s5)

    # we're just encoding 1 file, so ID=00
    print("DNA to oligos:")
    ID = '00'
    F = DNAtoOligos(s5, N, ID, s4_markers)
    
    print("Writing oligo FASTA file")
    outpath = os.path.join(os.path.dirname(opt.path), "fasta", os.path.basename(opt.path).split(".")[0] + ".fasta")
    SeqIO.write(F, outpath, "fasta")
    print("Done encoding!\n\n")
    return outpath


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required=True)
    opt = parser.parse_args()
    encode(opt)
