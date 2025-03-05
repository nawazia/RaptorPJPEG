import argparse

from decode import decode
from encode2 import encode

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required=True, help="path to file")
    parser.add_argument("--decode", action="store_true", required=True, help="decode file")
    opt = parser.parse_args()
    print()
    print("Encoding file:", opt.path)
    opt.path = encode(opt)
    if opt.decode:
        print("Decoding file:", opt.path)
        decode(opt)
