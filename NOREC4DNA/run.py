import argparse
import numpy as np
import pandas as pd
from encode import encode
from decode import decode

ID_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"
DEFAULT_CHUNK_SIZE = 71

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", metavar="file", type=str, help="the bmp file to encode. Will automatically be converted to JPEG")
    parser.add_argument("--scans", metavar="file", type=str, default="/Users/i/Downloads/DNA_STorage/code/libjpeg/scans2.txt", help="the scans.txt file to use for JPEG encoding")
    parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=DEFAULT_CHUNK_SIZE,
                        help="size of chunks to split the file into")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", action="store_true", required=False, default=False)
    parser.add_argument("--p_thr", required=False, type=float, default=0.0, help="p_thr (default=0)")
    parser.add_argument("--save_number_of_chunks", metavar="save_number_of_chunks", required=False, type=bool,
                        default=False)
    parser.add_argument("--drop_upper_bound", metavar="drop_upper_bound", required=False, type=float, default=0.5,
                        help="upper bound for calculated error probability of packet before dropping")
    parser.add_argument("--overhead", metavar="overhead", required=False, type=float, default=0.4,
                        help="desired overhead of packets")
    parser.add_argument("--header_crc_str", metavar="header_crc_str", required=False, type=str, default="B")
    parser.add_argument("--xor_by_seed", action="store_true", required=False)
    parser.add_argument("--no_mask_id", action="store_true", required=False, help="do not mask ids")
    parser.add_argument("--id_spacing", required=False, type=int, default=0, help="spacing between ids (default=0)")
    parser.add_argument("--file", metavar="file", required=False, type=str, default="")
    
    all_solved = []
    # df = pd.read_csv("/Users/i/Downloads/DNA_STorage/code/NOREC4DNA/auc_all.csv")
    # df.drop(columns=['Unnamed: 0'], inplace=True)
    # df.head()
    # all_solved = df["0"].to_list()
    for i in range(len(all_solved), 8):
        row_solved = []
        for j in range(5):
            opt = parser.parse_args()
            opt.p_thr = i/10

            outfile, config_filename = encode(opt)
            opt.file = config_filename
            solved, bmpPath = decode(opt)
            overhead = (len(solved) - solved[-1])/solved[-1]
            row_solved.append(overhead)
            # area = np.trapz(solved, dx=1)
            # area /= len(solved)
            # row_solved.append(area)
            pd.DataFrame(row_solved).to_csv("/Users/i/Downloads/DNA_STorage/code/NOREC4DNA/auc_row.csv")
        all_solved.append([round(np.mean(row_solved),2), round(np.std(row_solved),2)])
        pd.DataFrame(all_solved).to_csv("/Users/i/Downloads/DNA_STorage/code/NOREC4DNA/auc_all.csv")
    print(all_solved)




