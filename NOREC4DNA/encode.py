import argparse
import os
import subprocess
import math

from PIL import Image

from norec4dna.Encoder import Encoder
from norec4dna.RU10Encoder import RU10Encoder
from norec4dna.rules.FastDNARules import FastDNARules
from norec4dna.ErrorCorrection import nocode, get_error_correction_encode
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.helper.helper import readMarkers

ID_LEN_FORMAT = "H"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"
DEFAULT_CHUNK_SIZE = 71
cjpeg_path = "/Users/i/Downloads/ICL/DNA_STorage/code/libjpeg/cjpeg"
jpeg_path = "/Users/i/Downloads/ICL/DNA_STorage/code/jpeg/cmake-build-release/jpeg"

def average_base8_length(X):
    """
    Calculates the average length of numbers in the range 0 to X in base-8.

    Args:
        X (int): The upper limit of the range.

    Returns:
        float: The average length of the numbers in base-8.
    """
    total_digits = 0
    for n in range(X):
        if n == 0:
            total_digits += 1  # 0 requires 1 digit
        else:
            total_digits += math.floor(math.log(n, 8)) + 1

    return total_digits / (X)

def findRI(opt):
    """
    Find optimal restart interval for the given file, based on average Block (MCU) bit size in DC scan.

    Parameters
    ----------
    opt : Namespace
        args

    Returns
    -------
    None
    """
    print("Finding optimal restart interval")
    try:
        # compress JPEG
        # note: -sample 1x1 is used to avoid chroma subsampling
        # scans.txt is a file containing the number of MCUs in each scan. Defaults to 1 DC scan & 3 AC scans, but can be changed.
        command = [cjpeg_path, "-progressive", "-sample", "1x1", "-scans", opt.scans, opt.path]
        output_file = opt.path.split(".")[0] + ".jpg"
        print(os.path.exists(output_file), output_file)

        with open(output_file, "wb") as outfile:
            subprocess.run(command, stdout=outfile, check=True, shell=False)
    except subprocess.CalledProcessError as e:
        print(f"Error compressing JPEG: {e}")
        return False
    except FileNotFoundError as e:
        print(f"cjpeg or scan file not found. Make sure it is installed and the scan file exists: {e}")
        return False
    except Exception as generic_exception:
        print(f"An unexpected error occured: {generic_exception}")
        return False

    # first, run JPEG decoder without mapping.bin. Returns number of blocks in DC scan
    try:
        command = [jpeg_path, output_file]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except subprocess.CalledProcessError as e: 
        print("Error occurred while running the C program:")
        raise e
    
    # get average block size, in bytes
    avgBlockSize = result.returncode / 100
    print(f"Average block size: {avgBlockSize * 8} bits")

    # find number of blocks in DC scan
    with Image.open(output_file) as img:
        width, height = img.size

    numBlocksWidth = math.ceil(width / 8)
    numBlocksHeight = math.ceil(height / 8)
    numBlocks = numBlocksWidth * numBlocksHeight
    print(f"Number of blocks: {numBlocks}")

    # find optimal restart interval. Break RI when chunk > opt.chunk_size*0.66
    ri = None
    for i in range(1, 100, 2):
        numRI = math.floor(numBlocks / i)

        chunk = math.ceil(avgBlockSize * i) + 2*(average_base8_length(numRI))
        if chunk > 0.66*opt.chunk_size:
            print(f"Found optimal restart interval: {i}")
            ri = i
            break

    if ri is None:
        print("No optimal restart interval found.")
        raise ValueError("No optimal restart interval found.")

    print("Encoding into JPEG with optimal restart")
    try:
        # compress JPEG
        # note: -sample 1x1 is used to avoid chroma subsampling
        # scans.txt is a file containing the number of MCUs in each scan. Defaults to 1 DC scan & 3 AC scans, but can be changed.
        command = [cjpeg_path, "-progressive", "-restart", f"{ri}B", "-sample", "1x1", "-scans", opt.scans, opt.path]
        output_file = opt.path.split(".")[0] + ".jpg"

        with open(output_file, "wb") as outfile:
            subprocess.run(command, stdout=outfile, check=True, shell=False)

        opt.path = output_file
        print("Saved at:", output_file, "\n\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error compressing JPEG: {e}")
        return False
    except FileNotFoundError as e:
        print(f"cjpeg or scan file not found. Make sure it is installed and the scan file exists: {e}")
        return False
    except Exception as generic_exception:
        print(f"An unexpected error occured: {generic_exception}")
        return False


def _encode(file, chunk_size=DEFAULT_CHUNK_SIZE, error_correction=nocode, insert_header=False,
               save_number_of_chunks_in_packet=False, mode_1_bmp=False, prepend="", append="", upper_bound=0.5,
               overhead=0.40, checksum_len_str=None, xor_by_seed=False, mask_id=True, id_spacing=0, priority_chunks=None, p_thr=0.0):
        number_of_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(file, chunk_size)
        dist = RaptorDistribution(number_of_chunks, p_thr=p_thr)
        rules = FastDNARules()
        x = RU10Encoder(file, number_of_chunks, dist, insert_header=insert_header, pseudo_decoder=None,
                        chunk_size=0, rules=rules, error_correction=error_correction,
                        packet_len_format=PACKET_LEN_FORMAT,
                        crc_len_format=CRC_LEN_FORMAT, number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT,
                        id_len_format=ID_LEN_FORMAT, save_number_of_chunks_in_packet=save_number_of_chunks_in_packet,
                        mode_1_bmp=mode_1_bmp, prepend=prepend, append=append, drop_upper_bound=upper_bound,
                        checksum_len_str=checksum_len_str, xor_by_seed=xor_by_seed, mask_id=mask_id,
                        id_spacing=id_spacing)
        x.set_overhead_limit(overhead)
        if priority_chunks:
            x.set_priority_chunks(priority_chunks)
        x.encode_to_packets()
        outfile = x.save_packets_fasta(file_ending="_RU10", seed_is_filename=True)
        return x, outfile

def baseEight(n):
    if n == 0:
        return '0'
    nums = []
    while n:
        n, r = divmod(n, 8)
        nums.append(str(r))
    return ''.join(reversed(nums))

def addFFDX(opt):
    """
    Add global RI markers so allow for recovery of fragmentated data. Restart markers are converted into global markers with base-8.

    Parameters
    ----------
    opt : Namespace
        args

    Returns
    -------
    outpath : str
        Path to new FFDX.jpg
    """
    print("Adding FFDX markers")
    with open(opt.path, "rb") as f:
       markers, file = readMarkers(f)

    filesize = len(file)
    print("markers before:")
    print(markers ,"\n")

    assert markers[0][1] == b'\xd8' and markers[-1][1] == b'\xd9'

    # find first FFDA - start of scan
    firstDA = [marker for _, marker in markers].index(b'\xda')
    count = 0
    maxLen = 0
    offset = 0
    for index, marker in markers[firstDA:-1]:
        index = index + offset
        if marker in [b'\xc4']:
           continue
        # new scan, reset count
        if marker == b'\xda':
            count = 0
            continue

        rii = list(baseEight(count))
        assert eval("b'\\xd" + rii[-1] + "'") == marker, f"{rii[-1]} != {marker}"
        for ri in rii[-2::-1]:
            byte = eval("b'\\xd" + ri + "'")
            file.insert(index, byte)
            file.insert(index, b'\xff')
            # add 2 to all marker indices
            offset += 2

        count += 1
        maxLen = max(len(rii), maxLen)
    
    print("max len of new RIs",maxLen)
    outpath = os.path.join(os.path.dirname(opt.path), os.path.basename(opt.path).split(".")[0] + "_FFDX.jpg")
    print("Saving at:", outpath)
    increase = round((len(file) / filesize) * 100, 1)
    print(f"% of file size increase: {increase}%")
    with open(outpath, "wb") as f:
       f.write(b"".join(file))

    with open(outpath, "rb") as f:
       newMarkers, _ = readMarkers(f, dropFFDX=True)
       print(newMarkers)
    
    # find end of first scan
    firstDA = [marker for _, marker in newMarkers].index(b'\xda')
    eosIdx = newMarkers[firstDA+1]
    # ceiling division
    eosChunkIdx = -(eosIdx[0] // -opt.chunk_size)

    return outpath, eosChunkIdx

def encodeFile(opt, _file, eosChunkIdx):
    print("File to encode: " + str(_file))
    error_correction = get_error_correction_encode(opt.error_correction, opt.repair_symbols) # default from demo_raptor_encode.py
    encoder_instance, outfile = _encode(_file, chunk_size=opt.chunk_size, error_correction=error_correction,
                                    mode_1_bmp=False, insert_header=opt.insert_header,
                                    append="",
                                    save_number_of_chunks_in_packet=opt.save_number_of_chunks, upper_bound=opt.drop_upper_bound,
                                    overhead=opt.overhead, checksum_len_str=opt.header_crc_str, xor_by_seed=opt.xor_by_seed, mask_id=not opt.no_mask_id,
                                    id_spacing=opt.id_spacing, priority_chunks=eosChunkIdx, p_thr=opt.p_thr)
    conf = {'error_correction': opt.error_correction, 'repair_symbols': opt.repair_symbols, 'asdna': True,
            'number_of_splits': 1, 'read_all': True, "priority_chunks": eosChunkIdx, "p_thr": opt.p_thr, "jpg": _file}
    config_filename = encoder_instance.save_config_file(conf, add_dot_fasta=True)
    print("Saved config file: %s" % config_filename)
    return outfile, config_filename, encoder_instance.oligoLen

def encode(opt):
    if(not findRI(opt)):
        print("Error finding optimal restart interval.")
        exit(1)

    # add global base-8 RI markers
    file, eosChunkIdx = addFFDX(opt)
    
    # # encode FFDX file
    outfile, config_filename, oligoLen = encodeFile(opt, file, eosChunkIdx)
    print("oligoLen:", oligoLen)
    return outfile, config_filename


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", metavar="file", type=str, help="the bmp file to encode. Will automatically be converted to JPEG")
    parser.add_argument("--scans", metavar="file", type=str, default="/Users/i/Downloads/ICL/DNA_STorage/code/libjpeg/scans2.txt", help="the scans.txt file to use for JPEG encoding")
    parser.add_argument("--chunk_size", metavar="chunk_size", required=False, type=int, default=DEFAULT_CHUNK_SIZE,
                        help="size of chunks to split the file into")
    parser.add_argument("--error_correction", metavar="error_correction", required=False, type=str, default="nocode",
                        help="Error Correction Method to use; possible values: \
                        nocode, crc, reedsolomon, dna_reedsolomon (default=nocode)")
    parser.add_argument("--repair_symbols", metavar="repair_symbols", required=False, type=int, default=2,
                        help="number of repair symbols for ReedSolomon (default=2)")
    parser.add_argument("--insert_header", action="store_true", required=False, default=False)
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
    parser.add_argument("--p_thr", required=False, type=float, default=0.0, help="p_thr (default=0)")

    opt = parser.parse_args()

    encode(opt)
