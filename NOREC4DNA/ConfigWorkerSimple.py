import configparser
import typing

from norec4dna.Decoder import Decoder
from norec4dna.RU10Decoder import RU10Decoder
from norec4dna.ErrorCorrection import get_error_correction_decode, nocode

STATIC_NUM_CHUNKS = None  # 149
ID_LEN_FORMAT = "I"
NUMBER_OF_CHUNKS_LEN_FORMAT = "I"
PACKET_LEN_FORMAT = "I"
CRC_LEN_FORMAT = "L"
READ_ALL_BEFORE_DECODER = True

def decode2(file, error_correction=nocode, null_is_terminator=False,
               number_of_chunks=STATIC_NUM_CHUNKS, use_header_chunk=False, id_len_format=ID_LEN_FORMAT,
               number_of_chunks_len_format=NUMBER_OF_CHUNKS_LEN_FORMAT, packet_len_format=PACKET_LEN_FORMAT,
               crc_len_format=CRC_LEN_FORMAT, read_all=READ_ALL_BEFORE_DECODER,
               checksum_len_str=None, skip_solve=False, xor_by_seed=False,
               id_spacing=0, mask_id=True, p_thr=0.0):
        # print("Belief Propagation mode")
        x = RU10Decoder(file, use_headerchunk=use_header_chunk, error_correction=error_correction,
                        static_number_of_chunks=number_of_chunks, checksum_len_str=checksum_len_str,
                        xor_by_seed=xor_by_seed, mask_id=mask_id, id_spacing=id_spacing, p_thr=p_thr)
        x.read_all_before_decode = read_all
        x.decode(id_len_format=id_len_format,
                 number_of_chunks_len_format=number_of_chunks_len_format, packet_len_format=packet_len_format,
                 crc_len_format=crc_len_format)

        x.GEPP.insert_tmp()
        if not skip_solve:
            x.solve(partial=True)
        try:
            # solved_count = x.getSolvedCount()
            # if solved_count == 0:
            #     return False, 0, None, None

            x.filename = ""
            # x.saveDecodedFile(null_is_terminator=null_is_terminator, print_to_output=False,
            #                         return_file_name=True, partial_decoding=True)
        except:  # FileNotFoundError: #ValueError
            return False, x.getSolvedCount(), x.GEPP.result_mapping, x.filename
        
        if x.getSolvedCount() != (len(x.GEPP.result_mapping)):
            return False, x.getSolvedCount(), x.GEPP.result_mapping, x.filename
        else:
            return True, x.getSolvedCount(), x.GEPP.result_mapping, x.filename

class ConfigSimple:
    def __init__(self, filename):
        self.config: configparser.ConfigParser = configparser.ConfigParser()
        self.config.read(filename)
        assert len(self.config.sections()) == 1
        self.file = self.config.sections()[0]

    def execute(self, skip_solve=False) -> typing.Optional[typing.List[typing.Union[str, Decoder]]]:
        print("Decoding {}".format(self.file))
        decoder = self.__decode(self.file, self.config[self.file], skip_solve=skip_solve)
        return decoder

    def warn_unknown_items(self, config):
        known = ["error_correction", "repair_symbols", "as_mode_1_bmp", "number_of_splits", "split_index_position",
                 "split_index_length", "last_split_smaller", "is_null_terminated", "insert_header", "id_len_format",
                 "number_of_chunks_len_format", "packet_len_format", "crc_len_format", "algorithm", "number_of_chunks",
                 "read_all", "epsilon", "quality", "rules", "quality_len_format", "epsilon_len_format", "config_str",
                 "savenumberofchunks", "dropped_packets", "created_packets", "upper_bound", "asdna", "master_seed",
                 "mode_1_bmp", "chunk_size", "distribution", "checksum_len_str","xor_by_seed", "id_spacing", "mask_id",
                 "priority_chunks", "p_thr"]
        for cfg in config:
            if cfg not in known:
                print(f"[Warning] Config-entry '{cfg}' not known!")

    # @staticmethod
    def __decode(self, filename, decode_conf, skip_solve=False):
        self.warn_unknown_items(decode_conf)
        number_of_chunks = decode_conf.getint("number_of_chunks", None)
        e_correction = decode_conf.get("error_correction", "nocode")  # optional
        repair_symbols = decode_conf.getint("repair_symbols", 2)  # optional
        is_null_terminated = decode_conf.getboolean("is_null_terminated")  # optional, bool
        use_header_chunk = decode_conf.getboolean("insert_header")  # optional, bool
        id_len_format = decode_conf.get("id_len_format")  # optional, str
        number_of_chunks_len_format = decode_conf.get("number_of_chunks_len_format", "I")  # optional, str
        packet_len_format = decode_conf.get("packet_len_format", "I")  # optional, str
        crc_len_format = decode_conf.get("crc_len_format", "L")  # optional, str
        read_all_packets = decode_conf.getboolean("read_all", False)
        checksum_len_str = decode_conf.get("checksum_len_str", "")
        xor_by_seed = decode_conf.getboolean("xor_by_seed", False)
        id_spacing = decode_conf.getint("id_spacing", 0)
        mask_id = decode_conf.getboolean("mask_id", True)
        p_thr = decode_conf.getfloat("p_thr", 0.0)

        error_correction = get_error_correction_decode(e_correction, repair_symbols)
        print("File / Folder to decode: " + str(filename))
        try:
            decoded = decode2(filename, error_correction=error_correction, null_is_terminator=is_null_terminated,
                            id_len_format=id_len_format,
                            number_of_chunks_len_format=number_of_chunks_len_format,
                            packet_len_format=packet_len_format, crc_len_format=crc_len_format,
                            number_of_chunks=number_of_chunks,
                            use_header_chunk=use_header_chunk, read_all=read_all_packets,
                            checksum_len_str=checksum_len_str, skip_solve=skip_solve, xor_by_seed=xor_by_seed,
                            id_spacing=id_spacing, mask_id=mask_id, p_thr=p_thr)
        except Exception as ex:
            # raise ex
            print("Error decoding file: ", ex)
        return decoded
