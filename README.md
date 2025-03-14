# RaptorPJPEG

Official repo for *"Progressive decoding of DNA-stored JPEG data with on-the-fly error correction"*, 2025.

![Decoding cat GIF](data/cat.gif)

## Description

This repo combines a lot of packages and libraries into one full pipeline, to allow for full image encoding and decoding. We also use Badread as the storage channel simulator preceding clustering.

## Getting Started

### Dependencies

* CMake
* GCC v4.1+/Clang (Unix-like) or Microsoft Visual C++ 2005 or later (Windows)
* [Badread](https://github.com/rrwick/Badread)

### Installing

* First, [build libjpeg-turbo](https://github.com/libjpeg-turbo/libjpeg-turbo/blob/main/BUILDING.md):
```
git clone 'https://github.com/libjpeg-turbo/libjpeg-turbo.git'
mkdir {build_directory}
cd {build_directory}
cmake -G"Unix Makefiles" [additional CMake flags] {libjpeg-turbo}
make
```
* Modify cjpeg executable path in `NOREC4DNA/encode.py` to point to your own `{build_directory}/cjpeg` executable:
```
cjpeg_path = ".../libjpeg/cjpeg"
```
* Next, build JPEG decoder:
```
cd jpeg
mkdir build
cd build
cmake ..
make
```
* Modify jpeg executable path in `NOREC4DNA/encode.py` and `NOREC4DNA/decode.py` to point to your own `jpeg/build/jpeg` executable:
```
jpeg_path = ".../jpeg/build/jpeg"
```
* Next, create a conda env and [install NOREC4DNA](https://github.com/umr-ds/NOREC4DNA):
```
conda create -n "norec" python=3.10
conda activate norec
cd NOREC4DNA
pip install -r requirements.txt
python setup.py install
```
* Install [MUSCLE](https://github.com/rcedgar/muscle/releases) and change path in `clustering_dna_storage/strand_reconstruction.py`:
```
# Change the muscle path here
muscle_exe = ".../muscle-osx-arm64.v5.3"
```
## Usage
### Encoding and decoding

* Encoding an uncompressed bmp image:
```
python encode.py {file/to/image.bmp} [additional flags]
```
This generates a few files:
* image.jpg = standard JPEG
* image_FFDX.jpg = custom JPEG with base-8 RIs
* image_FFDX.jpg_RU10.fasta = FASTA containing generated oligo packets
* image_FFDX.jpg.ini = config file

* Decoding a FASTA:
```
python decode.py {file/to/image_FFDX.jpg.ini}
```
This creats a `tmp` dir, iteratively decodes & renders the image.
Images are saved as `tmp/patched_IDX.bmp`, and a GIF is generated at the end compiling all renders.

### Full pipeline

Raptor encoding an uncompressed image, cat.bmp:
```
python encode.py data/cat.bmp --chunk_size 47 --error_correction reedsolomon --insert_header --overhead 0.5 --p_thr 0.4
```

Add checksum bases to oligos for clustering, via `clustering_dna_storage/checksum.ipynb`:
```
# Read cat.fasta - add checksum encoding
# write fasta

RU10_path = "data/cat_FFDX.jpg_RU10.fasta"
checksum_path = os.path.join(os.path.dirname(RU10_path), "cat_FFDX_checksum.fasta")
recovered_path = os.path.join(os.path.dirname(RU10_path), "cat_FFDX_rec.fasta")
oligoLen = 204

checksum = CheckSum4(reference_length=oligoLen)
original_strands, original_strand_ids = read_synthesized_strands_from_file(file_path=RU10_path)
encoded_strands = checksum.encode(original_strands)
create_fasta_file(ids=original_strand_ids, strands=encoded_strands, output_filepath=checksum_path)
print("Checksum encoded strands saved to", checksum_path)
```

Running `data/cat_FFDX_checksum.fasta` thru [Badread](https://github.com/rrwick/Badread):
```
badread simulate --reference data/cat_FFDX_checksum.fasta --quantity 60x --identity 97,99,1.0 | gzip > reads.fastq.gz
```
Run noisy reads (unzipped at `data/reads.fastq`) thru agglomerative k-mer clustering, via `clustering_dna_storage/checksum.ipynb`:
```
# Raptor
records = get_fastq_records(fastq_filepath="data/reads.fastq")
reads = [str(i.seq) for i in records]
ids = [get_badread_strand_id(i) for i in records]

clustering = Clustering(strand_pool=reads, reference_length=oligoLen+4, original_strands=encoded_strands, strand_pool_ids=ids, distance_threshold=40)
clustering.run_pipeline(eval=True)

# We want to recover the original strands after checksum
decoded_strands = checksum.decode(candidates=clustering.candidates, n_reference_strands=len(original_strands), clustered_seqs=clustering.clustered_seqs, n_guesses=5, guesses=True)
ids = ["" for i in range(len(decoded_strands))]
create_fasta_file(ids, decoded_strands, output_filepath=recovered_path)
```
Running Raptor decoding on Badread + clustered strands:
```
python decode.py data/cat_FFDX.jpg.ini --badread data/cat_FFDX_rec.fasta
```

## Authors

[Ibrahim Nawaz](mailto:ibrahim.nawaz22@imperial.ac.uk)

## Version History

* 0.2
    * Updated README, LICENSE and priority encoding
* 0.1
    * Initial Release

## License

This project is licensed under the [MIT License](LICENSE)

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
