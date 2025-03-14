import argparse
import os
import re
import random
import subprocess
import numpy as np
import pandas as pd
import plotly.express as px
import shutil
import shlex
import imageio
from PIL import Image, ImageDraw, ImageFont

import ConfigWorkerSimple
from norec4dna import helper

from Bio import SeqIO

DEBUG = False
jpeg_path = "/Users/i/Downloads/ICL/DNA_STorage/code/jpeg/cmake-build-release/jpeg"

def create_random_strand(strand_length):
    base_choices = ['A', 'C', 'T', 'G']
    return "".join([random.choice(base_choices) for _ in range(strand_length)])

def create_tmp(src_filename, target_filename, badread=None):
    shutil.rmtree(os.path.dirname(target_filename))
    os.makedirs(os.path.dirname(target_filename))
    f = open(src_filename)
    first_line, remainder = f.readline(), f.read()
    f.close()

    t = open(target_filename,"w")
    if badread:
        replacement_line = "[" + os.path.join(os.path.dirname(first_line[1:-2]), "tmp", os.path.basename(badread)) + "]"
    else:
        replacement_line = "[" + os.path.join(os.path.dirname(first_line[1:-2]), "tmp", os.path.basename(first_line[1:-2])) + "]"
    t.write(replacement_line + "\n")
    t.write(remainder)
    t.close()

    return first_line[1:-2], replacement_line[1:-1]

def create_fasta_file(ids, strands, output_filepath='output.fasta'):
    with open(output_filepath, 'w') as f:
        for i, strand in enumerate(strands):
            f.write(f">{ids[i]}\n")
            f.write(strand + '\n')

def decode(opt):
    """
    Parameters
    ----------
    opt : Namespace
        args

    Returns
    -------
    solved : list
        array of numSolved
    bmpPath : str
        path to JPEG decoded bmps.
    """
    tmp_cpath = os.path.join(os.path.dirname(opt.file), "tmp", os.path.basename(opt.file))
    fpath, tmp_fpath = create_tmp(opt.file, tmp_cpath, badread=opt.badread)

    original_file = [str(sequence.seq) for sequence in SeqIO.parse(fpath, "fasta")]
    original_config = ConfigWorkerSimple.ConfigSimple(opt.file)
    catPath = original_config.config[original_config.config.sections()[0]].get("jpg")

    oligoLen = len(original_file[0])
    print("oligoLen:", oligoLen)

    solved = [0]
    totalChunks = None
    chunks = []

    with open(catPath, "rb") as f:
        markers, cat = helper.readMarkers(f)
        print(markers)
    
    cat = b"".join(cat)

    original_file = np.array(original_file)
    for numStrands in range(1, len(original_file)):
        strands_idx = np.random.choice(len(original_file), numStrands, replace=False)
        strands = original_file[strands_idx]

        create_fasta_file(["" for _ in strands], strands, tmp_fpath)

        x = ConfigWorkerSimple.ConfigSimple(tmp_cpath)
        chunk_size = x.config[x.config.sections()[0]].getint("chunk_size")
        if not totalChunks:
            totalChunks = x.config[x.config.sections()[0]].getint("number_of_chunks")
        try:
            (res, numSolved, mapping, filename) = x.execute()
        except:
            continue

        if isinstance(mapping, np.ndarray):
            assert len(mapping) == totalChunks
            chunks.append(mapping.reshape(-1))

        if filename and (solved[-1] < numSolved):
            print("Running JPEG decoding")
            outPath = os.path.join(os.path.dirname(tmp_fpath), f"patched_{numStrands:04}.jpg")
            mappingPath = os.path.join(os.path.dirname(tmp_fpath), f"{totalChunks}_{chunk_size}.bin")
            mapping.reshape(-1).astype(np.int32).tofile(mappingPath)

            with open(filename, 'rb') as in_file:
                with open(outPath, 'wb') as out_file:
                    inFile = in_file.read()
                    if mapping[0][0] == -1:
                        inFile = inFile[46:]

                    inFile = bytearray(inFile)
                    # add back the 'header' markers
                    for i in markers:
                        if i[1] in [ b'\xd0', b'\xd1', b'\xd2', b'\xd3', b'\xd4', b'\xd5', b'\xd6', b'\xd7']:
                            continue

                        if i[1] in [b'\xd8', b'\xd9']:
                            # no length
                            inFile[i[0]:i[0]+2] = b'\xff' + i[1]
                            if i[1] == b'\xd9':
                                inFile = inFile[:i[0]+2]
                            continue
                        
                        length = int.from_bytes(cat[i[0]+2:i[0]+4], "big")
                        inFile[i[0]:i[0]+length+2] = b'\xff' + i[1] + cat[i[0]+2:i[0]+length+2]
                    
                    out_file.write(inFile)
                    
            print("Saved patched.jpg at:", outPath)
            print("Saved mapping.bin at:", mappingPath)

            try:
                command = [jpeg_path, outPath, mappingPath]
                print(f"Running command: {shlex.join(command)}")  # Use shlex.join for quoting
                result = subprocess.run(command, capture_output=True, text=True, check=True) 
                print(result.stdout)
                print(result.stderr)
            except subprocess.CalledProcessError as e:
                print("Error occurred while running the C program:") 
                print(e.stdout)
                print(e.stderr)
            os.remove(outPath)
        solved.append(numSolved)

        print(f"numStrands: {numStrands} / {len(original_file)}\t\tnumber of solved: {numSolved} / {totalChunks}")
        if res != False:
            break
    
    assert solved[-1] == totalChunks

    if DEBUG:
        df = pd.DataFrame(np.array(chunks))
        df.to_csv("/Users/i/Downloads/ICL/DNA_STorage/code/NOREC4DNA/all.csv")
    return solved, os.path.dirname(tmp_fpath)

def create_gif_with_filenames(image_folder, output_gif, duration=0.1, font_path=None, font_size=20, font_color=(255, 255, 255)): #Added font parameters
    """Creates a GIF from a folder of images with filenames overlaid.

    Parameters
    -------
    image_folder : str
        Path to the folder containing the images.
    output_gif : str
        Path to the output GIF file.
    duration : float
        Duration of each frame in seconds (default: 0.1 seconds). Defaults to 0.1
    font_path : str
        Path to the font file. Defaults to None. If None, a default font is used.
    font_size : int
        Font size for the filename text. Defaults to 20.
    font_color : tuple
        Color of the font. Defaults to (255, 255, 255), white.
    """
    try:
        image_files = sorted([
            f for f in os.listdir(image_folder)
            if re.match(r"patched_.*\.bmp", f)  # Adjust regex if needed
        ], key=lambda x: int(x[8:-4]))  # Numerical sort (important!)

        if not image_files:
            raise ValueError(f"No image files found in folder: {image_folder}")

        images = []
        for filename in image_files:
            img_path = os.path.join(image_folder, filename)
            img = imageio.imread(img_path)

            # Convert to PIL Image (if needed)
            pil_img = Image.fromarray(img)

            # Create a drawing context
            draw = ImageDraw.Draw(pil_img)

            # Choose font
            if font_path:
                font = ImageFont.truetype(font_path, font_size)
            else:
                font = ImageFont.load_default() # Use default font if none provided

            # Add filename text
            text_position = (10, 10)  # Adjust position as needed
            draw.text(text_position, filename, font=font, fill=font_color)

            # Convert back to NumPy array (for imageio)
            img_with_text = np.array(pil_img)
            images.append(img_with_text)

        imageio.mimsave(output_gif, images, fps=1 / duration)  # Save as GIF

        print(f"GIF created successfully: {output_gif}")

    except (ValueError, FileNotFoundError, OSError, Exception) as e:
        print(f"Error creating GIF: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", metavar="file", type=str, help="config file")
    parser.add_argument("--badread", type=str, default=None, help="path to recovered (badread'd & clustered) file path")
    opt = parser.parse_args()
    solved, bmpPath = decode(opt)

    create_gif_with_filenames(bmpPath, os.path.join(os.path.dirname(bmpPath), "cat.gif"), duration=1.0)

    if DEBUG:
        fig = px.line(x=range(len(solved)), y=solved)
        fig.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
        xaxis=dict(title=dict(text="Number of oligos received")),
        yaxis=dict(title=dict(text="Number of chunks recovered")),
        plot_bgcolor='white',
        font_family="Arial",
        font_color="black",
        font_size=16,)

        fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
        )
        fig.update_yaxes(
            mirror=True,
            ticks='outside',
            showline=True,
            linecolor='black',
            gridcolor='lightgrey',
        )
        import plotly.io as pio
        pio.write_image(fig, "/Users/i/Desktop/solved.png", scale=3, width=1080, height=620)
        # fig.write_html(os.path.join(os.path.dirname(bmpPath), "solved.html"))
        fig.show()
        pd.DataFrame(solved).to_csv("/Users/i/Downloads/ICL/DNA_STorage/code/NOREC4DNA/solved.csv")
