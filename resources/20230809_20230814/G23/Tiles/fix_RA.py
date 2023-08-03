import numpy as np
import pandas as pd
from pathlib import Path
import shutil
from tqdm import tqdm

all_tile_files = list(Path(".").glob("G23_T*.fld"))
all_guide_files = list(Path(".").glob("G23_GT*.fld"))


def copy_to_backup(fname):
    backup_fname = f"{fname}.backup"
    shutil.copy(fname, backup_fname)


def fix_header(fname):
    with open(fname, "r") as f:
        contents = f.readlines()
    header = contents[:3]
    coordinates_line = header[1].strip().split(" ")
    ra = coordinates_line[1]
    dec = coordinates_line[2]

    new_ra = float(ra)
    if float(ra) < 0:
        new_ra += 360

    new_coordinates_line = f"# {new_ra} {dec}\n"
    new_header = [header[0], new_coordinates_line, header[2]]

    return new_header


def write_correct_tile_file(filename, df, header):
    with open(filename, "w") as g:
        for line in header:
            g.write(line)
        df.to_csv(g, index=False)


def write_correct_guide_file(filename, df, header):
    with open(filename, "w") as g:
        for line in header:
            g.write(line)
        df.to_csv(g, index=False, sep=" ")


for tile_file in tqdm(all_tile_files):
    if tile_file == "G23_T070.fld":
        continue

    new_header = fix_header(tile_file)
    df = pd.read_csv(tile_file, comment="#")
    df["RA"] = df["RA"] % 360
    copy_to_backup(tile_file)
    write_correct_tile_file(tile_file, df, new_header)

for guide_file in tqdm(all_guide_files):
    if guide_file == "G23_GT070.fld":
        continue

    new_header = fix_header(guide_file)
    df = pd.read_csv(guide_file, delim_whitespace=True, comment="#")
    df["RA"] = df["RA"] % 360
    copy_to_backup(guide_file)
    write_correct_guide_file(guide_file, df, new_header)
