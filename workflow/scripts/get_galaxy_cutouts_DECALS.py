import pandas as pd
import os
from tqdm import tqdm
import requests
from pathlib import Path
import shutil

df = pd.read_csv(snakemake.input.tile_file, skiprows=11)

base_folder = Path(snakemake.params.output_folder)
base_folder.mkdir(exist_ok=True)

for index, row in tqdm(df.iterrows(), total=len(df)):

    ra = row.RA
    dec = row.DEC
    name = str(row.Hexabundle)
    CATID = row.ID

    if name == 'nan':
        name = row.ID

    if "Sky" in name:
        folder = base_folder / "SkyFibres"
    elif "GS" in name:
        folder = base_folder / "GuideBundles"
    else:
        folder = base_folder / "Hexabundles"

    folder.mkdir(exist_ok=True)

    if not "Sky" in name:
        url = f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&layer=decals-dr7&pixscale=0.27&bands=grz"
    else:
        url = f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&layer=decals-dr7&pixscale=0.27&bands=grz&size=20"
    img_data = requests.get(url).content
    with open(f'{folder}/{name}.jpg', 'wb') as handler:
        handler.write(img_data)
    #print(f"wget --read-timeout=5 -O  https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&zoom=14&layer=ls-dr9")

Path(snakemake.output.finished_flag).touch()