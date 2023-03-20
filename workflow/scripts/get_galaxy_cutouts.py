import pandas as pd
import os
from tqdm import tqdm
import requests
from pathlib import Path
import shutil

df = pd.read_csv(snakemake.input.tile_file, skiprows=11)

base_folder = Path(snakemake.output.output_folder)
base_folder.mkdir(exist_ok=True)


print(f"Getting images from {snakemake.params.image_source}")

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
        image = Path(f"/Users/samvaughan/Science/Hector/Imaging/hsc_imaging/{CATID}_rgb.png")
        if image.exists():
            shutil.copy(image, f'{folder}/{name}.png')
            continue
        else:
            if snakemake.params.image_source == 'SDSS':
                url = f"https://test.preprod.skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra={ra:.5f}&dec={dec:.5f}&width=512&height=512&scale=0.4"
            elif snakemake.params.image_source == 'DECALS':
                 url = f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&layer=ls-dr9&pixscale=0.27&bands=grz"
            else:
                raise NameError(f"Image source must be one of SDSS or DECALS! Currently it's {snakemake.params.image_source}")
    else:
        if snakemake.params.image_source == 'SDSS':
            url = f"https://test.preprod.skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra={ra:.5f}&dec={dec:.5f}&width=300&height=300&scale=0.4&height=10"
        elif snakemake.params.image_source == 'DECALS':
             url = f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&layer=ls-dr9&pixscale=0.27&bands=grz&size=20"
        else:
            raise NameError(f"Image source must be one of SDSS or DECALS! Currently it's {snakemake.params.image_source}")

    img_data = requests.get(url).content
    with open(f'{folder}/{name}.jpg', 'wb') as handler:
        handler.write(img_data)
    #print(f"wget --read-timeout=5 -O  https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&zoom=14&layer=ls-dr9")

Path(snakemake.output.finished_flag).touch()