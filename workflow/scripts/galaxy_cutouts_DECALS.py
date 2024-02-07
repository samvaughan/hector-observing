import pandas as pd
from scipy import misc
from tqdm import tqdm
import requests
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import warnings
from astropy.utils.exceptions import AstropyWarning
import misc  # noqa
import tempfile

smk = snakemake  # noqa

# Read the tile file and set up the output directories
df = pd.read_csv(smk.input.tile_file, skiprows=11)
base_folder = Path(smk.output.output_folder)
base_folder.mkdir(exist_ok=True)
print(f"Getting images from {smk.params.image_source}")

for index, row in tqdm(df.iterrows(), total=len(df)):
    ra = row.RA
    dec = row.DEC
    name = str(row.Hexabundle)
    CATID = row.ID

    # The skyfibres have nan as their hexabundle name, so make this their
    # Skyfibre ID instead
    if name == "nan":
        name = row.ID

    if "Sky" in name:
        folder = base_folder / "SkyFibres"
    elif "GS" in name:
        folder = base_folder / "GuideBundles"
    else:
        folder = base_folder / "Hexabundles"

    folder.mkdir(exist_ok=True)

    if smk.params.image_source != "DECALS":
        raise NotImplementedError("Not yet sorted the SDSS images...")

    # URL for the jpg cutouts, for guides and skies
    jpg_url = f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra:.5f}&dec={dec:.5f}&layer=ls-dr9&pixscale=0.262&bands=grz"

    # Make the sky fibre cutouts very small
    if "Sky" in name:
        jpg_url = f"{jpg_url}&size=20"

    # URL for the fits images we need for the galaxies
    fits_url = f"https://www.legacysurvey.org/viewer/fits-cutout?ra={ra:.5f}&dec={dec:.5f}&layer=ls-dr9&pixscale=0.262&bands=grz"

    # If we're looking at a guide or a sky fibre, just get a jpg
    # If we're looking at a hexabundle, download the fits and add the hector bundle footprint
    if len(name) > 1:
        img_data = requests.get(jpg_url, verify=False).content
        with open(f"{folder}/{name}.jpg", "wb") as handler:
            handler.write(img_data)
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_fits_file = Path(f"{tmpdirname}/tmp.fits")
            tmp_jpg_file = Path(f"{tmpdirname}/tmp.jpg")
            # Save the jpg and the fits data
            img_data = requests.get(jpg_url, verify=False).content
            with open(tmp_jpg_file, "wb") as handler:
                handler.write(img_data)
            fits_data = requests.get(fits_url, verify=False).content
            with open(tmp_fits_file, "wb") as handler:
                handler.write(fits_data)

            # Now read it in again
            hdu = fits.open(tmp_fits_file)
            image = plt.imread(tmp_jpg_file)
        # Get the WCS from the fits header
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            wcs = WCS(hdu[0].header)
        # make the image
        fig, ax = misc.make_image_with_info(image, wcs, row)
        # And save
        fig.savefig(f"{folder}/{name}.png", bbox_inches="tight")
        plt.close("all")

Path(smk.output.finished_flag).touch()
