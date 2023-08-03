"""
Get image cutouts from Data Central. We have to get individual r, g and b fits images and then 
make them into a colour image ourselves using make_lupton_rgb
"""
import astropy.coordinates as coord
from astropy import units as u
import misc
import numpy as np
import pandas as pd
from pyvo.dal.sia2 import SIAService
import requests
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
from astropy.wcs import WCS
from pathlib import Path
import tempfile
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

smk = snakemake  # noqa

# Read the tile file and set up the output directories
df = pd.read_csv(smk.input.tile_file, skiprows=11)
base_folder = Path(smk.output.output_folder)
base_folder.mkdir(exist_ok=True)
print(f"Getting images from {smk.params.image_source}")

# URL for the AAO DataCentral SIA2 service
url = "https://datacentral.org.au/vo/sia2/query"
service = SIAService(url)

for index, row in tqdm(df.iterrows(), total=len(df)):
    # we request 1 arcmin square images for the galaxies and stars
    image_radius = 30.0 / 3600

    # The skyfibres have nan as their hexabundle name, so make this their
    # Skyfibre ID instead.
    # Also make their image very small
    name = str(row.Hexabundle)
    if name == "nan":
        name = row.ID
        image_radius = 0.014

    if "Sky" in name:
        folder = base_folder / "SkyFibres"
    elif "GS" in name:
        folder = base_folder / "GuideBundles"
    else:
        folder = base_folder / "Hexabundles"

    folder.mkdir(exist_ok=True)

    ra = coord.Angle(row.RA, unit=u.degree)
    dec = coord.Angle(row.DEC, unit=u.degree)
    # position and radius to query
    pos = (ra.degree * u.deg, dec.degree * u.deg, image_radius * u.deg)

    # Query for r, g and u filters
    filters = ["u", "g", "r"]
    custom = {}
    custom["FILTER"] = filters
    # Do the SIAP search and convert to a pandas dataframe
    results = service.search(pos=pos, **custom).to_table().to_pandas()
    # Only keep the rows which are from the GAMA pdr, and drop duplicates
    results = results.loc[results.obs_collection == "gama_pdr"]
    results = results.loc[
        ~results.duplicated(subset=["band_name", "em_min", "em_max"], keep="first")
    ]
    if len(results) == 0:
        print(f"No matches for {pos}! Skipping {name}...")
        continue

    hdus = {}
    usable_files = True
    for index, image in results.iterrows():
        access_url = image.access_url

        with tempfile.TemporaryDirectory() as tmpdirname:
            fname = f"{tmpdirname}/{name}_{image.facility_name}_{image.band_name}.fits"
            # print(f"Downloading {fname}...")
            fits_data = requests.get(image.access_url).content
            with open(fname, "wb") as handler:
                handler.write(fits_data)
            try:
                hdus[image.band_name] = fits.open(fname)
            except OSError:
                print(f"No file available for {name}! Skipping...")
                usable_files = False

    # If we don't have the fits files to download, just skip this row
    if not usable_files:
        continue
    # Make an RGB image
    r_hdu = hdus["r"]
    g_hdu = hdus["g"]
    b_hdu = hdus["u"]

    r = r_hdu[0].data
    g = g_hdu[0].data
    b = b_hdu[0].data

    # Get an average intensity image
    intensity = (r + g + b) / 3
    m = np.nanmean(intensity)
    minimum = np.nanpercentile(intensity / m, 0.01)
    rgb_image = make_lupton_rgb(
        r / m, 1.5 * g / m, 3 * b / m, stretch=20, Q=5, minimum=minimum
    )
    plt.imshow(rgb_image)
    # Get the wcs
    wcs = WCS(r_hdu[0].header).celestial

    # make the image
    if (row["type"] == 1) or (row["type"] == 0):
        fig, ax = misc.make_image_with_info(rgb_image, wcs, row)
    else:
        fig, ax = plt.subplots()

        if "Sky" in name:
            # If we're looking at a skyfibre, take a 5" cutout of the big image
            position = SkyCoord(ra, dec, frame="icrs")
            size = u.Quantity((1.5, 2.5), u.arcsec)
            # Also only plot the red image, and use a B/W colormap
            rgb_image = Cutout2D(
                rgb_image[..., 0], position, (5 * u.arcsec, 5 * u.arcsec), wcs=wcs
            ).data
            ax.imshow(rgb_image, cmap="binary")
        else:
            ax.imshow(rgb_image)
        ax.axis("off")
    fig.savefig(f"{folder}/{name}.png", bbox_inches="tight")
    plt.close("all")

Path(smk.output.finished_flag).touch()
