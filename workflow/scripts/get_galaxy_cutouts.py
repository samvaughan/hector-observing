import pandas as pd
from tqdm import tqdm
import requests
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from astropy import units as u
import warnings
from astropy.utils.exceptions import AstropyWarning


def make_image_with_info(image, wcs, row):
    # Do the plotting
    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs.celestial))
    ax.imshow(image, origin="lower")

    # Add the extra info
    ax.annotate(
        f"Hexabundle {row.Hexabundle}:\n{row.ID}\nMstar={row.Mstar:.2f}\nz={row.z:.4f}\nRe={row.Re:.2f}",
        xy=(0.03, 0.97),
        xycoords="axes fraction",
        bbox=dict(boxstyle="round", fc="0.4", ec="k", alpha=0.9),
        ha="left",
        va="top",
        color="w",
    )
    # Now add the bundle outline
    hexabundle_radius = hexabundle_radius_dict[row.Hexabundle]
    bundle_outline = SphericalCircle(
        center=(ra, dec) * u.deg,
        radius=hexabundle_radius * u.arcsec,
        edgecolor="white",
        facecolor="none",
        linewidth=1.5,
        transform=ax.get_transform("fk5"),
    )
    ax.add_patch(bundle_outline)

    ax.axis("off")

    return fig, ax


smk = snakemake # noqa

# Read the tile file and set up the output directories
df = pd.read_csv(smk.input.tile_file, skiprows=11)
base_folder = Path(smk.output.output_folder)
base_folder.mkdir(exist_ok=True)
print(f"Getting images from {smk.params.image_source}")

hexabundle_radius_dict = dict(
    A=13,
    B=13,
    C=11.2,
    D=19 / 2,
    E=15.5 / 2,
    F=15.5 / 2,
    G=15.5 / 2,
    H=12.1 / 2,
    I=19 / 2,
    J=15.5 / 2,
    K=15.5 / 2,
    L=15.5 / 2,
    M=15.5 / 2,
    N=15.5 / 2,
    O=15.5 / 2,
    P=15.5 / 2,
    Q=15.5 / 2,
    R=15.5 / 2,
    S=15.5 / 2,
    T=15.5 / 2,
    U=12.1 / 2,
)

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
        img_data = requests.get(jpg_url).content
        with open(f"{folder}/{name}.jpg", "wb") as handler:
            handler.write(img_data)
    else:
        tmp_fits_file = Path("tmp.fits")
        tmp_jpg_file = Path("tmp.jpg")
        # Save the jpg and the fits data
        img_data = requests.get(jpg_url).content
        with open(tmp_jpg_file, "wb") as handler:
            handler.write(img_data)
        fits_data = requests.get(fits_url).content
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
        fig, ax = make_image_with_info(image, wcs, row)
        # And save
        fig.savefig(f"{folder}/{name}.png", bbox_inches="tight")
        plt.close("all")
        tmp_fits_file.unlink()
        tmp_jpg_file.unlink()

Path(smk.output.finished_flag).touch()
