from astropy.io.votable import parse_single_table
import requests
from pathlib import Path
from tqdm import tqdm 
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
from astropy.wcs import WCS
import pandas as pd

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
    radius = 100 / 60 / 60
    url = f"https://datacentral.org.au/vo/sia2/query?FORMAT=fits&POS=CIRCLE {ra} {dec} {radius:.4f}"

    #filename to store the results of the query
    vot = Path('test_table.xml')
    #perform the query
    r = requests.get(url)
    #write the results of the query to a file
    # This seems dumb. Why is it necessary? I can't seem to get astropy to read this VOTable without writing it to a file, however...
    with open(vot,'wb') as f:
        f = open(vot,'wb')
        f.write(r.content)

    #open the table with astropy's votable package
    df = parse_single_table(vot).to_table().to_pandas()
    # Delete the temporary file
    vot.unlink()

    # Now donwload the u, g and r files and make a colour cutout. 
    required_colours = ['u', 'g', 'r']
    colours = df.loc[df.band_name.isin(required_colours)]

    for index, row in tqdm(colours.iterrows(), total=len(required_colours)):
        # obs_id = str(data['obs_id'][i])
        # obs_collection = str(data['obs_collection'][i])
        # facility_name =  str(data['facility_name'][i])
        # band_name = str(data['band_name'][i])
        band_name = row.band_name
        #URL of the cutout
        access_url = row.access_url
        #filename to download to
        ofname = Path(f"tmp_band_{band_name}.fits")
        #print out a wget command to download the data
        download_cmd = f"wget -O {ofname} \"{access_url}\""
        r = requests.get(access_url)
        with open(ofname, 'wb') as f:
            f.write(r.content)

    r_hdu = fits.open("tmp_band_r.fits")
    r_header = r_hdu[0].header
    r_img = r_hdu[0].data / 1e-12
    g_img = fits.open("tmp_band_g.fits")[0].data / 1e-12
    b_img = fits.open("tmp_band_u.fits")[0].data / 1e-12

    #Get the image WCS
    wcs = WCS(r_header)

    rgb = make_lupton_rgb(image_b=b_img, image_g=g_img, image_r=r_img, Q=10, stretch=5)
    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
    ax.imshow(rgb)
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_major_formatter('d.ddd')
    lat.set_major_formatter('d.ddd')

    # Pixel size
    pix_scale = np.abs(r_header['PC1_1'] * 60 * 60)
    ax.plot()

    trans = (fig.dpi_scale_trans +
         transforms.ScaledTranslation(xdata[0], ydata[0], ax.transData))