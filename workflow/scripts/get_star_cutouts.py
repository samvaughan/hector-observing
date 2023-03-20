"""ps1bulk.py: Get PS1 stack FITS cutout images at a list of positions"""
 
import numpy as np
from astropy.table import Table
import requests
from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt
from pathlib import Path
 
# ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
# fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
 
def getimages(ra,dec,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,filters=filters)
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"ra={ra}&dec={dec}&size={size}&format={format}")
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    url = geturl(ra,dec,size=size,filters=filters,output_size=output_size,format=format,color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="jpg"):
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = geturl(ra,dec,size=size,filters=filter,output_size=output_size,format=format)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))
    return im


if __name__ == "__main__":

    import pandas as pd
    import string
    import numpy as np

    tile_file = snakemake.input.tile_file
    hexabundle_filenames = snakemake.output.hexabundle_cutouts
    guidebundle_filenames = snakemake.output.guidestar_cutouts
    
    df = pd.read_csv(tile_file, skiprows=11)
    size_arcsec = 30
    pix_size = 0.25
    size = int(size_arcsec / pix_size)

    tick_marks = np.arange(0, size, 20)
    tick_labels = [f'{s * pix_size}"' for s in tick_marks]

    hexas = df.loc[df['type'].isin([0, 1])].sort_values('Hexabundle')
    guides = df.loc[df['type'] == 2]

    for hexabundle_filename, (index, row) in zip(hexabundle_filenames, hexas.iterrows()):
        
        ra = row.RA
        dec = row.DEC
        cim = getcolorim(ra,dec,size=size,filters="grz")

        # Plot the cutout
        fig, ax = plt.subplots(constrained_layout=True)
        ax.imshow(cim, origin='lower')
        ax.set_title(f"Hexabundle {row.Hexabundle}")
        ax.plot([20, 100], [20, 20], c='w', linewidth=3.0)
        ax.plot([20, 20], [20, 100], c='w', linewidth=3.0)
        ax.annotate(text='20"', xy=(60,20), xytext=(0, -10), textcoords='offset pixels', ha='center', va='top', c='w')
        ax.annotate(text='20"', xy=(20,60), xytext=(-10, 0), textcoords='offset pixels', ha='right', va='center', c='w')
        ax.set_axis_off()
        
        print(f"Saving Cutout {row.Hexabundle} to {hexabundle_filename}")
        fig.savefig(hexabundle_filename, bbox_inches='tight')
        plt.close('all')

    for guide_fname, (index, row) in zip(guidebundle_filenames, guides.iterrows()):

        ra = row.RA
        dec = row.DEC
        cim = getcolorim(ra,dec,size=size,filters="grz")

        # Plot the cutout
        fig, ax = plt.subplots(constrained_layout=True)
        ax.imshow(cim, origin='lower')
        ax.set_title(f"Hexabundle {row.Hexabundle}")
        ax.plot([20, 100], [20, 20], c='w', linewidth=3.0)
        ax.plot([20, 20], [20, 100], c='w', linewidth=3.0)
        ax.annotate(text='20"', xy=(60,20), xytext=(0, -10), textcoords='offset pixels', ha='center', va='top', c='w')
        ax.annotate(text='20"', xy=(20,60), xytext=(-10, 0), textcoords='offset pixels', ha='right', va='center', c='w')
        ax.set_axis_off()
        
        print(f"Saving Cutout {row.Hexabundle} to {guide_fname}")
        fig.savefig(guide_fname, bbox_inches='tight')
        plt.close('all')

    Path(snakemake.output.cutouts_finished_flag).touch()