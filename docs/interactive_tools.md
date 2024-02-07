# Interactive tools

There are two interactive tools which have been built to help with running these tools. 

## Selecting galaxies to create a tile from scratch

The script `interactive_apps/HectorTilingApp.py` launches a [dash/plotly](https://dash.plotly.com/) app to create a Hector tile from scratch by clicking on galaxies, guide stars and standard stars. These tiles can then be saved and run through the pipeline to make a set of files which can be observed at the telescope. This is particularly important for the beginning of every observing run, where two special "SNAFU" tiles are needed for the first and second halves of the first night.

### Starting the app

The app needs a catalogue of galaxy targets, standard stars and guide stars for a given region of sky. You also need to pass the name of the region you're running. At the command line, run:

```bash
python interactive_apps/HectorTilingApp.py /path/to/target_catalogue.csv /path/to/standard_star_catalogue.csv /path/to/guide_star_catalogue.csv field_name
```

Note that the code expects that the standard star and guide star catalogues name have "standard_star" and "guide_star" in their filenames, and will error if this isn't the case. This was added to catch the bug where you pass these two catalogues to the code in the wrong order!

You should get the following in the terminal:

![Terminal output from running HectorTilingApp.py](img/tile_app_terminal.png)

Copy the web address (in this case `http://127.0.0.1:8051/`) into a browser. Note that the exact port number you see might be different!

You should see the following:

