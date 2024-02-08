# Making Hector Tiles for Observing 

Documentation for the code which makes the tile and robot files used at the telescope when observing for the Hector Galaxy Survey.

The steps to make a tile are below. These are performed are performed automatically using the workflow tool [`snakemake`](https://snakemake.github.io/). See [Making a Hector Tile](Making_A_Hector_Tile.md) for more details!

## Quick Overview of Pipeline Steps

1. Remove all guide stars which are outside the central 70% of the field. ([`select_guides_within_set_radius.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/select_guides_within_set_radius.py)).
2. Copy the tile and guide tile to the `results/{run_start_date}_{run_end_date}/TilingOutputs/{region_name}/Tiles` folder.
3. Run the distortion correction code, which turns right-ascension and declination into plate $x$ and $y$ coordinates and finds the appropriate sky-fibre positions ([`run_distortion_correction_code.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/apply_DC_correction.py)).
4. Run the configuration code, which finds magnet positions across the plate ([`HECTOR_ClusterFieldsTest.R`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/HECTOR_ClusterFieldsTest.R)).
5. Fix the file headers into the correct format ([`fix_header_after_configuration.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/fix_header_after_configuration.py)).
6. Run the allocation code, which decides which hexabundle is placed on which galaxy. This code also creates the final output files in the format required by the computers at the AAT ([`allocate_tile_for_galaxies.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/allocate_tile_for_galaxies.py)).
7. Download the galaxy cutouts for each Hexabundle, Guide bundle and sky fibre ([`get_galaxy_cutouts_DECALS.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/get_galaxy_cutouts_DECALS.py) or [`get_galaxy_cutouts_KIDS.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/get_galaxy_cutouts_KIDS.py)).
8. Run a series of sanity-checks to catch a number of simple mistakes (e.g. we check the number of guide bundles in each tile is 6, etc). Also use the Hector observing database to see how many galaxies in this tile have been observed before ([`verify_tile_properties.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/verify_tile_properties.py)).
9. Add the tile to our database of configured tiles ([`add_configured_tile_to_database.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/add_configured_tile_to_database.py)).
10. Make the plot the observers used at the top-end of the telescope ([`make_updated_hexabundle_diagrams.py`](https://github.com/samvaughan/hector-observing/blob/main/workflow/scripts/make_hexabundle_plot.py)).
11. Package everything up in a nice format which can be easily uploaded to the cloud.

You can now upload the `{tile_ID}.tar.gz` file to the Data Central cloud.

## Tips, Tricks and things to keep in mind

- The `snakemake` workflow system depends on the _filenames_ of the inputs and outputs to each rule. A lot of the headaches I've had with `snakemake` are when I think I've placed a file in one directory but it's actually in another, or when a file I've saved has a typo, etc. 
- If you _really_ need to, you can fool some of the `snakemake` steps by seeing which files a given rule creates and making that file yourself. There _might_ be cases when this is necessary (e.g. to bypass a check we have in place) but you should really understand what you're doing in this case.
