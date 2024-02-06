# Making Hector Tiles for Observing 

A first go at making docs for the process of making tiles for the Hector Survey.

## Tiling Steps

1. Removing all guide stars which are outside the central 70% of the field (`select_guides_within_set_radius.py`).
2. Copying the tile and guide tile to the `results/{run_start_date}_{run_end_date}/TilingOutputs/{region_name}/Tiles` folder.
3. Run the distortion correction code, which turns right-ascension and declination into plate $x$ and $y$ coordinates and finds the appropriate sky-fibre positions (`run_distortion_correction_code.py`).
4. Run the configuration code, which finds magnet positions across the plate (`HECTOR_ClusterFieldsTest.R`).
5. Fix the file headers into the correct format (`update_header.py`).
6. Run the allocation code, which decides which hexabundle is placed on which galaxy. This code also creates the final output files in the format required by the computers at the AAT (`allocate_tile_for_galaxies.py`).
7. Download the galaxy cutouts for each Hexabundle, Guide bundle and sky fibre (`get_galaxy_cutouts.py`).
8. Run a series of sanity-checks to catch a number of simple mistakes (e.g. we check the number of guide bundles in each tile is 6, etc). Also use the Hector observing database to see how many galaxies in this tile have been observed before (`scripts/verify_tile_properties.py`).
9. Add the tile to our database of configured tiles (`add_configured_tile_to_database.py`).
10. Make the plot the observers used at the top-end of the telescope (`make_updated_hexabundle_diagrams.py`).
11. Package everything up in a nice format which can be easily uploaded to the cloud.