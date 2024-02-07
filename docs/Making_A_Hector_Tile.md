# Making a Hector Tile

Below is a detailed set of instructions for making the files we need at the telescope. 



## At the beginning of a run

- Make a new folder with the format `{run_start_date}_{run_end_date}` in the ```resources``` folder: e.g. ```20240306_20240318``` for an observing run between March 3rd and March 18th in 2024.
- Take the tiles from the Hector Tiling code and transfer them into this directory.
- Make a new file called ```{run_start_date}_{run_end_date}_galaxy_tiles.csv```, e.g. ```20240306_20240318_galaxy_tiles.csv```. This must have the columns named ```field,tile_number,filename,guide_filename,image_source```. Each tile you want to run the pipeline on must be included here on a separate row.
- Pick a tile to make. Fill out the row in the ```{run_start_date}_{run_end_date}_galaxy_tiles.csv``` file with the correct information. The `image_source` column should be '`DECALS`' for all regions other than the G23 field where it should be '`KIDS`'. All filenames are relative to the main observing folder, i.e. they should begin with `resources/{run_start_date}_{run_end_date}/` etc.
- Now update the config file, which lives in `config/`. TODO: Write up explanation for this file.


## Making a tile

- Make sure the `hector` environment has been activated. TODO: Add page about the Hector environment. 
- Make sure the sky masks are available (i.e. plug in your external hard drive!)
- From the main folder, run ```snakemake -npr --cores 1 --configfile config/{run_start_date}_{run_end_date}.yaml -- results/{run_start_date}_{run_end_date}/Upload/{tile_file_name}.tar.gz```. This will show you the commands which `snakemake` is about to run.
- To actually execute these, run `snakemake --cores 1 --configfile config/{run_start_date}_{run_end_date}.yaml -- results/{run_start_date}_{run_end_date}/Upload/{tile_file_name}.tar.gz` (i.e. remove the `-npr` bit).


All going well, the pipeline will continue and the `R` plotting window will appear like so:

![Configuration code plot from the R code](img/R_configuration_code.png)

This can take 30-60 minutes to complete, so probably best to leave it running in the background while you do something else!

When the R code is complete, the allocation code will run and assign which Hexabundle (A through to U) will be placed on which target. 
