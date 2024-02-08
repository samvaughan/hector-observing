# FAQs and troubleshooting

Below are some frequently asked questions or bugs that have arisen during our use of the pipeline.

## The allocation code fails and says "Some magnets are fully blocked!"

This means that there is no way for the Robot to place the magnets in the desired positions. The robot arm itself has a set clearance when it picks up and puts down magnets, and so some magnet configurations which don't clash aren't actually possible to put on the plate.

In this case, you'll need to go back and use the interactive configuration tool to change the magnet positions to have a bit more clearance. The error message will let you know which magnets are affected.

## Snakemake has an error which mentions _"`The files below seem to be incomplete`"_.

This happens when a `snakemake` task doesn't finish properly- often because it's been killed by someone pressing `Ctrl-C`. It's basically saying "I don't know if the files below are finished or not!".

If you want to re-run the task which `snakemake` was in the middle of before it was killed, add `--rerun-incomplete` to the beginning of your command (e.g `snakemake --rerun-incomplete --configfile ...` etc). 

If you _don't_ want to re-run the command that was killed, there are two options. I've found the best thing to do is delete the `.snakemake` folder: ```rm -r .snakemake```. This deletes `snakemake's` history and it just forgets that it was halfway through a step. You can also follow the commands that `snakemake` suggests, about "cleaning up metadata", but I've found that this often doesn't work for some reason... (see e.g. [this](https://github.com/snakemake/snakemake/issues/1497) github issue).

## Snakemake complains about _"`Index 0 is out of bounds for axis with size 0`"_ and mentions the _`_get_tile_filename`_ function.

This happens because Snakemake is trying to select a row from your `{start_date}_{end_date}_galaxy_tiles.csv` but there aren't any matches. Check that you've typed the tile filename exactly right on the command line (it's very easy to type `G12_tile_220` instead of `G15_tile_220`, for example!)

## Snakemake says "there is no rule to make {filename}"

This happens when `snakemake` can't figure out how to create the output file you've asked for. 

The most common cause of this is a typo in the filename you've passed it (e.g. asking it to make `G155_tile_220.tar.gz` will fail, for example). Another possible cause would be if you're asking it to make a tile which you haven't entered in the `resources/{run_dates}/{run_dates}_galaxy_tiles.csv` file. Check that the tile is located there and doesn't have any typos.

## Snakemake says "Looks like we have a mismatch between the tile numbers and the filenames somewhere"

This occurs when one (or more) of the rows in the `resources/{run_dates}/{run_dates}_galaxy_tiles.csv` is incorrect. If you tell `snakemake` to use a tile file from one tile number and a guide tile from another (e.g. one row references two different numbers), or if your tile number column in that file doesn't match the filenames, you'll see this error. Go back and check that all of the rows in that file are correct!

## The distortion correction code complains about "Reading past the end of a file"

This means a sky mask fits file the code is trying to read is corrupted or damaged. 

The distortion correction code reads _all_ fits files in the folder path you give it, so make sure there aren't any files in there which aren't sky masks. We've seen this error before when the distortion correction code tried to read a hidden file called `._segmap_waves_18.0_-34.1_1.41_1.17.fits`. After deleting this file and all the others like it, everything worked as expected.

## Snakemake complains about "Error Tokenizing data: Expecting 5 lines of data, saw {n}"

This is almost certainly caused by a typo in the `resources/{run_dates}/{run_dates}_galaxy_tiles.csv`. Check that you haven't added or forgotten a comma somewhere.