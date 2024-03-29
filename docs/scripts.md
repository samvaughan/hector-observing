# Scripts

There are several scripts included in the repository which perform helpful tasks. Running any of the scripts below with the single argument `-h` will print out some useful text describing each argument and the script's purpose.

## Swapping Hexabundles

It is often important to swap the Hexabundles which are assigned by the allocation code between galaxies. You can do this using the script [`swap_hexabundles.py`](https://github.com/samvaughan/hector-observing/blob/main/swap_hexabundle.py).

The inputs for this script are the "tile" and "robot" files in the "FinalOutputs" folder. You then supply the Hexabundles you wish to swap as a pair of letters separated by a dash. Multiple pairs of Hexabundles may be passed at the same time. For example, to swap the galaxies placed on Hexabundles J and K for the tile G15 tile 220, you
d run:

```bash
python swap_hexabundle.py results/20240306_20240318/TilingOutputs/G15/FinalOutputs/Tile_G15_T220.csv results/20240306_20240318/TilingOutputs/G15/FinalOutputs/Robot_G15_T220.csv J-K
```

The script makes a backup of the original file, which is saved to the original filename with `.backup` appended to the end. If you'd like to undo your changes, just run 

```bash
mv path/to/Tile_G15_T220.csv.backup path/to/Tile_G15_T220.csv
```

## Selecting Guide Stars for a tile

It's often useful to be able to select which guide stars you'd like to include in a tile. We can do that using [`select_guides_from_file.py`](https://github.com/samvaughan/hector-observing/blob/main/select_guides_from_file.py).

The input is a guide file _which hasn't been configured by the R code_, i.e. it should have `NOT_CONFIGURED` in the filename. You also pass the integer IDs of the 6 guide stars you'd like to select (i.e. the numbers which are shown in the interactive configuration tool).

For example, to select the guide stars 1, 3, 5, 7, 9 and 11 from G15 tile 220, you'd run:

```bash
python select_guides_from_file.py results/20240306_20240318/TilingOutputs/G15/Configuration/Guides_G15_tile_220_NOT_CONFIGURED.csv 1 3 5 7 9 11
```

The script makes a backup by copying the original file to a file with the same name with `.backup` appended to the end. If you'd like to undo your changes, you can just run:

```bash
mv path/to/Guides_G15_tile_220_NOT_CONFIGURED.csv.backup pat/to/Guides_G15_tile_220_NOT_CONFIGURED.csv
```

## Selecting Standard Stars for a tile

Occasionally you'll need to select which standard stars go into a tile (although much less regularly than for the guides!). You can do this by using the `select_standards_from_file.py` script.

The input is a hexabundle file _which hasn't been configured by the R code_, i.e. it should have `NOT_CONFIGURED` in the filename. You then pass two numbers as integers on the command line. These numbers refer to the standard stars in order of distance from the centre, _not_ their ID number. We use python 0-based indexing, so passing `0` and `1` would select the two stars closest to the field centre.

So to select the two standards closest to the field centre for G15 tile 220, you'd run:

```bash
python select_standards_from_file.py results/20240306_20240318/TilingOutputs/G15/Configuration/Hexas_G15_tile_220_NOT_CONFIGURED.csv 0 1
```

This script makes a backup of the `Hexas_` file and save it to the original filename with `.backup` appended to the end.

Once you've selected your stars, you can then use the interactive configuration tool to load this edited file and check that these two stars aren't clashing with each other or another galaxy. If you need to pick different stars and undo your changes, you can run:

```bash
mv path/to/Hexas_G15_tile_220_NOT_CONFIGURED.csv.backup pat/to/Hexas_G15_tile_220_NOT_CONFIGURED.csv
```

and use the script again.