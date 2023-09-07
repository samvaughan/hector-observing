"""
This is a script to select a subsample of 2 standard stars from a "...NOT_CONFIGURED.csv" file. It edits the file, deletes the unnecessary stars and saves it to the same name (after making a backup). We can then load the same NOT_CONFIGURED.csv file into the configuration app and choose the correct angles to make a valid Hector configuration.

Note that the stars you select here are **sorted by their distances from the centre!!**. I.e. selecting stars 0 and 1 will give you the two stars which are closest to the middle. 
"""
import pandas as pd
import argparse
from pathlib import Path
import shutil
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("hexa_filename")
parser.add_argument("standard_IDs", nargs=2, type=int)

args = parser.parse_args()
hexa_filename = Path(args.hexa_filename)
standard_IDs = args.standard_IDs

backup_fname = hexa_filename.parent / (hexa_filename.name + ".backup")

assert hexa_filename.name.startswith(
    "Hexas_"
), "This code is only supposed to run on a Guide file, not anything else!"
assert (
    "NOT_CONFIGURED" in hexa_filename.name
), "This code must be run on a NOT_CONFIGURED file!"

stringIDs = [str(i) for i in standard_IDs]
print(
    f"Selecting standard stars with IDs {', '.join(stringIDs)} (1-based indexing) from file {hexa_filename}"
)

print("\tSelecing the standards")
df = pd.read_csv(hexa_filename, comment="#")
df_targets = df.loc[df["type"] == 1]

df_stars = df.loc[df["type"] == 0].copy()
df_stars["radius"] = np.sqrt(df_stars["MagnetX"] ** 2 + df_stars["MagnetY"] ** 2)
df_stars_sorted = df_stars.sort_values("radius")


python_standard_IDs = [index for index in standard_IDs]
# python_standard_IDs_adjusted = np.array(python_standard_IDs) + len(df_targets)

df_stars_selected = df_stars_sorted.iloc[python_standard_IDs, :]
assert np.all(
    df_stars_selected["type"] == 0.0
), "We have an error with the standard stars- they should have type 0"

df_skies = df.loc[~np.isfinite(df["type"])]
df_stars_selected = df_stars_selected.drop("radius", axis=1)

selected = pd.concat((df_targets, df_stars_selected, df_skies))

print(f"\tMaking a backup at {backup_fname}")
shutil.copy(hexa_filename, backup_fname)
selected.to_csv(hexa_filename, index=True)
print("Done!")
