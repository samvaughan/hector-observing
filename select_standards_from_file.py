"""
This is a script to select a subsample of 2 standard stars from a "...NOT_CONFIGURED.csv" file. It edits the file, deletes the unnecessary stars and saves it to the same name (after making a backup). We can then load the same NOT_CONFIGURED.csv file into the configuration app and choose the correct angles to make a valid Hector configuration.

Note that the indexing is R-based, so starts from 1! I.e. you should pass the numbers that you see in the configuration app into this file. It then subtracts one from each of them to get the correct rows in the file.
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

assert hexa_filename.name.startswith("Hexas_"), "This code is only supposed to run on a Guide file, not anything else!"
assert 'NOT_CONFIGURED' in hexa_filename.name, "This code must be run on a NOT_CONFIGURED file!"

stringIDs = [str(i) for i in standard_IDs]
print(f"Selecting standard stars with IDs {', '.join(stringIDs)} (1-based indexing) from file {hexa_filename}")
print(f"\tMaking a backup at {backup_fname}")
shutil.copy(hexa_filename, backup_fname)

print("\tSelecing the standards")
df = pd.read_csv(hexa_filename, comment='#')
df_targets = df.loc[df['type'] == 1]

python_standard_IDs = [index - 1 for index in standard_IDs]
python_standard_IDs_adjusted = np.array(python_standard_IDs) + len(df_targets)

df_stars = df.loc[python_standard_IDs_adjusted, :]
assert np.all(df_stars['type'] == 0.0), "We have an error with the standard stars- they should have type 0"

df_skies = df.loc[~np.isfinite(df['type'])]

selected = pd.concat((df_targets, df_stars, df_skies))

selected.to_csv(hexa_filename, index=True)
print("Done!")
