"""
This is a script to select a subsample of 6 guide stars from a "...NOT_CONFIGURED.csv" file. It edits the file, deletes the unnecessary stars and saves it to the same name (after making a backup). We can then load the same NOT_CONFIGURED.csv file into the configuration app and choose the correct angles to make a valid Hector configuration.

Note that the indexing is R-based, so starts from 1! I.e. you should pass the numbers that you see in the configuration app into this file. It then subtracts one from each of them to get the correct rows in the file.
"""
import pandas as pd
import argparse
from pathlib import Path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("guide_filename")
parser.add_argument("guide_IDs", nargs=6, type=int)

args = parser.parse_args()
guide_filename = Path(args.guide_filename)
guide_IDs = args.guide_IDs

python_guide_IDs = [index - 1 for index in guide_IDs]

backup_fname = guide_filename.parent / (guide_filename.name + ".backup")

assert guide_filename.name.startswith("Guides_"), "This code is only supposed to run on a Guide file, not anything else!"
assert 'NOT_CONFIGURED' in guide_filename.name, "This code must be run on a NOT_CONFIGURED file!"

stringIDs = [str(i) for i in guide_IDs]
print(f"Selecting guide stars with IDs {', '.join(stringIDs)} (1-based indexing) from file {guide_filename}")
print(f"\tMaking a backup at {backup_fname}")
shutil.copy(guide_filename, backup_fname)

print("\tSelecing the guides")
df = pd.read_csv(guide_filename, comment='#')
selected = df.iloc[python_guide_IDs]

selected.to_csv(guide_filename, index=True)
print("Done!")
