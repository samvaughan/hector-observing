import pandas as pd
import argparse
from pathlib import Path
import shutil

"""
Code to edit the Tile and Robot files to swap a pair of hexabundles
"""

# parser = argparse.ArgumentParser()
# parser.add_argument("Run_date")
# parser.add_argument("Field")
# parser.add_argument("Tile_Number")
# parser.add_argument("Hexabundle1")
# parser.add_argument("Hexabundle2")

# args = parser.parse_args()
# master_folder = args.Run_date
# field = args.field
# tile_number = args.tile_number
# hexabundle_1 = args.Hexabundle1
# hexabundle_2 = args.Hexabundle2


# tile_file = Path(
#     f"results/{master_folder}/TilingOutputs/{field}/FinalOutputs/Tile_FinalFormat_{field}_tile_{tile_number}_CONFIGURED_correct_header.csv"
# )
# robot_file = Path(
#     f"results/{master_folder}/TilingOutputs/{field}/FinalOutputs/Robot_FinalFormat_{field}_tile_{tile_number}_CONFIGURED_correct_header.csv"
# )
def pair(arg):
    # For simplity, assume arg is a pair of letters
    # separated by a dash. If you want to do more
    # validation, raise argparse.ArgumentError if you
    # encounter a problem.
    return [x for x in arg.split("-")]


parser = argparse.ArgumentParser(prog="Swap Hexabundles for a Hector tile")
parser.add_argument("tile_file", help="The Hector Tile filename")
parser.add_argument("robot_file", help="The Hector Robot filename")
parser.add_argument(
    "swaps_list",
    type=pair,
    nargs="+",
    help="A list of hexabundle pairs to swap. These are two letters (A through U) separated by a dash, e.g. A-B would swap the galaxies in Hexabundles A and B. You can list multiple pairs",
)
args = parser.parse_args()

tile_file = Path(args.tile_file)
robot_file = Path(args.robot_file)

swaps_list = args.swaps_list

n_elements_per_pair = [len(x) for x in swaps_list]
assert all(
    [n == 2 for n in n_elements_per_pair]
), "The swaps list must be made up of pairs of Hexabundles separated by a dash (-). Looks like some letters aren't in pairs. To do three-way swaps, enter multiple pairs. You can repeat a letter in sa many pairs as you like, e.g 'A-B B-C' is allowed, but A-B-C isn't."

print("\tMaking a backup of the files")
backup_tile_fname = tile_file.parent / (tile_file.name + ".backup")
backup_robot_fname = robot_file.parent / (robot_file.name + ".backup")

shutil.copy(tile_file, backup_tile_fname)
shutil.copy(robot_file, backup_robot_fname)


# Swap the hexabundles in the tile file
tile_df = pd.read_csv(tile_file, skiprows=11)
tile_header = []
with open(tile_file, "r") as f:
    lines = f.readlines()
    tile_header = lines[:11]


# Now do the Robot file
robot_df = pd.read_csv(robot_file, skiprows=6)
robot_header = []
with open(robot_file, "r") as f:
    lines = f.readlines()
    robot_header = lines[:6]

for hexabundle_1, hexabundle_2 in swaps_list:
    print(f"\tSwapping the Hexabundles {hexabundle_1} and {hexabundle_2}")
    hexabundle_1_mask = tile_df["Hexabundle"] == hexabundle_1
    hexabundle_2_mask = tile_df["Hexabundle"] == hexabundle_2
    tile_df.loc[hexabundle_1_mask, "Hexabundle"] = hexabundle_2
    tile_df.loc[hexabundle_2_mask, "Hexabundle"] = hexabundle_1

    hexabundle_1_mask = robot_df["Hexabundle"] == hexabundle_1
    hexabundle_2_mask = robot_df["Hexabundle"] == hexabundle_2
    robot_df.loc[hexabundle_1_mask, "Hexabundle"] = hexabundle_2
    robot_df.loc[hexabundle_2_mask, "Hexabundle"] = hexabundle_1

with open(tile_file, "w") as g:
    for line in tile_header:
        g.write(line)
    tile_df.to_csv(g, index=False)


with open(robot_file, "w") as g:
    for line in robot_header:
        g.write(line)
    robot_df.to_csv(g, index=False)

print("Done!")
