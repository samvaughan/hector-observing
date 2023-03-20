import pandas as pd
import numpy as np
import argparse

"""
Code to edit the Tile and Robot files to swap a pair of hexabundles
"""

parser = argparse.ArgumentParser()
parser.add_argument("Run_date")
parser.add_argument("Field")
parser.add_argument("Tile_Number")
parser.add_argument("Hexabundle1")
parser.add_argument("Hexabundle2")

args = parser.parse_args()
master_folder = args.Run_date
field = args.field
tile_number = args.tile_number

tile_file = f'results/{master_folder}/TilingOutputs/{field}/FinalOutputs/Tile_FinalFormat_{field}_tile_{tile_number}_CONFIGURED_correct_header.csv'
robot_file = f'results/{master_folder}/TilingOutputs/{field}/FinalOutputs/Robot_FinalFormat_{field}_tile_{tile_number}_CONFIGURED_correct_header.csv'


# Swap the tile file
tile_df = pd.read_csv(tile_file)
tile_header = []
with open(tile_file, 'r') as f:
    if f.readline
    tile_header 