import pandas as pd
import numpy as np
from pathlib import Path
import string

# Things to check against
N_hexabundles = 19
N_guide_stars = 6
N_standard_stars = 2
N_lines_in_tile_file = 118

N_magnets = 27

required_tile_header_keywords = ['PROXIMITY',
 'TILING_DATE',
 'OBS_TEMP',
 'ROBOT_TEMP',
 'MDLPARS',
 'EQUINOX',
 'CENTRE',
 'UTTIME',
 'UTDATE',
 'PLATEID',
 'LABEL']

required_robot_header_keywords = ['Label',
'Date_and_Time_file_created',
'Radial_Offset_Adjustment',
'RobotTemp',
'ObsTemp',
'Radial_Scale_factor']

all_hexabundles = list(string.ascii_uppercase[:21]) + ['GS1', 'GS2', 'GS3', 'GS4', 'GS5', 'GS6'] 


tile_file = snakemake.input.tile_file
robot_file = snakemake.input.robot_file

df_tile = pd.read_csv(tile_file, skiprows=11)
df_robot = pd.read_csv(robot_file, skiprows=6)

# Check the type of each hexabundle which has been allocated
hexabundle_type = df_tile['type'].values

print("Checking tile file contents:")
print("\tChecking number of hexabundles...")
assert np.sum(hexabundle_type == 1) == N_hexabundles, f"The tile file doesn't have {N_hexabundles} hexabundles allocated! It only has {np.sum(hexabundle_type == 1)}"
print("\tChecking number of guide stars...")
assert np.sum(hexabundle_type == 2) == N_guide_stars, f"The tile file doesn't have {N_guide_stars} guide stars allocated! It only has {np.sum(hexabundle_type == 2)}"
print("\tChecking number of standard stars...")
assert np.sum(hexabundle_type == 0) == N_standard_stars, f"The tile file doesn't have {N_standard_stars} standard stars allocated! It only has {np.sum(hexabundle_type == 0)}"
print("\tChecking file length...")
assert len(df_tile) == N_lines_in_tile_file, f"The tile file should have {N_lines_in_tile_file} rows! It currently has {len(df_tile)}"

print("Passed!")


# Check the robot file has the right number of lines...
print("\nChecking robot file contents:")
print("\tChecking we have the correct number of circular and rectangular magnets...")
circular_magnets = df_robot.loc[df_robot['#Magnet'] == 'circular_magnet']
rectangular_magnets = df_robot.loc[df_robot['#Magnet'] == 'rectangular_magnet']

assert len(circular_magnets) == N_magnets, f"We should have {N_magnets} circular magnets but in fact we have {len(circular_magnets)}"
assert len(rectangular_magnets) == N_magnets, f"We should have {N_magnets} circular magnets but in fact we have {len(rectangular_magnets)}"

print("\tChecking file length...")
assert len(df_robot) == 2*N_magnets, f"The robot file should have {2*N_magnets} rows! It currently has {len(df_robot)}"
print("Passed!")

def check_header(header, required_header_keywords):

    keywords = []
    values = []
    for line in header:
        split = line.split(',', 1)
        keywords.append(split[0])
        values.append(split[1])

    assert set(keywords) == set(required_header_keywords), "The tile file header doesn't contain all the required keywords!"



print("\nChecking the tile file header properties:")
with open(tile_file, 'r') as f:
    tile_file_contents = f.readlines()

tile_header = tile_file_contents[:11]

print("\tChecking for the presence of required header keywords...")
check_header(tile_header, required_tile_header_keywords)
print("\tChecking the column headings line starts with #...")
assert tile_file_contents[11].startswith("#"), "The line with the column names must start with #!"
print("Passed!")


print("\nChecking the robot file header properties:")
with open(robot_file, 'r') as f:
    robot_file_contents = f.readlines()

robot_header = robot_file_contents[:6]


print("\tChecking for the presence of required header keywords...")
check_header(robot_header, required_robot_header_keywords)
print("\tChecking the column headings line starts with #...")
assert robot_file_contents[6].startswith("#"), "The line with the column names must start with #!"
print("Passed!")

print(f"\nMaking the flag file: {snakemake.output.verification_passed_file}")
Path(snakemake.output.verification_passed_file).touch()