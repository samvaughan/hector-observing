import pandas as pd
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename")
args = parser.parse_args()
filename = args.filename

df = pd.read_csv(filename)

# We're going to check for probes which are near to a radial direction and then make sure that their angles are exactly radial
# This margin for error value sets how close to the correct angle a probe can be
margin_for_error = 10

print("Fixing the angles to be exactly radial...")

radial_angles = np.degrees(np.arctan(df.MagnetY / df.MagnetX))

n_changes = 0
for ((index, row), radial_angle) in zip(df.loc[np.isfinite(df['type'])].iterrows(), radial_angles):

    angs = row.angs
    angs_degrees = np.degrees(angs)

    if np.abs(angs_degrees - radial_angle) < 10:
        print(f"Fixing Hexabundle {index} to be radial: {angs_degrees:.3f} --> {radial_angle:.3f}")

        azAngs = row.azAngs
        angs_azAng = np.radians(radial_angle) - azAngs

        if angs_azAng > np.pi:
            angs_azAng = angs_azAng - 2 * np.pi
        elif angs_azAng < 0:
            angs_azAng = angs_azAng + 2 * np.pi

        df.at[index, 'angs'] = np.radians(radial_angle)
        df.at[index, 'angs_azAng'] = angs_azAng

    # Now check if we're aligned radially but pointing away
    if np.abs(angs_degrees - (radial_angle + 180)) < 10:
        print(f"Fixing Hexabundle {index} to be radial: {angs_degrees:.3f} --> {radial_angle + 180:.3f}")

        azAngs = row.azAngs
        angs_azAng = np.radians(radial_angle + 180) - azAngs

        if angs_azAng > np.pi:
            angs_azAng = angs_azAng - 2 * np.pi
        elif angs_azAng < 0:
            angs_azAng = angs_azAng + 2 * np.pi

        df.at[index, 'angs'] = np.radians(radial_angle + 180)
        df.at[index, 'angs_azAng'] = angs_azAng


print("Done!")
# Save the file
df.to_csv(filename)
