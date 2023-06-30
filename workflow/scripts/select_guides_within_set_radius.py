import pandas as pd
import numpy as np

smk = snakemake # noqa
df = pd.read_csv(smk.input.guide_file, comment='#', delim_whitespace=True)
outer_radius_for_selection = smk.params.outer_radius_for_selection

print(f"\nOnly keeping guide stars within a radius of {outer_radius_for_selection} degrees\n")

with open(smk.input.guide_file, 'r') as f:
    header = []
    for line in f:
        if line.startswith("#"):
            header.append(line)

field_centre_RA = float(header[1].split()[1])
field_centre_DEC = float(header[1].split()[2])

# Get the radius of each star from the centre, accounting for the cos(DEC) term
# This is in degrees
cos_dec_correction = np.cos(np.radians(np.abs(field_centre_DEC)))
radii = np.sqrt((df['RA'] - field_centre_RA)**2 * cos_dec_correction ** 2 + (df["DEC"] - field_centre_DEC)**2)

mask = radii < outer_radius_for_selection

# Now write the output file in the correct format
with open(smk.output.guide_file, 'a') as f:
    for line in header:
        f.write(line)
    df.loc[mask].to_csv(f, index=False, sep=' ')
