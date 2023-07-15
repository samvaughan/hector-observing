import pandas as pd
import sqlite3
import misc
from pathlib import Path

smk = snakemake # noqa

con = sqlite3.connect(smk.input.database)
cur = con.cursor()

region = smk.wildcards["field"]
tile_number = int(smk.wildcards["tile_number"])

# Make the tile ID
tile_ID = f"{region}_{tile_number:05}"

print(f"Adding tile {tile_ID} to the configured tiles database...")

# Load the tile
df_tile = pd.read_csv(smk.input.tile_file, skiprows=11)
df_tile = df_tile.loc[df_tile["type"] >= 0]
df_tile["Spectrograph"] = df_tile.apply(misc.get_spectrograph_from_hexabundle, axis=1)
df_tile["tile_ID"] = tile_ID

# Add the letters to the start of the IDs if they're not there already
# Check if we need to first- see if it will be coerced to a numeric column or not
add_prefixes = True
try:
    pd.to_numeric(df_tile["ID"])
except ValueError:
    add_prefixes = False

if add_prefixes:
    # Guide stars
    df_tile.loc[df_tile["type"] == 2, "ID"] = df_tile.loc[df_tile["type"] == 2].apply(
        lambda x: f"G{x.ID}", axis=1
    )
    # Standard Stars
    df_tile.loc[df_tile["type"] == 0, "ID"] = df_tile.loc[df_tile["type"] == 0].apply(
        lambda x: f"S{x.ID}", axis=1
    )
    # Galaxies
    if region in ["G12", "G15", "G23", "H01", "H03"]:
        prefix = "W"
    elif region.startswith("A"):
        prefix = "C"
    elif region.startswith("Commissioning"):
        if region.split("_")[1].startswith("A"):
            prefix = "C"
        else:
            prefix = "W"
    else:
        prefix = "G"

    df_tile.loc[df_tile["type"] == 1, "ID"] = df_tile.loc[df_tile["type"] == 1].apply(
        lambda x: f"{prefix}{x.ID}", axis=1
    )

# Now add the row IDs
df_tile["row_ID"] = df_tile.apply(lambda row: f"{row.ID}_{row.tile_ID}", axis=1)

# Now we'll check if we're re-making a previous tile or adding a new one.
# If we're remaking a previous tile, check to make sure that it is an old tile just with some edits (i.e. a swapped Hexabundle)
# If it's a tile with entirely new galaxies but the same tile ID then something has gone wrong and we should raise an error!

previously_configured_tiles = pd.read_sql("SELECT * FROM configured_tiles", con)
old_tile = previously_configured_tiles.loc[
    previously_configured_tiles.tile_ID == tile_ID
]

# Delete any rows in the database which already have this tile ID
if len(old_tile) > 0:
    print("\tFound a previous tile with the same tile ID")
    assert set(df_tile.loc[df_tile['type'] == 1, 'row_ID']) == set(
        old_tile.loc[old_tile['type'] == 1, 'row_ID']
    ), "Looks like the previous tile in the database tile with matching tile ID has different galaxies to this tile! This shouldn't happen- you'll need to sort it out!"
    print("\tDeleting the old tile from the database...")

    cur.executescript(
        f"""
        DELETE FROM configured_tiles
        WHERE configured_tiles.tile_ID='{tile_ID}';
        """
    )
    print("\tDone")

configured_tile = (
    df_tile.loc[
        :,
        ["row_ID", "ID", "tile_ID", "Hexabundle", "Spectrograph", "RA", "DEC", "type"],
    ]
    .sort_values("Hexabundle")
    .reset_index(drop=True)
)

# Save it to the database
configured_tile.to_sql("configured_tiles", con, index=False, if_exists="append")
Path(smk.output.database_added_flag).touch()
print("Done!")
