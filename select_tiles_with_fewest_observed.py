import pandas as pd
import sqlite3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tiles_dataframe")

args = parser.parse_args()
tiles_df_filename = args.tiles_dataframe

con = con = sqlite3.connect(
    "/Users/samvaughan/Science/Hector/Database/HectorDB/hector.db"
)
observed_galaxies = pd.read_sql("select * from galaxies_observed", con)

tile_df = pd.read_csv(tiles_df_filename)


def n_already_observed(group):
    merged = pd.merge(
        group,
        observed_galaxies,
        left_on="ID",
        right_on="ID",
        suffixes=("", "_2"),
        how="inner",
    )

    return len(merged)


summary = tile_df.groupby("tile_number").agg(
    N=("ID", "size"),
    RA=("RA", "mean"),
    DEC=("DEC", "mean"),
    n_repeats=("ID", n_already_observed),
    MW_analogues=("MW_analogue", "sum"),
    EO_winds=("EO_wind_galaxy", "sum"),
)

summary = summary.loc[summary["N"] == 19]
print(summary.sort_values("EO_winds"))
