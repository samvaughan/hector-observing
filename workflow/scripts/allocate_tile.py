from hop import pipeline

# import argparse
from pathlib import Path
import string

# parser = argparse.ArgumentParser()
# parser.add_argument('Guide_Fname')
# parser.add_argument('Hexas_fname')
# parser.add_argument('--field')

# args = parser.parse_args()
smk = snakemake  # noqa
field = smk.config["field"]
fname_Guides = smk.input.configured_guides_fname
fname_Hexas = smk.input.configured_tile_fname

# This is the output folder which I've made to put everything in in the 'Configuration folder'
output_folder = smk.config["output_folder"]
# assert output_folder in ['Evenly_Spread', 'Guides_Central', 'Hexas_Horizontal', 'Hexas_Vertical',
# 'Guides_Around_Edge', 'Hexas_Central', 'Hexas_Radial', 'Not_Configured', 'Guides_Central_Extra_Bright'], 'Unknown Output folder!'

robot_temperature = smk.config["robot_temp"]
obs_temperature = smk.config["obs_temp"]

fname_stem = Path(fname_Hexas).stem.replace("Hexas_", "")

base_allocation_folder = Path(f"results/TilingOutputs/{field}/Allocation/")
base_plot_folder = Path(f"results/TilingOutputs/{field}/Plots/")
base_final_files_folder = Path(f"results/TilingOutputs/{field}/FinalOutputs/")

output_folder_for_Allocation_folder = (
    base_allocation_folder / "tile_outputs" / f"{output_folder}"
)
output_folder_for_Plots_folder = base_plot_folder / f"{output_folder}"

# Make the folders
output_folder_for_Allocation_folder.mkdir(exist_ok=True)
output_folder_for_Plots_folder.mkdir(exist_ok=True)

config_filename = smk.config["pipeline_config_file"]
skyfiles_location = smk.config["skyfiles_location"]

try:
    P_and_Q_offset_filename = smk.config["P_and_Q_offset_filename"]
    print(
        f"Using a non-standard P and Q filename!! Currently {P_and_Q_offset_filename}"
    )
except KeyError:
    P_and_Q_offset_filename = "Hexa_final_prism_gluing_PQ_table.xlsx"

HP = pipeline.HectorPipe(
    config_filename=config_filename,
    Profit_files_location=skyfiles_location,
    P_and_Q_offset_filename=P_and_Q_offset_filename,
)

HP.allocate_hexabundles_for_single_tile(
    fileNameGuides=fname_Guides,
    fileNameHexa=fname_Hexas,
    robot_temperature=robot_temperature,
    obs_temperature=obs_temperature,
    plot=True,
)


print("Cleaning up the files...")
# Now move the output files into the output folder

# Move the Fibre_coordData_hexabundle_ files
for hexabundle in list(string.ascii_uppercase)[:21]:
    Path(
        base_allocation_folder
        / "tile_outputs"
        / f"Fibre_coordData_hexabundle_{hexabundle}.txt"
    ).rename(
        output_folder_for_Allocation_folder
        / f"Fibre_coordData_hexabundle_{hexabundle}.txt"
    )

# Move the HECTOROutput_Guides and HECTOROutput_Hexas files
Path(
    base_allocation_folder / "tile_outputs" / f"HECTOROutput_Guides_{fname_stem}.txt"
).rename(output_folder_for_Allocation_folder / f"HECTOROutput_Guides_{fname_stem}.txt")
Path(
    base_allocation_folder / "tile_outputs" / f"HECTOROutput_Hexas_{fname_stem}.txt"
).rename(output_folder_for_Allocation_folder / f"HECTOROutput_Hexas_{fname_stem}.txt")

# And now the slitInfo file
Path(
    base_allocation_folder / "tile_outputs" / f"Fibre_slitInfo_{fname_stem}.csv"
).rename(output_folder_for_Allocation_folder / f"Fibre_slitInfo_{fname_stem}.csv")

# And now the Plots

for plot_fname in [
    "fibre_slitletAAOmega",
    "fibre_slitletSpector",
    "hexabundlePlot",
    "robotPlot",
    "skyfibre_slitletAAOmega",
    "skyfibre_slitletSpector",
]:
    Path(base_plot_folder / f"{plot_fname}_{fname_stem}.pdf").rename(
        output_folder_for_Plots_folder / f"{plot_fname}_{fname_stem}.pdf"
    )

print("...Done!")
# Now make the final files for Tony
final_tile_file = (
    output_folder_for_Allocation_folder / f"HECTOROutput_Hexas_{fname_stem}.txt"
)
final_robot_file = base_allocation_folder / "robot_outputs" / f"Robot_{fname_stem}.txt"

print("Making the final files in the correct format...")
HP.make_output_file_for_Tony(final_tile_file, final_robot_file)

# Now move this file
Path(base_final_files_folder / f"{output_folder}").mkdir(exist_ok=True)
for prefix in ["Tile_FinalFormat", "Robot_FinalFormat"]:
    Path(base_final_files_folder / f"{prefix}_{fname_stem}.csv").rename(
        base_final_files_folder / f"{output_folder}" / f"{prefix}_{fname_stem}.csv"
    )

Path(smk.output.flag_file).touch()
print("Done!")
