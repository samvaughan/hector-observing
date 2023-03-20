from hop import pipeline
from hop.misc import misc_tools
from pathlib import Path

DC_tile_file = Path(snakemake.input.DC_tile_file)
configured_tile_file = Path(snakemake.input.configured_field)
output_tile_file = Path(snakemake.output.configured_tile_correct_header)

DC_guide_file = Path(snakemake.input.DC_guide_file)
configured_guide_file = Path(snakemake.input.configured_guide_field)
output_guide_file = Path(snakemake.output.configured_guide_file_correct_header)


field = snakemake.wildcards['field']
pipeline_config_dictionary = dict(output_folder=f"results/{snakemake.config['folder']}/TilingOutputs/{field}/", proximity=220, output_filename_stem='_')

skyfiles_location = Path(snakemake.config['skyfiles_location'])
if field.startswith('A'):
    print("\nUsing the Cluster sky masks\n")
    skyfiles_location = skyfiles_location / "Clusters"
elif field.startswith("G"):
    print("\nUsing the GAMA sky masks\n")
    skyfiles_location = skyfiles_location / "GAMA"
elif field.startswith("H"):
    print("\nUsing the WAVES sky masks\n")
    skyfiles_location = skyfiles_location / "WAVES"
elif field.startswith("custom_H"):
    print("\nUsing the WAVES sky masks\n")
    skyfiles_location = skyfiles_location / "WAVES"
elif field.startswith("S"):
    print("\nUsing the StarFields sky masks\n")
    skyfiles_location = skyfiles_location / "SkyMasks"
else:
    raise NameError("Not sure which part of the survey this field is from")
HP = pipeline.HectorPipe(config_dictionary=pipeline_config_dictionary, Profit_files_location=skyfiles_location)

# Get the correct header
tile_header = HP.read_header_dictionary_from_file(filename=DC_tile_file)
# Copy our other file to its new name
output_tile_file.write_text(configured_tile_file.read_text())
misc_tools.update_header(output_tile_file, tile_header)

guide_header = HP.read_header_dictionary_from_file(filename=DC_guide_file)
output_guide_file.write_text(configured_guide_file.read_text())
misc_tools.update_header(output_guide_file, guide_header)
