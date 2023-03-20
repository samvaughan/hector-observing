
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/opt/anaconda3/envs/snakemake/lib/python3.9/site-packages', '/Users/samvaughan/Science/Hector/Observing/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xdc\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c)resources/20221019/H01/Tiles/tile_001.fld\x94\x8c/resources/20221019/H01/Tiles/guide_tile_001.fld\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\ttile_file\x94K\x00N\x86\x94\x8c\nguide_file\x94K\x01N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x15\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1b)}\x94\x8c\x05_name\x94h\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh\x0fh\nh\x11h\x0bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8cJresults/20221019/TilingOutputs/H01/DistortionCorrected/DC_H01_tile_001.csv\x94\x8cPresults/20221019/TilingOutputs/H01/DistortionCorrected/guide_DC_H01_tile_001.csv\x94e}\x94(h\r}\x94(\x8c\x0coutput_fname\x94K\x00N\x86\x94\x8c\x12output_guide_fname\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh-h)h/h*ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94\x8c\x05{H01}\x94a}\x94(h\r}\x94\x8c\x05field\x94K\x00N\x86\x94sh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhAh>ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x03H01\x94\x8c\x03001\x94e}\x94(h\r}\x94(\x8c\x05field\x94K\x00N\x86\x94\x8c\x0btile_number\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhAhP\x8c\x0btile_number\x94hQub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/bd/ym7_px1s7_53zt92h0rfdb000000gp/T\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhjK\x01hlK\x01hnhgub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x06folder\x94J[\x8c4\x01\x8c\x0btiles_table\x94\x8c%resources/20221019/20221019_tiles.csv\x94\x8c\x11skyfiles_location\x94\x8c2/Volumes/OtherFiles/Science/Hector/SegMaps_Hector/\x94\x8c\x08obs_date\x94\x8c\n2022 10 27\x94\x8c\nrobot_temp\x94K\x0e\x8c\x08obs_temp\x94K\x0e\x8c\x16Hector_Linear_location\x94\x8c\x97/Users/samvaughan/Science/Hector/HectorObservationPipeline/hop/distortion_correction/HectorTranslationSoftware/DataFiles/April2022_Pos/HectorLinear.sds\x94\x8c\x1aHector_Distortion_Location\x94\x8c\x9b/Users/samvaughan/Science/Hector/HectorObservationPipeline/hop/distortion_correction/HectorTranslationSoftware/DataFiles/April2022_Pos/HectorDistortion.sds\x94u\x8c\x04rule\x94\x8c\x1erun_distortion_correction_code\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c;/Users/samvaughan/Science/Hector/Observing/workflow/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/samvaughan/Science/Hector/Observing/workflow/scripts/apply_DC_correction.py';
######## snakemake preamble end #########
from hop import pipeline
import pandas as pd
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
from hop.misc import misc_tools
import astroplan
import astropy.units as u
from astroplan import Observer
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimezoneInfo
from dateutil.parser import parse

# import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument('tile_fname')
# parser.add_argument('guide_fname')
# parser.add_argument('--field')
# args = parser.parse_args()

master_folder = snakemake.config['folder']
field = snakemake.wildcards['field']
tile_number = snakemake.wildcards['tile_number']

pipeline_config_dictionary = dict(output_folder=f"results/{master_folder}/TilingOutputs/{field}/", proximity=220, output_filename_stem='_')

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

tile_fname = Path(snakemake.input.tile_file)
guide_fname = Path(snakemake.input.guide_file)

output_fname = snakemake.output.output_fname
output_guide_fname = snakemake.output.output_guide_fname

output_fname_for_config = f"{pipeline_config_dictionary['output_folder']}/Configuration/Hexas_tile_{field}_{tile_number}_NOT_CONFIGURED.csv"
output_guide_fname_for_config = f"{pipeline_config_dictionary['output_folder']}/Configuration/Guides_tile_{field}_{tile_number}_NOT_CONFIGURED.csv"
plot_save_filename = f"{pipeline_config_dictionary['output_folder']}/DistortionCorrected/DC_{tile_fname.stem}.pdf"


# Get the centre of the tile
with open(tile_fname, 'r') as f:
    lines = f.readlines()

tile_RA, tile_DEC = lines[1].strip().lstrip('#').split()


guide_stars_for_tile = pd.read_csv(guide_fname, comment='#', delim_whitespace=True)

# Configure the field for the time when it is at its highest point
date = snakemake.config['obs_date']

# Get a SkyCoord of the centre of the tile
target = FixedTarget(name=f'{field}', coord=SkyCoord(ra=float(tile_RA) * u.degree, dec=float(tile_DEC)*u.degree))

# Get the time corresponding to the middle of the night in australia
australia_UTC = TimezoneInfo(utc_offset=10*u.hour)
time = parse(date) # Not giving a time makes it midnight
# Now move this time to the Australian time zone by applying the +10hours UT offset
time = time.replace(tzinfo=australia_UTC)
# Now get this in UTC time by using the Astropy Time package (see https://docs.astropy.org/en/stable/time/index.html#timezones; Time(datetimeobject) makes it lose its timezone properties)
utc_time = Time(time)

# Now find the closest transit on this date at the AAT
AAT = Observer.at_site("sso")
meridian_time = AAT.target_meridian_transit_time(utc_time, target)
meridian_time_formatted = meridian_time.strftime("%Y %m %d %H %M %S")

#Do some basic checks- that this target is up and that it's actually night time!
#assert AAT.is_night(utc_time) & AAT.target_is_up(time, target), f"A time of {meridian_time_formatted} fails AAT.is_night() and AA.target_is_up() !!!"

print(f"Configuring this field for a UT time of {meridian_time_formatted}")
current_date = datetime.datetime.now().strftime("%Y %m %d")

robot_temperature = snakemake.config['robot_temp']
obs_temperature = snakemake.config['obs_temp']


# Set up the correct Hector distortion and Hector linearity files
HP.Hector_distortion_file_location= Path(snakemake.config['Hector_Distortion_Location']) #Path("/Users/samvaughan/Science/Hector/HectorObservationPipeline/hop/distortion_correction/HectorTranslationSoftware/DataFiles/HectorDistortion.sds") #Path(snakemake.config['Hector_Distortion_Location'])
HP.Hector_linearity_file_location = Path(snakemake.config['Hector_Linear_location']) #Path("/Users/samvaughan/Science/Hector/HectorObservationPipeline/hop/distortion_correction/HectorTranslationSoftware/DataFiles/HectorLinear.sds") #Path(snakemake.config['Hector_Linear_location'])

if HP.Hector_linearity_file_location.parent.name == "fit2-b-pos":
    correction_direction = "POSITIVE"
elif HP.Hector_linearity_file_location.parent.name == "fit2-b-neg":
    correction_direction = "NEGATIVE"
elif HP.Hector_linearity_file_location.parent.name == 'PreFeb2022':
    correction_direction = "Pre_Feb2022"
elif HP.Hector_linearity_file_location.parent.name == 'April2022_Pos':
    correction_direction = "April 2022 POSITIVE"
else:
    correction_direction = "April2022"

print(f"\nOptical Model Correction is {correction_direction}\n")
print(f"Tdf linearity file is {HP.Hector_linearity_file_location}, Optical Model is {HP.Hector_distortion_file_location}\n")

print("Running the DC code...")
HP.apply_distortion_correction(tile_fname, guide_fname, output_fname, output_guide_fname, tile_RA=tile_RA, tile_Dec=tile_DEC, guide_stars_for_tile=guide_stars_for_tile, plot_save_filename=plot_save_filename, date=meridian_time_formatted, robot_temp=robot_temperature, obs_temp=obs_temperature, label=f"test", plateID=1, distortion_file=HP.Hector_distortion_file_location, linearity_file=HP.Hector_linearity_file_location, sky_fibre_file=HP.Hector_sky_fibre_location, profit_file_dir=HP.Profit_files_location, check_sky_fibres=True, verbose=False, apply_PM_corr=True)

plt.close('all')
print("...Done!")

# Now make the file which we can configure manually if the automatic version fails

# Get the header so we can add the lines we expect
tile_header, tile_header_numbers = HP._read_file_get_header(output_fname)
guide_header, guide_header_numbers = HP._read_file_get_header(output_guide_fname)

# Add some new rows to the header
tile_header = ["#PROXIMITY,220 # tiling proximity value in arcseconds\n"] + [f"#TILING_DATE,{current_date} # Date the tile was created/configured\n"] + tile_header


guide_header = ["#PROXIMITY,220 # tiling proximity value in arcseconds\n"] + [f"#TILING_DATE,{current_date} # Date the tile was created/configured\n"] + guide_header



# read the files in again
tile_df = pd.read_csv(output_fname, comment='#')
guide_df = pd.read_csv(output_guide_fname, comment='#')



def save_file_with_header(outfilename, header, dataframe):

    # First delete the contents of the file
    with open(outfilename, 'w') as g:
        pass
    # Now write to it
    with open(outfilename, 'a') as f:

        for line in header:
            f.write(line)

        # Now write a comment character before the column names get written
        # Tony requires us to do this 
        dataframe.to_csv(f, index=False)

save_file_with_header(output_fname, tile_header, tile_df)

# Write the guide file
save_file_with_header(output_guide_fname, guide_header, guide_df)

# Now make the files in case we have to hand configure the tile
# We need to add in the extra columns for the config code: "x","y","rads","angs","azAngs","angs_azAng"
tile_df['x'] = tile_df['MagnetX'] / 1000.0
tile_df['y'] = tile_df['MagnetY'] / 1000.0
tile_df.loc[:, ["rads","angs","azAngs","angs_azAng"]] = 0.0


guide_df['x'] = guide_df['MagnetX'] / 1000.0
guide_df['y'] = guide_df['MagnetY'] / 1000.0
guide_df.loc[:, ["rads","angs","azAngs","angs_azAng"]] = 0.0

#save_file_with_header(output_fname_for_config, tile_header, tile_df)
print(f"Saving 'NOT_CONFIGURED' outputs to {output_fname_for_config} and {output_guide_fname_for_config}")
tile_df.to_csv(output_fname_for_config)
guide_df.to_csv(output_guide_fname_for_config)
# Write the guide file
#save_file_with_header(output_guide_fname_for_config, guide_header, guide_df)
