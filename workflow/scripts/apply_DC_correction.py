from hop import pipeline
import pandas as pd
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import astropy.units as u
from astroplan import Observer
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimezoneInfo
from dateutil.parser import parse

smk = snakemake # noqa
master_folder = smk.config["folder"]
field = smk.wildcards["field"]
tile_number = smk.wildcards["tile_number"]

label_string = smk.params.label_string

pipeline_config_dictionary = dict(
    output_folder=f"results/{master_folder}/TilingOutputs/{field}/",
    proximity=220,
    output_filename_stem="_",
)

skyfiles_location = Path(smk.config["skyfiles_location"])
if field.startswith("A"):
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
elif field.startswith("Commissioning"):
    print("\nUsing the GAMA sky masks\n")
    skyfiles_location = skyfiles_location / "GAMA"
else:
    raise NameError("Not sure which part of the survey this field is from")

HP = pipeline.HectorPipe(
    config_dictionary=pipeline_config_dictionary,
    Profit_files_location=skyfiles_location,
)

tile_fname = Path(smk.input.tile_file)
guide_fname = Path(smk.input.guide_file)

output_fname = smk.output.output_fname
output_guide_fname = smk.output.output_guide_fname

output_fname_for_config = f"{pipeline_config_dictionary['output_folder']}/Configuration/Hexas_{field}_tile_{tile_number}_NOT_CONFIGURED.csv"
output_guide_fname_for_config = f"{pipeline_config_dictionary['output_folder']}/Configuration/Guides_{field}_tile_{tile_number}_NOT_CONFIGURED.csv"
plot_save_filename = f"{pipeline_config_dictionary['output_folder']}/DistortionCorrected/DC_{tile_fname.stem}.pdf"


# Get the centre of the tile
with open(tile_fname, "r") as f:
    lines = f.readlines()

tile_RA, tile_DEC = lines[1].strip().lstrip("#").split()


guide_stars_for_tile = pd.read_csv(guide_fname, comment="#", delim_whitespace=True)

# Configure the field for the time when it is at its highest point
date = smk.config["obs_date"]

# Get a SkyCoord of the centre of the tile
target = FixedTarget(
    name=f"{field}",
    coord=SkyCoord(ra=float(tile_RA) * u.degree, dec=float(tile_DEC) * u.degree),
)

# Get the time corresponding to the middle of the night in australia
australia_UTC = TimezoneInfo(utc_offset=10 * u.hour)
time = parse(date)  # Not giving a time makes it midnight
# Now move this time to the Australian time zone by applying the +10hours UT offset
time = time.replace(tzinfo=australia_UTC)
# Now get this in UTC time by using the Astropy Time package (see https://docs.astropy.org/en/stable/time/index.html#timezones; Time(datetimeobject) makes it lose its timezone properties)
utc_time = Time(time)

# Now find the closest transit on this date at the AAT
AAT = Observer.at_site("sso")
meridian_time = AAT.target_meridian_transit_time(utc_time, target)
meridian_time_formatted = meridian_time.strftime("%Y %m %d %H %M %S")

# Do some basic checks- that this target is up and that it's actually night time!
# assert AAT.is_night(utc_time) & AAT.target_is_up(time, target), f"A time of {meridian_time_formatted} fails AAT.is_night() and AA.target_is_up() !!!"

print(f"Configuring this field for a UT time of {meridian_time_formatted}")
current_date = datetime.datetime.now().strftime("%Y %m %d")

robot_temperature = smk.config["robot_temp"]
obs_temperature = smk.config["obs_temp"]


# Set up the correct Hector distortion and Hector linearity files
HP.Hector_distortion_file_location = Path(
    smk.config["Hector_Distortion_Location"]
)  # Path("/Users/samvaughan/Science/Hector/HectorObservationPipeline/hop/distortion_correction/HectorTranslationSoftware/DataFiles/HectorDistortion.sds") #Path(smk.config['Hector_Distortion_Location'])
HP.Hector_linearity_file_location = Path(
    smk.config["Hector_Linear_location"]
)  # Path("/Users/samvaughan/Science/Hector/HectorObservationPipeline/hop/distortion_correction/HectorTranslationSoftware/DataFiles/HectorLinear.sds") #Path(smk.config['Hector_Linear_location'])

if HP.Hector_linearity_file_location.parent.name == "fit2-b-pos":
    correction_direction = "POSITIVE"
elif HP.Hector_linearity_file_location.parent.name == "fit2-b-neg":
    correction_direction = "NEGATIVE"
elif HP.Hector_linearity_file_location.parent.name == "PreFeb2022":
    correction_direction = "Pre_Feb2022"
elif HP.Hector_linearity_file_location.parent.name == "April2022_Pos":
    correction_direction = "April 2022 POSITIVE"
else:
    correction_direction = "April2022"

print(f"\nOptical Model Correction is {correction_direction}\n")
print(
    f"Tdf linearity file is {HP.Hector_linearity_file_location}, Optical Model is {HP.Hector_distortion_file_location}\n"
)

print("Running the DC code...")
HP.apply_distortion_correction(
    tile_fname,
    guide_fname,
    output_fname,
    output_guide_fname,
    tile_RA=tile_RA,
    tile_Dec=tile_DEC,
    guide_stars_for_tile=guide_stars_for_tile,
    plot_save_filename=plot_save_filename,
    date=meridian_time_formatted,
    robot_temp=robot_temperature,
    obs_temp=obs_temperature,
    label=1,
    plateID=label_string,
    distortion_file=HP.Hector_distortion_file_location,
    linearity_file=HP.Hector_linearity_file_location,
    sky_fibre_file=HP.Hector_sky_fibre_location,
    profit_file_dir=HP.Profit_files_location,
    check_sky_fibres=True,
    verbose=False,
    apply_PM_corr=True,
)

plt.close("all")
print("...Done!")

# Now make the file which we can configure manually if the automatic version fails

# Get the header so we can add the lines we expect
tile_header, tile_header_numbers = HP._read_file_get_header(output_fname)
guide_header, guide_header_numbers = HP._read_file_get_header(output_guide_fname)

# Add some new rows to the header
tile_header = (
    ["#PROXIMITY,220 # tiling proximity value in arcseconds\n"]
    + [f"#TILING_DATE,{current_date} # Date the tile was created/configured\n"]
    + tile_header
)


guide_header = (
    ["#PROXIMITY,220 # tiling proximity value in arcseconds\n"]
    + [f"#TILING_DATE,{current_date} # Date the tile was created/configured\n"]
    + guide_header
)

# read the files in again
tile_df = pd.read_csv(output_fname, comment="#")
guide_df = pd.read_csv(output_guide_fname, comment="#")


def save_file_with_header(outfilename, header, dataframe):
    # First delete the contents of the file
    with open(outfilename, "w"):
        pass
    # Now write to it
    with open(outfilename, "a") as f:
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
tile_df["x"] = tile_df["MagnetX"] / 1000.0
tile_df["y"] = tile_df["MagnetY"] / 1000.0
tile_df.loc[:, ["rads", "angs", "azAngs", "angs_azAng"]] = 0.0


guide_df["x"] = guide_df["MagnetX"] / 1000.0
guide_df["y"] = guide_df["MagnetY"] / 1000.0
guide_df.loc[:, ["rads", "angs", "azAngs", "angs_azAng"]] = 0.0

# save_file_with_header(output_fname_for_config, tile_header, tile_df)
print(
    f"Saving 'NOT_CONFIGURED' outputs to {output_fname_for_config} and {output_guide_fname_for_config}"
)
tile_df.to_csv(output_fname_for_config)
guide_df.to_csv(output_guide_fname_for_config)
# Write the guide file
# save_file_with_header(output_guide_fname_for_config, guide_header, guide_df)
