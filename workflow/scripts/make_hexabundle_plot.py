import matplotlib.pyplot as plt
import pandas as pd
import hop.hexabundle_allocation.hector.constants as constants
import misc_functions
import numpy as np
from matplotlib.patches import Circle, Rectangle
from matplotlib.transforms import Affine2D
import string
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


smk = snakemake  # noqa

# Load the Robot and Tile files
robot = pd.read_csv(
    smk.input.robot_file,
    skiprows=6,
)
tile = pd.read_csv(
    smk.input.tile_file,
    skiprows=11,
)

# Just select the skyfibres
skyfibres = tile.loc[tile["type"] < 0].set_index("ID")

plt.rc("font", size=30)  # controls default text sizes
plt.rc("axes", titlesize=30)  # fontsize of the axes title
plt.rc("axes", labelsize=30)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=30)  # fontsize of the tick labels
plt.rc("ytick", labelsize=30)  # fontsize of the tick labels
plt.rc("legend", fontsize=30)  # legend fontsize
plt.rc("figure", titlesize=30)  # fontsize of the figure title

# Some constants
rm_length = constants.rectangle_magnet_length * 1.7
rm_breadth = constants.rectangle_magnet_width * 1.5
circular_magnet_radius = constants.circular_magnet_radius * 2
robot_centre_in_mm = [constants.robot_center_x, constants.robot_center_y]
plate_radius = constants.HECTOR_plate_radius


# Find the distance to an exit
def get_exit_distance(x, y, exit_x, exit_y):
    return np.sqrt((x - exit_x) ** 2 + (y - exit_y) ** 2)


# Get the exit gaps of each bundle
def get_exit(row, angle_buffer=15, distance_to_exit_cutoff=200, radius_cutoff=180):
    """Try and find the appropriate exit for each bundle. This is surprisingly hard to get 100% right!
    This function mainly uses the angle of the hexabundle tail. We divide up the plate into thirds, with
    divisions at angles of 90, 210 and 330 degrees. A probe with an angle in the first quadrant goes to
    exit 1, etc. This is right ~90% of the time.
    However there are edge cases where a magnet is right next to an exit but its angle suggest it's actually
    for another one. These bundles have angles very close to our boundaries. If a probe has an angle within
    'angle_buffer' degrees of a boundary, and it's also within 'distance_to_exit_cutoff' mm from an exit,
    we instead pick that exit instead. E.g. if a probe has an angle of 82 degrees, it would normally go to exit
    1. However if it's also less than 100 mm from exit 2 then it will go to exit 2 instead.
    """

    boundary_points = np.array([90, 210, 330])

    x = row.Center_y - constants.robot_center_x
    y = row.Center_x - constants.robot_center_y

    # Get the distances to each exit.
    # Treat each exit as a point
    # This is the same location we plot their numbers later
    e1_distance = get_exit_distance(x, y, exit_x=250, exit_y=-150)
    e2_distance = get_exit_distance(x, y, exit_x=-250, exit_y=-150)
    e3_distance = get_exit_distance(x, y, exit_x=0, exit_y=280)
    distances = np.array([e1_distance, e2_distance, e3_distance])

    bundle_radius = np.sqrt(x**2 + y**2)

    # Select by angle
    if (0 < row.tail_angle < 90) | (330 < row.tail_angle < 360):
        exit = 1
    elif 90 < row.tail_angle < 210:
        exit = 2
    elif 210 < row.tail_angle < 330:
        exit = 3

    # Catch edge cases which have angles close to our boundary

    if (
        (np.any(np.abs((row.tail_angle - boundary_points)) < angle_buffer))
        & np.any(distances < distance_to_exit_cutoff)
        & (bundle_radius > radius_cutoff)
    ):
        # Choose the closest exit (and add 1 for 0-based indexing)
        exit = np.argmin([e1_distance, e2_distance, e3_distance]) + 1

    print(
        f"{row.Hexabundle}: exit {exit}, angle {row.tail_angle:.3f}, distance {distances}"
    )
    return exit


def exit_colour(row):
    if np.squeeze(row.exit_gap.values) == 1:
        return "#E4502D"
    elif np.squeeze(row.exit_gap.values) == 2:
        return "#4F86F6"
    elif np.squeeze(row.exit_gap.values) == 3:
        return "#C45FF6"


# Get the rectangular magnets only
rectangular_magnets = robot.loc[robot["#Magnet"] == "rectangular_magnet"].copy()
# Find their tail angles
tail_angles = (
    540 - rectangular_magnets.rot_holdingPosition - rectangular_magnets.rot_platePlacing
)
rectangular_magnets["tail_angle"] = tail_angles
exit_gaps = rectangular_magnets.apply(get_exit, axis=1)

# repeat to be the right length (so we can join back to the robot file)
tail_angles = pd.concat((tail_angles, tail_angles)).reset_index(drop=True)
exit_gaps = pd.concat((exit_gaps, exit_gaps)).reset_index(drop=True)
robot["tail_angle"] = tail_angles
robot["exit_gap"] = exit_gaps


def sort_distance(row):
    # Get the coordinate we need for the sorting
    # This is the y value for exit gaps 1 and 2, and the x axis for exit 3
    # Note the flipped coordinates here (i.e. x --> row.Center_y)
    if row.exit_gap == 3:
        return -(row.Center_y - constants.robot_center_y)
    else:
        return row.Center_x - constants.robot_center_x


robot["sort_distance"] = robot.apply(sort_distance, axis=1)

# Work out the order to plug these
guides = robot.loc[
    (robot["#Magnet"] == "circular_magnet") & (robot["Hexabundle"].str.startswith("GS"))
]
bundles = robot.loc[
    (robot["#Magnet"] == "circular_magnet")
    & ~(robot["Hexabundle"].str.startswith("GS"))
]

sorted_guides = guides.sort_values(["exit_gap", "sort_distance"])
sorted_bundles = bundles.sort_values(["exit_gap", "sort_distance"])
all_sorted = pd.concat((sorted_guides, sorted_bundles)).reset_index(drop=True)

# radii of various things for the skyfibre wedges
outer_wedge_radius = 333
title_radius = 270
bullet_radius = 308
skyfibre_position_radius = 322
skyfibre_number_radius = 342

# Get the list of Hexabundles
Hexabundles = np.unique(robot.Hexabundle)

# Add the plate
plate_circle = Circle(
    xy=(0, 0),
    radius=plate_radius,
    facecolor="None",
    edgecolor="k",
    linewidth=3.0,
    zorder=10,
)

# Make a plot
fig, ax = plt.subplots(figsize=(12, 12))
ax.set_aspect(1)
ax.axis("off")
ax.set_xlim(-350, 350)
ax.set_ylim(-350, 350)

ax.add_patch(plate_circle)
ax = misc_functions.draw_telecentricity_rings(ax)
ax.scatter(0, 0, marker="x", color="k")

# Add the exit gap labels
ax.annotate(
    "1",
    xy=(250, -150),
    ha="center",
    va="center",
    color="#E4502D",
    bbox=dict(boxstyle="round", alpha=0.3, facecolor="0.9"),
)
ax.annotate(
    "2",
    xy=(-250, -150),
    ha="center",
    va="center",
    color="#4F86F6",
    bbox=dict(boxstyle="round", alpha=0.3, facecolor="0.9"),
)
ax.annotate(
    "3",
    xy=(0, 280),
    ha="center",
    va="center",
    color="#C45FF6",
    bbox=dict(boxstyle="round", alpha=0.3, facecolor="0.9"),
)

# Loop through the hexabundles
for i, Hexabundle in enumerate(Hexabundles):
    # For each hexabundle, get its circular and rectangular magnet
    cm = robot.loc[
        (robot["Hexabundle"] == Hexabundle) & (robot["#Magnet"] == "circular_magnet")
    ]
    rm = robot.loc[
        (robot["Hexabundle"] == Hexabundle) & (robot["#Magnet"] == "rectangular_magnet")
    ]

    # Get the centres of each magnet
    # Note that as you look at the plate, the robot axes are: poitive x is TOWARDS you and positive y is to the RIGHT
    # To plot this such that the plot is the same as the configured field, we flip the x and y axes.
    # Note we also swap the trig functions below- so sin goes to cos an vice versa
    # This makes this code confusing and there's probably a better way to do this.
    # We also need to swap the Y axis at the end for things to look right- but if I add some minus signs *here* instead
    # Then everything doesn't work and the rectangular magnets get plotted in the wrong positions... Weird!
    # Anyway, this works so I'm going to leave it
    cm_x = cm.Center_y - constants.robot_center_y
    cm_y = cm.Center_x - constants.robot_center_x
    rm_x = rm.Center_y - constants.robot_center_y
    rm_y = rm.Center_x - constants.robot_center_x

    # Make the label for each bundle
    # This is the letter for the galaxy and standard star bundles
    # And the number for the guide bundles
    if len(Hexabundle) == 1:
        Hexabundle_label = Hexabundle
    else:
        Hexabundle_label = f"G{Hexabundle[-1]}"

    # Label for the rectangular magnets is just the index
    # rectangular_magnet_label = f"{cm.Index.values[0]:02}"
    rectangular_magnet_label = f"{np.squeeze(rm.exit_gap)}"

    # Colour AAOmega bundles white, Spector and guide bundles black
    if Hexabundle in string.ascii_uppercase[:8]:
        hexabundle_colour = "white"
        hexabundle_text_colour = "black"
    else:
        hexabundle_colour = "black"
        hexabundle_text_colour = "white"
    # Now check to see if we're looking at a guide bundle
    # If so, change the text colour to yellow
    if len(Hexabundle) > 1:
        hexabundle_text_colour = "gold"

    # The angle of the rectangular magnet- 270 minus the robot holding angle minus the robot placing angle
    rm_angle = np.squeeze(
        np.radians(270 - rm.rot_holdingPosition - rm.rot_platePlacing).values
    )

    # Plot the circular magnet
    c = Circle(
        xy=(cm_x, cm_y),
        radius=circular_magnet_radius,
        facecolor=hexabundle_colour,
        edgecolor="k",
        alpha=1,
        zorder=7,
    )
    ax.add_patch(c)

    # Plot the rectangular magnet
    r = Rectangle(
        xy=(rm_x - rm_length / 2, rm_y - rm_breadth / 2),
        width=rm_length,
        height=rm_breadth,
        transform=Affine2D().rotate_deg_around(*(rm_x, rm_y), 90 - np.degrees(rm_angle))
        + ax.transData,
        facecolor=hexabundle_colour,
        edgecolor="k",
        alpha=1,
        zorder=6,
    )
    ax.add_patch(r)

    # Add the hexabundle name, rotated to match the angle of the hexabundle
    text_rotation_angle = ((90 + np.degrees(rm_angle))) % 90

    # Add the label in the centre of the circular magnet
    ax.annotate(
        xy=(cm_x, cm_y),
        text=f"{Hexabundle_label}",
        xytext=(0, 0),
        textcoords="offset points",
        weight="bold",
        color=hexabundle_text_colour,
        fontsize=11,
        rotation=text_rotation_angle,
        ha="center",
        va="center",
        zorder=10,
    )

    # And do the same for the label in the centre of the rectangular magnet
    ax.annotate(
        xy=(rm_x, rm_y),
        text=rectangular_magnet_label,
        xytext=(0, 0),
        textcoords="offset points",
        weight="bold",
        color=exit_colour(rm),
        fontsize=10,
        rotation=text_rotation_angle,
        ha="center",
        va="center",
        zorder=10,
    )


# Now we need to draw the skyfibre wedges around the edge
skyfibreTitles_top = ["H3", "A3", "H4", "A4"]
skyfibreTitles_left = ["A1", "H1", "H2", "A2"]
skyfibreTitles_right = ["H5", "H6", "A5", "H7"]

# Loop through the three blocks
for angle, titles in zip(
    [0, 120, 240], [skyfibreTitles_right, skyfibreTitles_left, skyfibreTitles_top]
):
    # Loop through each wedge within a block
    for i, title in enumerate(titles):
        # Adjust the colours to match the original plots
        if title.startswith("H"):
            alpha = 0.4
        else:
            alpha = 0.7

        # Get the angles of the overall wedge
        # The wedges are 18 degrees in angular size
        theta_start = (angle + i * 20) - 9
        theta_end = (angle + i * 20) + 9
        central_angle = np.mean([theta_start, theta_end])

        # Plot the overall wedge
        skyfibre_wedge = misc_functions.wedge_patch(
            outer_wedge_radius, theta_start, theta_end, alpha
        )
        ax.add_artist(skyfibre_wedge)

        # Plot the title
        ax = misc_functions.annotate_sky_fibre(
            ax, title_radius, central_angle, title, fontsize=14, weight="bold"
        )

        # Find the fibres which belong to this wedge
        wedge_fibres = skyfibres.filter(regex=f"Sky-{title}-*", axis=0)

        # Now loop through individual fibres within a wedge
        for j, ((index, skyfibre), delta_angle) in enumerate(
            zip(wedge_fibres.iterrows(), [-7, -5, -3, -1, 1, 3, 5, 7])
        ):
            # Get the right angle for the text/symbols
            bullet_angle = central_angle + delta_angle

            # Check to see if it's one of the blocked off fibres
            colour = "k"
            skyfibre_position = skyfibre.SkyPosition
            if index in [
                "Sky-H3-8",
                "Sky-H4-8",
                "Sky-H5-8",
                "Sky-H6-7",
                "Sky-H7-8",
                "Sky-A2-2",
            ]:
                colour = "red"
                skyfibre_position = 0

            # Annotate the bullet symbol
            ax = misc_functions.annotate_sky_fibre(
                ax, bullet_radius, bullet_angle, "▮", fontsize=7, color=colour
            )
            # Annotate the Skyfibre position (0, 1, 2, or 3)
            ax = misc_functions.annotate_sky_fibre(
                ax,
                skyfibre_position_radius,
                bullet_angle,
                f"{skyfibre_position}",
                fontsize=9,
                weight="bold",
                color=colour,
            )
            # Annotate the fibre number
            ax = misc_functions.annotate_sky_fibre(
                ax,
                skyfibre_number_radius,
                bullet_angle,
                f"{index[-1]}",
                fontsize=9,
            )

info_string = r"$\it Suggested$ plugging order. Use your best judgement!"
guide_string = r"$\bf Guides$"
e1_string = r"$\bf Exit~1$"
e2_string = r"$\bf Exit~2$"
e3_string = r"$\bf Exit~3$"

order_string = f"""
{info_string}\n
{guide_string}
{" --> ".join(np.squeeze(sorted_guides['Hexabundle'].values))}\n
{e1_string}
{" --> ".join(np.atleast_1d(np.squeeze(sorted_bundles.loc[sorted_bundles.exit_gap == 1,'Hexabundle'].values)))}\n
{e2_string}
{" --> ".join(np.atleast_1d(np.squeeze(sorted_bundles.loc[sorted_bundles.exit_gap == 2,'Hexabundle'].values)))}\n
{e3_string}
{" --> ".join(np.atleast_1d(np.squeeze(sorted_bundles.loc[sorted_bundles.exit_gap == 3,'Hexabundle'].values)))}\n
"""
ax.annotate(
    order_string,
    xy=(1.05, 0.5),
    fontsize=14,
    ha="left",
    va="center",
    xycoords="axes fraction",
    bbox=dict(boxstyle="round", alpha=0.3, facecolor="0.9"),
)

# Here's the invert y-axis which doesn't make sense to me...
ax.invert_yaxis()
print(f"Saving to {smk.output.final_plot}...")
fig.savefig(smk.output.final_plot, bbox_inches="tight")
print("Done!")
