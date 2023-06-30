import matplotlib.pyplot as plt
import pandas as pd
import hop.hexabundle_allocation.hector.constants as constants
import misc_functions
import numpy as np
from matplotlib.patches import Circle, Rectangle
from matplotlib.transforms import Affine2D
import string

smk = snakemake # noqa

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
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect(1)
ax.axis("off")
ax.set_xlim(-350, 350)
ax.set_ylim(-350, 350)

ax.add_patch(plate_circle)
ax = misc_functions.draw_telecentricity_rings(ax)
ax.scatter(0, 0, marker="x", color="k")

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
        Hexabundle_label = Hexabundle[-1]

    # Label for the rectangular magnets is just the index
    rectangular_magnet_label = f"{cm.Index.values[0]:02}"

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
        hexabundle_text_colour = "yellow"

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
        color=hexabundle_text_colour,
        fontsize=11,
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
            if index in ["Sky-H3-8", "Sky-H4-8", "Sky-H5-8", "Sky-H6-7", "Sky-H7-8"]:
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

# Here's the invert y-axis which doesn't make sense to me...
ax.invert_yaxis()
print(f"Saving to {smk.output.final_plot}...")
fig.savefig(smk.output.final_plot, bbox_inches="tight")
print("Done!")
