import matplotlib.patches as patches
from hop.hexabundle_allocation.hector import constants
import numpy as np
from matplotlib.patches import Wedge


def draw_telecentricity_rings(ax):
    """Draw the concentric circles corresponding to our telecentricity boundaries

    Args:
        ax (_type_): _description_

    Returns:
        _type_: _description_
    """
    red = patches.Circle(
        xy=(0, 0),
        radius=constants.HECTOR_plate_radius,
        fill=True,
        color="#E78BE7",
        alpha=0.4,
    )
    yellow = patches.Circle(xy=(0, 0), radius=196.05124, fill=True, color="#f6f93b")
    green = patches.Circle(xy=(0, 0), radius=147.91658, fill=True, color="#60fb3d")
    blue = patches.Circle(xy=(0, 0), radius=92.71721, fill=True, color="#add8e6")

    ax.add_patch(red)
    ax.add_patch(yellow)
    ax.add_patch(green)
    ax.add_patch(blue)
    return ax


def wedge_patch(radius, theta_start, theta_end, alpha):
    skyfibre_wedge = Wedge(
        (0, 0),
        r=radius,
        theta1=theta_start,
        theta2=theta_end,
        width=90,
        facecolor="gray",
        edgecolor="black",
        alpha=alpha,
    )
    return skyfibre_wedge


def annotate_sky_fibre(ax, radius, angle, text, color="k", **kwargs):
    x = radius * np.cos(np.radians(angle))
    y = radius * np.sin(np.radians(angle))
    ax.annotate(
        text,
        (x, y),
        color=color,
        rotation=270 - angle,
        ha="center",
        va="center",
        **kwargs,
    )
    return ax
