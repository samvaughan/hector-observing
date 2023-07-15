import string
import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u


def make_image_with_info(image, wcs, row):
    hexabundle_radius_dict = get_hexabundle_radius_dict()
    # Do the plotting
    fig, ax = plt.subplots(subplot_kw=dict(projection=wcs.celestial))
    ax.imshow(image, origin="lower")

    # Add the extra info
    ax.annotate(
        f"Hexabundle {row.Hexabundle}:\n{row.ID}\nMstar={row.Mstar:.2f}\nz={row.z:.4f}\nRe={row.Re:.2f}",
        xy=(0.03, 0.97),
        xycoords="axes fraction",
        bbox=dict(boxstyle="round", fc="0.4", ec="k", alpha=0.9),
        ha="left",
        va="top",
        color="w",
    )
    # Now add the bundle outline
    hexabundle_radius = hexabundle_radius_dict[row.Hexabundle]
    bundle_outline = SphericalCircle(
        center=(row.RA, row.DEC) * u.deg,
        radius=hexabundle_radius * u.arcsec,
        edgecolor="white",
        facecolor="none",
        linewidth=1.5,
        transform=ax.get_transform("fk5"),
    )
    ax.add_patch(bundle_outline)

    ax.axis("off")

    return fig, ax


def get_hexabundle_radius_dict():
    hexabundle_radius_dict = dict(
        A=13,
        B=13,
        C=11.2,
        D=19 / 2,
        E=15.5 / 2,
        F=15.5 / 2,
        G=15.5 / 2,
        H=12.1 / 2,
        I=19 / 2,
        J=15.5 / 2,
        K=15.5 / 2,
        L=15.5 / 2,
        M=15.5 / 2,
        N=15.5 / 2,
        O=15.5 / 2,
        P=15.5 / 2,
        Q=15.5 / 2,
        R=15.5 / 2,
        S=15.5 / 2,
        T=15.5 / 2,
        U=12.1 / 2,
        GS1=11,
        GS2=11,
        GS3=11,
        GS4=11,
        GS5=11,
        GS6=11,
    )
    return hexabundle_radius_dict


def get_spectrograph_from_hexabundle(row):
    if row.Hexabundle in (string.ascii_uppercase[:8]):
        return "AAOmega"
    elif row.Hexabundle in (string.ascii_uppercase[8:21]):
        return "Spector"
    elif row.Hexabundle.startswith("GS"):
        return "Guider"
    else:
        raise NameError(
            f"Don't know which spectrograph Hexabundle {row.Hexabundle} is in!"
        )
