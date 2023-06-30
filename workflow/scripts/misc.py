import string

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
