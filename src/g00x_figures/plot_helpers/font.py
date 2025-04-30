import matplotlib as mpl
import matplotlib.pyplot as plt


def apply_global_font_settings(fontsize: int = 16):
    mpl.rcParams["font.sans-serif"] = [
        "Arial",
        "Liberation Sans",
        "Bitstream Vera Sans",
    ]
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.size"] = fontsize
    mpl.rcParams["axes.titlesize"] = fontsize
    mpl.rcParams["axes.labelsize"] = fontsize
    mpl.rcParams["xtick.labelsize"] = fontsize
    mpl.rcParams["ytick.labelsize"] = fontsize
    mpl.rcParams["legend.fontsize"] = fontsize
    # table font size
    mpl.rcParams["figure.titlesize"] = fontsize
