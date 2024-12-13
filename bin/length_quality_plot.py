#!/usr/bin/env python3
# Makes a nice length vs quality plot
# Diego Alvarez <dialvarezs@gmail.com>
# 2021-10-12
import argparse
import gzip
from dataclasses import dataclass
from typing import List, Optional, Tuple
import pandas as pd
import numpy as np

import seaborn as sns
from Bio import SeqIO
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt


@dataclass
class Args:
    input: str
    output: str
    min_length: int
    max_length: int
    min_quality: int
    max_quality: int
    figure_width: int
    figure_height: int
    bins_x: int
    bins_y: int


def main():
    args = parse_arguments()

    qualities, lenghts = load_quality_and_length(
        args.input,
        min_length=args.min_length,
        max_length=args.max_length,
        min_quality=args.min_quality,
    )

    p = plot(
        lenghts,
        qualities,
        bins_x=args.bins_x,
        bins_y=args.bins_y,
        xlimits=(args.min_length, args.max_length),
        ylimits=(args.min_quality, args.max_quality),
        xticks=np.arange(args.min_length, args.max_length + 1, 50),
        yticks=np.arange(args.min_quality, args.max_quality + 1, 3),
    )

    #set figure size and save
    p.fig.set_figwidth(args.figure_width)
    p.fig.set_figheight(args.figure_height)
    p.fig.savefig(args.output, bbox_inches="tight")


def plot(
    lenghts: np.ndarray,
    qualities: np.ndarray,
    bins_x: int,
    bins_y: int,
    xlimits=Optional[Tuple[int, int]],
    ylimits=Optional[Tuple[int, int]],
    xticks=Optional[List[int]],
    yticks=Optional[List[int]],
) -> sns.JointGrid:
    """
    Make the plot
    """
    palette = get_continuous_cmap(["#ffffff", "#5dba9c", "#2d6b57"])

    sns.set_context("paper", font_scale=1.5)
    with sns.axes_style(
        {
            "figure.facecolor": "white",
            "axes.labelcolor": "#2a3f5f",
            "axes.edgecolor": "#2a3f5f",
            "xtick.color": "#2a3f5f",
            "ytick.color": "#2a3f5f",
            "font.sans-serif": ["Lato"] + sns.axes_style()["font.sans-serif"],
        }
    ):
        g = sns.JointGrid(
            x=lenghts,
            y=qualities,
        )
        g.plot_joint(
            plt.hexbin, edgecolors=None, cmap=palette, gridsize=(bins_x, bins_y)
        )
        sns.histplot(
            x=lenghts,
            bins=bins_x + 1,
            color="#4cb391",
            lw=0,
            ax=g.ax_marg_x,
        )
        sns.histplot(
            y=qualities,
            bins=bins_y * 2,
            color="#4cb391",
            lw=0,
            ax=g.ax_marg_y,
        )

        if xlimits is not None:
            g.ax_joint.set_xlim(xlimits)
        if ylimits is not None:
            g.ax_joint.set_ylim(ylimits)
        if xticks is not None:
            g.ax_joint.set_xticks(xticks)
        if yticks is not None:
            g.ax_joint.set_yticks(yticks)

        g.ax_joint.spines.left.set_position(("axes", -0.012))
        g.ax_joint.spines.bottom.set_position(("axes", -0.03))
        g.ax_marg_x.axis("off")
        g.ax_marg_y.axis("off")

        g.ax_joint.set(xlabel="Largo de lecturas (pb)", ylabel="Calidad Phred promedio de lectura")

    return g


def load_quality_and_length(
    file_name: str,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    min_quality: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Loads quality scores and lengths from the reads on a FASTQ file.

    :param file_name: input file path
    :param min_length: minimum length of reads to include
    :param max_length: maximum length of reads to include
    :return: a tuple with two numpy arrays, one with the quality scores and one with the lengths
    """
    # handle gzipped files
    handle = (
        gzip.open(file_name, "rt")
        if file_name.endswith(".gz")
        else open(file_name, "r")
    )

    qualities = []
    lengths = []
    for record in SeqIO.parse(handle, "fastq"):
        # if read length is not within the specified range, skip it
        if min_length is not None and len(record.seq) < min_length:
            continue
        if max_length is not None and len(record.seq) > max_length:
            continue

        mean_quality = get_quality_mean(record.letter_annotations["phred_quality"])
        if min_quality is not None and mean_quality < min_quality:
            continue

        qualities.append(mean_quality)
        lengths.append(len(record))

    qualities = np.array(qualities)
    lenghts = np.array(lengths)

    return qualities, lenghts


def get_quality_mean(qualities: List[int]) -> int:
    """
    Gets the mean quality of a list of quality scores considering the logarithmic scale
    """
    return -10 * np.log10(np.mean(10 ** (-1 * np.array(qualities, dtype=int) / 10)))


# source: https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
def get_continuous_cmap(hex_list: List[str], float_list: Optional[List[float]] = None):
    """
    Creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

    :param hex_list: list of hex code strings
    :param float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
    :return: colour map
    """
    if float_list is None:
        float_list = []

    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(["red", "green", "blue"]):
        col_list = [
            [float_list[i], rgb_list[i][num], rgb_list[i][num]]
            for i in range(len(float_list))
        ]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap("my_cmp", segmentdata=cdict, N=256)
    return cmp


def hex_to_rgb(value: str) -> Tuple[int, int, int]:
    """
    Converts hex to rgb colours

    :param value: string of 6 characters representing a hex colour
    :return: list length 3 of RGB values
    """
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i : i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value: int) -> List[int]:
    """
    Converts rgb to decimal colours (i.e. divides each value by 256)

    :param value: list (length 3) of RGB values
    :return: list (length 3) of decimal values
    """
    return [v / 256 for v in value]


def parse_arguments() -> Args:
    parser = argparse.ArgumentParser(
        description="Plot the distribution of quality scores and read lengths in a FASTQ file"
    )

    parser.add_argument("--input", "-i", help="Input FASTQ file", required=True)
    parser.add_argument(
        "--output", "-o", help="Ouput plot file (PDF)", default="quality_plot.pdf"
    )
    parser.add_argument(
        "--min-length", "-l", help="Minimum read length", type=int, default=1550
    )
    parser.add_argument(
        "--max-length", "-L", help="Maximum read length", type=int, default=1750
    )
    parser.add_argument(
        "--min-quality", "-q", help="Minimum quality score", type=int, default=9
    )
    parser.add_argument(
        "--max-quality", "-Q", help="Maximum quality score", type=int, default=28
    )
    parser.add_argument(
        "--figure-width", "-W", help="Figure width", type=int, default=12
    )
    parser.add_argument(
        "--figure-height", "-H", help="Figure height", type=int, default=6
    )
    parser.add_argument(
        "--bins-x", "-x", help="Number of bins in X axis", type=int, default=66
    )
    parser.add_argument(
        "--bins-y", "-y", help="Number of bins in Y axis", type=int, default=20
    )

    return Args(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()
