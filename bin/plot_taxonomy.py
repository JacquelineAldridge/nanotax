#!/usr/bin/env python
"""Provide a command line tool to generate stacked plot from taxonomies."""

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from pathlib import Path 
import seaborn as sns
import glob as glob
import pandas as pd
import numpy as np
import argparse
import logging 
import sys
import re
import os

sns.set_style("white")
plot_pallete = plt.get_cmap("Spectral_r")

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge blast out file with lineages and summary",
        epilog="Example: python TODO ",
    )

    parser.add_argument(
        "-dir",
        "--dir_csv",
        required=True,
        metavar="FILE",
        help="Input file",
    )
        
    return parser.parse_args(argv)


def plot(csv_file, type):
    level = re.search(r'.*/(.*).csv',csv_file).group(1)
    df = pd.read_csv(csv_file) 
    df_melt = df.melt(id_vars=[level]).rename(columns={level:"tax"})
    df_melt = df.melt(id_vars=[level]).rename(columns={level:"tax"})
    
    df_melt = (df_melt               
               .assign(tax=np.where(df_melt['value'] < 1 , "Others", df_melt["tax"]))
               )
    df_melt = df_melt.groupby(["variable","tax"]).sum()
    df_melt = df_melt.reset_index().pivot(index="tax", columns="variable", values="value").reset_index()
    if "Others" in df_melt.tax.unique():
        taxid_order = {
                x: i
                for i, x in enumerate(
                    [x for x in sorted(df_melt.tax.unique()) if x != "Others"] + ["Others"]
                )

        }
    else:
        taxid_order = {
                x: i
                for i, x in enumerate(
                    [x for x in sorted(df_melt.tax.unique()) if x != "Others"] 
                )

        } 
    df_taxs = df_melt.set_index("tax").T
    sorted_series = (df_taxs.sum()).sort_values(ascending=False)
    cmap = ListedColormap(
            sns.color_palette(
                plot_pallete(np.linspace(0, 0.98, num=len(taxid_order) - 1))
            ).as_hex()
            + ["#bbbbbb"]
        )
    plot_kwargs = (
            dict(
                kind="bar",
                stacked=True,
                colormap=cmap,
                figsize=(6,4),
                xlabel=level,
                ylabel="percentage",
                width=1,#0.95,  
                alpha=0.8,
            )
        )
    fig, ax = plt.subplots()
    df_taxs = df_taxs[list(taxid_order.keys())]
    df_taxs.plot(**plot_kwargs, ax=ax)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")
    ax.tick_params(left=True, axis="y", colors="#888888")
    plt.xticks(rotation=30, ha="right")
    ax.legend(
            bbox_to_anchor=(0.5, -0.20),
            ncol=3,
            title="",
            prop={"style": "italic"},
            loc="upper center",
            frameon=False,
        )
    sns.despine(bottom=True, ax=ax)
    fig.savefig(
            f"{type}s/{level}.pdf",
            # f"{outdir}{tax_level}_{suffix_image}.{format_images}",
            format="pdf",
            bbox_inches="tight",
            dpi = 300
    )
    


def main(argv=None):
    args = parse_args(argv)
    os.mkdir(f"{args.dir_csv}s")
    for file in glob.glob(f"{args.dir_csv}/*"):
        plot(file,args.dir_csv)
    
if __name__ == "__main__":
    sys.exit(main())

