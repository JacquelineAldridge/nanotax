#!/usr/bin/env python
"""Provide a command line tool to generate stacked plot from taxonomies."""
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import argparse
import logging
import sys


logger = logging.getLogger()


## J: Add metavar in arguments
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate stacked plot from taxonomies.",
        epilog="Example: python stacked_plot.py family.csv <PARAMS TODO>",
    )

    parser.add_argument(
        "-i",
        "--lefse_input",
        required=True,
        metavar = "FILE",
        help="Input file",
    )

    ## ----------------------------------
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )

    return parser.parse_args(argv)


def plot_lefse(df:pd.DataFrame, conditions: list, name_plot:str):
    # group1 = re.search(r"non_rep.*_(A.*)(C.*)",name_plot).group(1)
    # group2 = re.search(r"non_rep.*_(A.*)(C.*)",name_plot).group(2)
    for var_to_plot in ["feature"]: #['description',"pathway"]:
        ##plt.style.use('seaborn-white')

        fig, ax = plt.subplots(figsize=(6,12))
        # plot a bar chart
        sns.barplot(
            x="value", 
            y=var_to_plot, 
            data=df, 
            estimator=sum, 
            ax=ax,
            palette=["#69b3a2" if x>0 else '#4374B3' for x in df["value"]])



        ax.set_xlabel("LDA score (log 10)", fontsize=12, color="#6b6b6b")
        ax.set_ylabel("", fontsize=12)
        ax.tick_params(axis='x', colors='#7C7E7C') 
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)


        ax.spines['bottom'].set_color('#7C7E7C')
        for i, (value, description) in enumerate(zip(df["value"], df[var_to_plot])):
            x = value + 1 if value > 0 else value - 2  # Ajuste la posición x según el valor
            ha = 'left' if value < 0 else 'right' 
            extra = 0.1 if value < 0 else -0.1 
            ax.annotate(description, (0+extra, i), va='center', fontsize=10, ha=ha, color="#6b6b6b")

        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_xlim(-3, 3)
        patch_1 = mpatches.Patch(color='#69b3a2', label=conditions[0])
        patch_2 = mpatches.Patch(color='#4374B3', label=conditions[1])

        ax.legend(loc='upper right', bbox_to_anchor=(0.1, 1),handles=[patch_1,patch_2])

        plt.savefig(f'{name_plot}_de.pdf',format="pdf",bbox_inches='tight')

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    df = pd.read_csv(f"{args.lefse_input}/{args.lefse_input}_lefse_results.csv", sep="\t", header=None, names=["feature","log","condition","LDA effect size","pvalue"])
    conditions =  df["condition"].unique()


    df = df.query("condition.notnull()")

    df["pvalue"] = df["pvalue"].astype(float)
    if(len(df)>0):

        df = df.assign(value=np.where(df['condition'] == conditions[0], df['LDA effect size'], -df['LDA effect size']),
                        color=np.where(df["condition"] == conditions[1], "blue","red"),
                        pathway = lambda x: x["feature"].str.replace("_","-"),
                        ).sort_values(by="value", ascending=False)
        
        plot_lefse(df, conditions, args.lefse_input)
        df.filter(items=["feature","condition","value"]).to_csv(f"{args.lefse_input}_lefse_df.csv", index=False)

if __name__ == "__main__":
    sys.exit(main())

