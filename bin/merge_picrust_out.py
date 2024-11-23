#!/usr/bin/env python
"""Provide a command line tool to generate stacked plot from taxonomies."""
from pathlib import Path
import pandas as pd
import argparse
import logging
import sys
import re
import glob
from collections import Counter
logger = logging.getLogger()


## J: Add metavar in arguments
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge blast out file with lineages and summary",
        epilog="Example: python TODO ",
    )

    parser.add_argument(
        "--groups",
        metavar = "",
        help="Groups",
        default=14,
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

def create_lefse_input(df, groups,groups_to_compare, func):
    df_picrustT = (df.T.reset_index())
    df_groups = (pd.DataFrame(groups.items(), columns=['index','group']))
    df_groups["index"] = df_groups["index"].apply(lambda x: x.strip())

    df_groups = df_groups.query("group in @groups_to_compare")
    df_groups["group"] = df_groups["group"].apply(lambda x: x.replace(" ","_"))

    df = (df_groups.merge(df_picrustT,on='index', how="left")).rename(columns={"index":"sample"})

    df = df.set_index("group").T.reset_index()
    df.to_csv(f"lefseinput_{func}.tsv",index=False, sep="\t")    

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    
    groups = args.groups[1:len(args.groups)-1].split(", ")
    groups_dict = ({item.split(":")[0]: item.split(":")[1] for item in groups})
    groups = set((groups_dict.values()))
    groups.discard('false')
    group_counts = Counter(groups_dict.values())
    groups_to_compare = [group for group, count in group_counts.items() if count >= 3]
    
    if(len(groups_to_compare) == 1):
        print("Nada que comparar, solo hay un grupo")
        
    for f in ["EC","KO"]:
        dfs = []
        df_picrust_out = pd.concat(
            (pd.read_csv(file, sep="\t", compression='gzip', index_col=0) 
            for file in glob.glob(f"*/{f}_metagenome_out/pred_metagenome_unstrat.tsv.gz")),
            axis=1, join='outer'
        ).fillna(0).round(2)

        df_picrust_out.to_csv(f"{f}.csv",index=True)
        if(len(groups_to_compare) == 2):
            create_lefse_input(df_picrust_out, groups_dict,groups_to_compare, f)

    ## Pathways
    dfs = []
    df_picrust_out = pd.concat(
            (pd.read_csv(file, sep="\t", compression='gzip', index_col=0) 
            for file in glob.glob(f"*/pathways_out/path_abun_unstrat.tsv.gz")),
            axis=1, join='outer'
    ).fillna(0).round(2)
    df_picrust_out.to_csv(f"pathways.csv",index=True)
    if(len(groups_to_compare) == 2):
        create_lefse_input(df_picrust_out, groups_dict,groups_to_compare, "pathways")

    # ToDo: More than two groups

if __name__ == "__main__":
    sys.exit(main())