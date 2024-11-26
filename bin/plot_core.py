#!/usr/bin/env python
"""Provide a command line tool to generate stacked plot from taxonomies."""

import numpy as np
import pandas as pd
import argparse
import logging 
import plotly.graph_objects as go
import plotly.express as px 
import sys
import json
import glob as glob

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge blast out file with lineages and summary",
        epilog="Example: python TODO ",
    )

    parser.add_argument(
        "-tsv",
        "--mmseqs_tsv",
        required=True,
        metavar="FILE",
        help="Input file",
    )
    parser.add_argument(
        "-gr",
        "--groups",
        required=False,
        help="groups dict",
    ) 
    parser.add_argument(
        "-tax",
        "--taxlineage_csv",
        required=True,
        metavar="FILE",
        help="Input file",
    )  
    return parser.parse_args(argv)


def generate_core_plot(df_input:str,columns_to_core:list=[] ,group:str = " " ,colors:list = ["#FFFFFF","#88c3b5" ,"#7da0cf","#e4b7a8","#cfbfcb"],categories:list = ["phylum","class_","order","family","genus","species"]):
    df_data= pd.DataFrame(columns = ["taxonomie","parent","value"])
    sum_all = 0
    if(columns_to_core == "all"):
        columns_to_core = list(df_input["sample"].unique())
    df_input = df_input.rename(columns = {"class":"class_"})
    for i,level in enumerate(categories):
        tax_up1 = list(set((df_input.filter(items=[level,"value","sample"]).query("sample in @columns_to_core").query("value >= 1"))[level]))
        df_subset = (df_input.query(f"{level} in @tax_up1")
                 .filter(items = [level, "sample","value"])
                 .groupby([level, "sample"])
                 .agg("sum")
                 .pivot_table(index=level, columns='sample', values='value')
                 .fillna(0)
                 .reset_index()
                 )

        c = [level] + columns_to_core
        df_subset = df_subset.filter(items = [level] + columns_to_core)
        shared = set(df_subset[~(df_subset.iloc[:, 1:] < 1).any(axis=1)][level])

        df_subset = df_subset.query(f"{level} in @shared")
        df_subset= df_subset.groupby(level).sum().reset_index()
        df_subset["value"] = df_subset.set_index(level).mean(axis=1).reset_index()[0] 
        df_subset = df_subset.filter(items=[level,"value"])
        if(level == "species"):
            df_subset = df_subset.query("~species.str.contains('mixed')")
        if(level == "phylum"):
            df_subset = df_subset.rename(columns= {level:"taxonomie"})
            df_subset["parent"] = f"Core"

        else:
            df_p =  df_subset.merge(df_input.filter(items=[level,parent]),how="left").filter(items=[level,parent,"value"]).drop_duplicates(keep="first")
            df_subset = df_p.rename(columns= {level:"taxonomie", parent:"parent"})

        df_data = pd.concat([df_data, df_subset])
        parent = level
    taxonomies = list(df_data["taxonomie"])
    taxonomies_parent = list(df_data["parent"])
    taxonomies_percentage = list(df_data["value"])
    sum_all = df_data.query("parent=='Core'")["value"].sum()
    df_data.to_csv(f"core_plot_{group}{level}.csv")
    taxonomies_percentage = [sum_all] + taxonomies_percentage
    df = pd.DataFrame({'node_names':["Core"] + taxonomies,
                    'node_parent':  [""] + taxonomies_parent,
                    'node_counts': taxonomies_percentage
                    })
    n_colors =  len(df_data.query("parent == 'Core'"))
    colors = colors[0:n_colors+1]
    fig=go.Figure(
        data=go.Sunburst(
            ids=df["node_names"],
            labels=df["node_names"], 
            parents=df["node_parent"],
            values=df["node_counts"],
            branchvalues="total",
            marker=dict(colors=colors,
                )
        ),
        
    )
    fig.update_layout(
        autosize=False,
        width=800,
        height=800)
    fig.update_traces(textfont=dict(
            size=14,color="white"
        ))
    
    fig.write_image(f"core_{group}.pdf")


def main(argv=None):
    args = parse_args(argv)
    df_lineage = pd.concat([pd.read_csv(file)
                            for file in glob.glob("*taxlineage.csv")])
    df_species = pd.read_csv("last_assignment.csv")
    if("species" in df_species.columns):
        df_melt = pd.melt(df_species, id_vars=['species'],var_name='sample')
        df_lineage=df_lineage.rename(columns={"taxname":"species"})
        cat = ["phylum","class_","order","family","genus","species"]
    else:
        df_melt = pd.melt(df_species, id_vars=['genus'],var_name='sample')
        cat = ["phylum","class_","order","family","genus"]


    df_merge = pd.merge(df_melt,df_lineage,how="left").drop_duplicates(keep="first")

    ## Core all samples
    generate_core_plot(df_merge,"all","all",categories=cat)

    ## Core by groups
    groups = args.groups[1:len(args.groups)-1].split(", ")
    groups_dict = ({item.split(":")[0]: item.split(":")[1] for item in groups})
    groups = set((groups_dict.values()))
    groups.discard('false')
    for group in groups:
        samples_group = {key: value for key, value in groups_dict.items() if value == group}
        samples_group = samples_group.keys()
        generate_core_plot(df_merge,list(samples_group),group,categories=cat)

    

if __name__ == "__main__":
    sys.exit(main())

