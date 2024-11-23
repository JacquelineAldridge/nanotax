#!/usr/bin/env python
"""Provide a command line tool to generate stacked plot from taxonomies."""

import polars as pl
import argparse
import sys
import glob as glob
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge blast out file with lineages and summary",
        epilog="Example: python TODO ",
    )

    parser.add_argument(
        "-csv",
        "--csv",
        required=True,
        help="Input file",
    )
    
    parser.add_argument(
        "-db",
        "--db",
        required=True,
        help="database",
    )


    return parser.parse_args(argv)



def main(argv=None):
    args = parse_args(argv)
    df = pl.LazyFrame()
    categories = ["species","genus","family","order","class","phylum"] if args.db == "genbank" else ["genus","family","order","class","phylum"]
    for category in categories:
        df = pl.concat([(pl.scan_csv(f) 
                       # .assign(sample=re.search(r"(.*)_T1",f) .group(1))
                        )
                        for f in glob.glob(f"*{category}.csv")])
        df.collect().pivot("sample", index=category, values="perc").fill_null(0).write_csv(f"sample/{category}.csv")
        if(category == 'species'):
            df.collect().pivot("sample", index=category, values="count").fill_null(0).write_csv(f"{category}_nreads.csv")

        if("group" in df.collect().columns):
            df_gby = (df
                    .group_by([category,"group"])
                    .agg(pl.col("count").sum().alias("count"))
                    .join(  
                            df.group_by("group")
                            .agg(pl.col("count").sum().alias("total_count")),
                            on="group",
                            how="left"
                        )
                        .with_columns(  
                            perc=(pl.col("count") * 100 / pl.col("total_count")).round(2)
                        )
                    .filter(pl.col("perc") >= 0.01)
                    .with_columns( perc=(pl.col("count") * 100 / pl.col("total_count")).round(2))
                    .collect()
            )   
            df_gby.pivot("group", index=category, values="perc").fill_nan(0).write_csv(f"group/{category}.csv")
            #df.collect().pivot("group", index=category, values="perc").fill_nan(0).write_csv(f"sample/{category}.csv")
            #df_cat.pivot(index=category, columns="sample", values="count").fillna(0).to_csv(f"{category}_by_sample_count.csv",sep=",")


if __name__ == "__main__":
    sys.exit(main())

