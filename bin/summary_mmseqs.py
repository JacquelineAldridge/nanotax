#!/usr/bin/env python
"""Provide a command line tool to generate stacked plot from taxonomies."""

import polars as pl
import argparse
import sys

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
        "-db",
        "--db",
        required=True,
        help="database",
    )

    parser.add_argument(
        "-id",
        "--min_identity",
        required=False,
        default=0.95,
        metavar="FLOAT",
        type=float,
        help="Input file",
    )

    parser.add_argument(
        "-aln",
        "--min_aln",
        required=False,
        default=1000,
        metavar="FLOAT",
        type=float,
        help="Input file",
    )
    parser.add_argument(
        "-sample",
        "--sample",
        required=True,
        help="sample",
    )
    parser.add_argument(
        "-group",
        "--group",
        required=True,
        help="group",
    )
        
    return parser.parse_args(argv)



def main(argv=None):
    args = parse_args(argv)
    cutoff = 0.6
    ## ToDo para  Silva
    if(args.db == 'silva'):
        print("silva")
        select_colums = ["query","genus","family","order","class","phylum"]
    if(args.db == 'genbank'):
        print("genbank",args)
        select_colums = ["query","taxname","genus","family","order","class","phylum"]

    
    names = "query,target,pident,tcov,alnlen,taxname,taxlineage".split(",")
    df = (pl.scan_csv(args.mmseqs_tsv,separator="\t",has_header=False, new_columns=names)
          .filter((pl.col("pident")>=(args.min_identity)*100) & (pl.col("alnlen")>=args.min_aln))
          .with_columns(         
            pl.col("taxlineage").str.extract(r".*;g_([^;]*);._.*", 1).alias("genus"),   
            pl.col("taxlineage").str.extract(r".*;f_([^;]*);._.*", 1).alias("family"),
            pl.col("taxlineage").str.extract(r".*;o_([^;]*);._.*", 1).alias("order"),
            pl.col("taxlineage").str.extract(r".*;c_([^;]*);._.*", 1).alias("class"),
            pl.col("taxlineage").str.extract(r".*;p_([^;]*);._.*", 1).alias("phylum"),
              pident = pl.col("pident").round(1),
              tcov = pl.col("tcov").round(3))
          .sort(["pident","tcov"], descending=True)
          )
    ## If a read have more than one aln with the same specie (genbank) or genus (silva), keep only the first aln 
    df = df.unique(subset=["query","taxname"], keep="first")

    
    ## Case 1: uniq assignment
    df_uniq = (df
        .unique(subset="query", keep="none")
         .with_columns(
            pl.col("taxlineage").str.extract(r".*;g_([^;]*);._.*", 1).alias("genus"),
            pl.col("taxlineage").str.extract(r".*;f_([^;]*);._.*", 1).alias("family"),
            pl.col("taxlineage").str.extract(r".*;o_([^;]*);._.*", 1).alias("order"),
            pl.col("taxlineage").str.extract(r".*;c_([^;]*);._.*", 1).alias("class"),
            pl.col("taxlineage").str.extract(r".*;p_([^;]*);._.*", 1).alias("phylum")
         )
       .select(["query","taxname","genus","family","order","class","phylum"])  
    )
    uniq_ids = set(df_uniq.select(["query"]).collect().to_series())
    
    ## Reads with secondary alns
    df_first = (df
                .filter(~pl.col("query").is_in(uniq_ids))
                .sort(["pident","tcov"], descending=True)
                .unique(subset=["query"], keep="first")
                .rename(mapping={"target":"target_first","pident":"pident_first","taxname":"taxname_first","genus":"genus_first","family":"family_first","order":"order_first","class":"class_first","phylum":"phylum_first"})
    )
    df_merge =  (df_first
                 .join(df, on=["query"], how = "left")
                 .with_columns(pident_diff = (pl.col("pident_first") - pl.col("pident")).round(1))
                 .filter(pl.col("pident_diff") < cutoff)
    )


    ## Caso 2: Secondary alns have a difference with fist aln higher than 0.6 so we select first aln
    uniq_ids_after_cutoff = (df_merge .group_by(["query"])
                        .agg(pl.count()
                        .alias("count"))
                        .filter(pl.col("count") == 1)
                        .select(["query"])
                        .collect()#["query"]
                        # .unique()
                        # .to_list()
    )
    df_non_valid_sec_alns = (df_first
                             .filter(pl.col("query").is_in(uniq_ids_after_cutoff))
                             .rename(mapping={"target_first":"target","pident_first":"pident","taxname_first":"taxname","genus_first":"genus","family_first":"family","order_first":"order","class_first":"class","phylum_first":"phylum"})
                             .select(["query","taxname","genus","family","order","class","phylum"])
                             )
    
    ## Case 3 and 4: Reads with valid secondary alns
    df_valid_sec_alns = df_merge.filter(~pl.col("query").is_in(uniq_ids_after_cutoff))
                         
    ## Case 3: Reads with secondary alns: Alns with different genus
    ids_dif_genus = set(df_valid_sec_alns.filter(pl.col("genus") !=pl.col("genus_first") ).select(["query"]).collect().to_series())
    df_dif_genus = df_valid_sec_alns.filter(pl.col("query").is_in(ids_dif_genus))


    if(len(df_dif_genus.collect()) == 0):
        df_dif_genus_final = pl.LazyFrame({ "query": [], "taxname": [], "genus": [], "family": [], "order": [], "class": [], "phylum": []})
    else:
        df_dif_genus_final = ( df_dif_genus.group_by("query").agg([
            pl.concat_str([
                pl.col("genus").unique().sort(),
                pl.col("genus_first").unique().sort()
            ], separator=" ").alias("genus"),
            pl.concat_str([
                pl.col("family").unique().sort(),
                pl.col("family_first").unique().sort()
            ], separator=" ").alias("family"),

            pl.concat_str([
                pl.col("order").unique().sort(),
                pl.col("order_first").unique().sort()
            ], separator=" ").alias("order"),

            pl.concat_str([
                pl.col("class").unique().sort(),
                pl.col("class_first").unique().sort()
            ], separator=" ").alias("class"),
            pl.concat_str([
                pl.col("phylum").unique().sort(),
                pl.col("phylum_first").unique().sort()
            ], separator=" ").alias("phylum")

            ])
            .with_columns(pl.col("genus").list.join(" ").alias("genus") )
            .with_columns(pl.col("genus").str.split(" ")  
                    .list.unique()              
                    .list.join(" / ")         
                    .alias("genus")             
                    )
            .with_columns((pl.lit("mixed: ") + pl.col("genus")).alias("genus"))

            .with_columns((pl.lit("mixed genus: ") + pl.col("genus")).alias("taxname"))
            ##
            .with_columns(pl.col("family").list.join(" ").alias("family") )
            .with_columns(pl.col("family").str.split(" ")  
                    .list.unique()              
                    .list.join(" / ")         
                    .alias("family")             
                    )
            ##
            .with_columns(pl.col("order").list.join(" ").alias("order") )
            .with_columns(pl.col("order").str.split(" ")  
                    .list.unique()              
                    .list.join(" / ")         
                    .alias("order")             
                    )
            ##
            .with_columns(pl.col("class").list.join(" ").alias("class") )
            .with_columns(pl.col("class").str.split(" ")  
                    .list.unique()              
                    .list.join(" / ")         
                    .alias("class")             
                    )
            ##
            .with_columns(pl.col("phylum").list.join(" ").alias("phylum") )
            .with_columns(pl.col("phylum").str.split(" ")  
                    .list.unique()              
                    .list.join(" / ")         
                    .alias("phylum")             
                    )

        )
        ## Remove those reads that have alns with more than two genus
        df_dif_genus_final = df_dif_genus_final.filter(pl.col("genus").str.count_matches("/").lt(2))
        df_dif_genus_final = df_dif_genus_final.select(["query","taxname","genus","family","order","class","phylum"])  

        

    ## Case 4: Reads with secondary alns: Alns with same genus (always)
    id_same_genus = set(df_valid_sec_alns
                     .filter(~pl.col("query").is_in(ids_dif_genus))
                     .select(["query"])
                     .collect()
                     .to_series()
                     )
    df_same_genus = (df_first
                      .filter(pl.col("query").is_in(id_same_genus))
                      .rename(mapping={"target_first":"target","pident_first":"pident","taxname_first":"taxname","genus_first":"genus","family_first":"family","order_first":"order","class_first":"class","phylum_first":"phylum"})
                      .select(["query","taxname","genus","family","order","class","phylum"])
                    )
    ## Concat everything

    df_final = pl.concat([df_uniq,df_non_valid_sec_alns,df_dif_genus_final,df_same_genus])
    df_final = df_final.rename(mapping = {"taxname":"species"})
    for category in ["species","genus","family","order","class","phylum"]:
        df_tax = (df_final
                    .group_by([category])
                    .agg(pl.count().alias("count"))
                    .sort(["count"], descending=True)
                    .with_columns(perc = (pl.col("count")/pl.col("count").sum()*100).round(2))
                    .collect()
                    .filter(pl.col("perc") >= 0.01)
                    .with_columns(perc = (pl.col("count")/pl.col("count").sum()*100).round(2))
        )

        if(args.group!='false'):
            df_tax = df_tax.with_columns( pl.lit(args.group).alias("group"),
                                          pl.lit(args.sample).alias("sample"))
        else:
            df_tax = df_tax.with_columns(pl.lit(args.sample).alias("sample"))
        df_tax.write_csv(f"{args.sample}_{category}.csv")

if __name__ == "__main__":
    sys.exit(main())

