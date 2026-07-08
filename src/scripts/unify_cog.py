import pandas as pd

wp_cog = pd.read_csv("/Data/COG/cog-24.cog.csv", header=None, names=["1", "GCF", "GCA", "2", "3", "4", "COG no.", "5", "6", "7", "8", "9", "10"])

cog_cat = pd.read_csv("/Data/COG/cog-24.mapping.tab", sep="\t")

# Create e mapping from the GCA in wp_cog to the Func columns in cog_cat going through COG no. column in both dataframes
wp_cog_cat = pd.merge(wp_cog, cog_cat, left_on="COG no.", right_on="COG no.", how="left")

wp_cog_cat = wp_cog_cat[["GCA", "COG no.", "Func"]]

import pdb; pdb.set_trace()