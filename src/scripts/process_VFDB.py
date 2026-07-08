import pandas as pd

columns = ["Gene", "Description", "Organism", "Category", "VFID"]

df = pd.read_csv("/Data/VFDB/VFs.csv", sep="\t", header=0)

