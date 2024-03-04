import pandas as pd
import plotly.express as px

table1_file = '/Users/simplexdna/Library/CloudStorage/GoogleDrive-christopher.hempel@kaust.edu.sa/.shortcut-targets-by-id/1c91orCwfstL7NmlFbhdJ8Ssz87pDMHx-/SIREN project/table1.xlsx'

df1 = pd.read_excel(table1_file, header=None, names=["metric", "Species count"])

df1["Category"] = 3 * ["Standardization"] + 5 * ["Databases"]

px.bar(df1, x='metric', y='Species count', color='Category')