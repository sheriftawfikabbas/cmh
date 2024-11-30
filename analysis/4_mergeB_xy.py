import pandas as pd
import glob
import numpy as np

subfolder = "1"

g = glob.glob(subfolder+'/Bxy_*.csv')
df = pd.read_csv(g[0],header=0,index_col=0)
l = df.values
for i in g[1:]:
    df = pd.read_csv(i,header=0,index_col=0)
    l = np.concatenate((l,df.values))
df = pd.DataFrame(l)
df.sort_values(by=[0,1], inplace=True)
df.to_csv('Bxy.csv')
