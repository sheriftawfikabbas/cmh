import pandas as pd
import glob
import pickle

g = glob.glob('lookup_table_*.pkl')

d = {}

for gg in g:
    f = open(gg,'rb')
    dd = pickle.load(f)
    d.update(dd)

a_file = open("lookup_table.pkl", "wb")

pickle.dump(d, a_file)

a_file.close()


bad = open('bad','w')
for k in d.keys(): 
     if type(k) is float: 
         bad.write(str(k)+'\n') 
bad.close()