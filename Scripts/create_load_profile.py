#I create the load file
from create_network import n
import numpy as np
import pandas as pd

n_load=35 #Number of loads

p_set=np.random.rand(len(n.snapshots),(n_load))

array=p_set*100

df=pd.DataFrame(array)

filepath = 'my_file.xlsx'

df.to_excel(filepath, index=False)


