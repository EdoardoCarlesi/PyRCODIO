import pandas as pd
import tools as t
import numpy as np

df = pd.read_csv('/home/edoardo/CLUES/DATA/2MRS/groups.csv')

center = np.zeros((3, 1))
x_cols = ['SGX', 'SGY', 'SGZ']

all_x = df[x_cols].T.values
all_d = np.sum((all_x - center) ** 2.0, axis=0)
df['R'] = np.sqrt(all_d)
#df['R'] = df[x_cols].T.apply(lambda x: t.distance(center, x))
print(df.head())
x_vals = df[x_cols].values

triax = t.inertia_tensor(x=x_vals)

print(triax)

r_min = 2
r_max = 12
n_shuffle = 10

t_r = np.zeros((r_max-r_min, n_shuffle, 3))
tw_r = np.zeros((r_max-r_min, n_shuffle, 3))
pca_r = np.zeros((r_max-r_min, n_shuffle, 3))

for r in range(r_min, r_max):

    rad = r * 1.0

    df_tmp = df[df['R'] < rad]
    n_tmp = len(df_tmp)
    
    for j in range(0, n_shuffle):
        n_sample = int(n_tmp/ 2.0)
        x_tmp = df_tmp[x_cols].sample(n_sample, random_state=j).T.values
        ml_tmp = df_tmp['mass_lum'].sample(n_sample, random_state=j).T.values
        t_tmp = np.sort(t.inertia_tensor(x=x_tmp))
        tw_tmp = np.sort(t.inertia_tensor(x=x_tmp, w=ml_tmp))
        pca_tmp = np.sort(t.std_pca(x=x_tmp))

        t_r[r-r_min, j] = t_tmp
        tw_r[r-r_min, j] = tw_tmp
        pca_r[r-r_min, j] = pca_tmp

        print(rad, t_tmp, tw_tmp, pca_tmp)
        print('-')
