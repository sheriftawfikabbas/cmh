import pandas as pd

b2a = 1
positions_df = pd.read_csv(
    'results.lammpstrj', sep='\s+', lineterminator='\n', header=None, skiprows=0, index_col=0)
print(len(positions_df))
positions_labelled_df = pd.read_csv(
    'initial.lammpstrj', sep='\s+', header=None, skiprows=0, index_col=0)
positions_labelled_df = positions_labelled_df[4]
positions_df = pd.merge(positions_df, positions_labelled_df,
                        left_index=True, right_index=True)


positions_df = positions_df.loc[(positions_df['4_y'] == 1.5) | (positions_df['4_y'] == 0)]
print(positions_df.head())


positions_df.loc[positions_df.iloc[:, 0] == 1, 'Atom'] = 'Li'
positions_df.loc[positions_df.iloc[:, 0] == 2, 'Atom'] = 'e'
export_df = pd.DataFrame()
export_df['Atom'] = positions_df.Atom
positions_df.iloc[:, 2] *= b2a
positions_df.iloc[:, 3] *= b2a
positions_df.iloc[:, 4] *= b2a
positions_df.iloc[:, 5] *= b2a
positions_df.iloc[:, 6] *= b2a
positions_df.iloc[:, 7] *= b2a
export_df['x'] = positions_df.iloc[:, 2]
export_df['y'] = positions_df.iloc[:, 3]
export_df['z'] = positions_df.iloc[:, 4]
export_df['vx'] = positions_df.iloc[:, 5]
export_df['vy'] = positions_df.iloc[:, 6]
export_df['vz'] = positions_df.iloc[:, 7]
export_df['conduction'] = positions_df.iloc[:, 8]
electrons_df = export_df.loc[export_df.Atom == 'e']

print('e drift =', electrons_df.loc[(electrons_df.x < 0)].vz.mean())
print('e drift =', electrons_df.loc[(electrons_df.x > 0)].vz.mean())
print('all electrons: average velocity along x:', electrons_df.vx.mean())
print('all electrons: average velocity along y:', electrons_df.vy.mean())
print('all electrons: average velocity along z:', electrons_df.vz.mean())

electrons_cond_df = export_df.loc[(export_df.Atom == 'e') & export_df.conduction]
print('conduction electrons: average velocity along z:', electrons_cond_df.vz.mean())

electrons_surface_df = export_df.loc[(export_df.Atom == 'e') & (
    (abs(export_df.x) > 20) & (abs(export_df.y) > 20) & (abs(export_df.z) > 20))]
print('e surface drift =', electrons_surface_df.vx.mean(),
      electrons_surface_df.vy.mean(), electrons_surface_df.vz.mean())


positions_df.reset_index(inplace=True)
positions_df.columns = ['original_index','atom_type','charge','x','y','z','vx','vy','vz','take','atom']
positions_original_index = positions_df.original_index
del positions_df['original_index']
positions_df.to_csv('positions.csv',index=True,header=True)
positions_original_index.to_csv('positions_original_index.csv',index=True,header=True)