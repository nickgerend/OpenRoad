# Written by: Nick Gerend, @dataoutsider
# Viz: "Open Road", enjoy!

import pandas as pd
import os

df = pd.read_csv(os.path.dirname(__file__) + '/USI.csv')

df_out = df.loc[df['PRIME_NAME'].str.contains('Interstate')]

#df_out = df.dropna(subset=['PRIME_NAME'])

# print(df1)
# df2 = pd.read_csv(os.path.dirname(__file__) + '/CA_Counties_Density.csv')

# scale = 100000.
# df2['x'] = df2['x']/scale
# df2['y'] = df2['y']/scale
# df2['z'] = 0.
# df2['side'] = 'base'
# df2['index'] = 0
# df1['item2'] = 0
# df2.rename(columns={'NAME': 'item', 'vertex_part': 'item2', 'vertex_part_index': 'path'}, inplace=True)
# df1['type'] = 'density'
# df2['type'] = 'map'
# df2 = df2[['index', 'item', 'item2', 'side', 'path', 'x', 'y', 'z', 'type']]
# df_out = pd.concat([df1,df2], axis=0)
df_out.to_csv(os.path.dirname(__file__) + '/USI_i_data.csv', encoding='utf-8', index=False)