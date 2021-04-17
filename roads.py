# Written by: Nick Gerend, @dataoutsider
# Viz: "Open Road", enjoy!

import pandas as pd
pd.set_option('display.max_rows', None)
import os
import math

def haversine_miles(lat1, lon1, lat2, lon2):
    degrees_to_radians = math.pi/180.0
    dLat = (lat2 - lat1) * degrees_to_radians
    dLon = (lon2 - lon1) * degrees_to_radians
    lat1 = (lat1) * degrees_to_radians
    lat2 = (lat2) * degrees_to_radians
    h = (pow(math.sin(dLat / 2), 2) + pow(math.sin(dLon / 2), 2) * math.cos(lat1) * math.cos(lat2))
    earth_radius = 3958.8
    return  2 * earth_radius * math.asin(math.sqrt(h)) 

def heading(lat1, lon1, lat2, lon2):
    degrees_to_radians = math.pi/180.0
    radians_to_degrees = 180.0/math.pi
    lat1 = (lat1) * degrees_to_radians
    lat2 = (lat2) * degrees_to_radians
    dLon = (lon2 - lon1) * degrees_to_radians
    X = math.cos(lat2) * math.sin(dLon)
    Y = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dLon)
    bearing_degrees = math.atan2(X,Y) * radians_to_degrees
    return (bearing_degrees + 360) % 360

df = pd.read_csv(os.path.dirname(__file__) + '/Interstates.csv')

#region data prep
df['CF_Miles'] = df.apply(lambda x: haversine_miles(x['S_W_Lat'], x['S_W_Long'], x['N_E_Lat'], x['N_E_Long']), axis=1)
df['CF_Miles'] = [450. if x == 'I-69' else y for x, y in zip(df['Interstate'], df['CF_Miles'])]
df['CF_Miles'] = [450. if x == 'I-49' else y for x, y in zip(df['Interstate'], df['CF_Miles'])]
df['CF_Miles'] = [800. if x == 'I-74' else y for x, y in zip(df['Interstate'], df['CF_Miles'])]
df['Diff_Miles'] = df['Length_Miles'] - df['CF_Miles']
df['Diff_Miles'] = [0.5 if x < 0. else x for x in df['Diff_Miles']]
df['Heading'] = df.apply(lambda x: heading(x['S_W_Lat'], x['S_W_Long'], x['N_E_Lat'], x['N_E_Long']), axis=1)
df['arc_ratio'] = (df['Length_Miles']-(df['Length_Miles']-df['Diff_Miles']))/(df['Length_Miles']-df['Diff_Miles'])
#endregion

df.to_csv(os.path.dirname(__file__) + '/Interstates_Clean.csv', encoding='utf-8-sig', index=False)

ratio_scale = 5.

