# Written by: Nick Gerend, @dataoutsider
# Viz: "Open Road", enjoy!

import pandas as pd
import os

data = pd.read_csv(os.path.dirname(__file__) + '/Metro_Interstate_Traffic_Volume.csv.gz',compression='gzip', error_bad_lines=False)
print(data)
data.to_csv(os.path.dirname(__file__) + '/i_volume.csv', encoding='utf-8-sig', index=False)