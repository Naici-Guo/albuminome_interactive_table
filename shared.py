from pathlib import Path
import numpy as np
import pandas as pd

app_dir = Path(__file__).parent
#df = pd.read_csv(app_dir / "penguins.csv")

data_index = pd.read_csv(app_dir / "data_index.csv", encoding="macroman")
data_index = data_index.replace(r'^\s*$',np.nan, regex=True)
data_main= pd.read_csv(app_dir/ "data_main.csv")
data_main = data_main.replace(r'^\s*$',np.nan, regex=True)