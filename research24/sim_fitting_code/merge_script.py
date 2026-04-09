from glob import glob
import pandas as pd
import numpy as np

all_files = sorted(glob("./results/chunk_*.pkl"))
merged = pd.concat([pd.read_pickle(f) for f in all_files]).sort_index()
merged.to_pickle("./results/sims_fitted.pkl")
print(f"Done. {len(merged)} rows saved to sims_fitted.pkl")
