import sys
import pickle as pkl
import pandas as pd
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
from empirical_functions import ModelFitting
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

line_mask_dict = {
    1453.1:   200,
    1555.89:  [75, 175],
    1463.83:  200,
    1613.72:  [75, 200],
    1415.33:  180,
    1407.29:  [75, 250],
    1468.39:  [75, 299],
    1435.05:  [80, 299],
    1636.34:  [75, 190],
}
save_dir = '/blue/ast7939/pkottapalli/results'
def fit_worker(args):
    chunk_df, disp_df, lsf_df, worker_id = args
    fitter = ModelFitting(chunk_df, line_mask_dict, disp_df, lsf_df, verbose=False)
    fitter.fit_lines()

    # Each worker writes its own file immediately on completion
    out_path = save_dir+f"/chunk_{worker_id:04d}.pkl"
    fitter.df.to_pickle(out_path)
    print(f"[Worker {worker_id}] done, saved to {out_path}", flush=True)

    return fitter.df

def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else cpu_count()
    chunk_size = int(sys.argv[2]) if len(sys.argv) > 2 else 5  # rows per chunk

    print("Loading data...", flush=True)
    sims_df = pd.read_pickle('data/sims_df.pkl')
    disp_df = pd.read_pickle('data/BC_conf_dispdf.pkl')
    lsf_df  = pd.read_pickle('data/BC_conf_LSFdf.pkl')

    # Check for already-completed chunks and skip them
    import os, glob
    completed_files = glob.glob(save_dir+f"/chunk_*.pkl")
    completed_indices = set()
    for f in completed_files:
        try:
            completed_indices.update(pd.read_pickle(f).index)
        except Exception:
            pass  # skip corrupted files

    if completed_indices:
        print(f"Skipping {len(completed_indices)} already-completed rows.", flush=True)
        sims_df = sims_df[~sims_df.index.isin(completed_indices)]

    n_rows = len(sims_df)
    if n_rows == 0:
        print("All rows already completed!")
        return

    # Split into small chunks — pool will hand these out dynamically
    indices = range(0, n_rows, chunk_size)
    chunks = [sims_df.iloc[i:i+chunk_size].copy() for i in indices]
    n_chunks = len(chunks)
    print(f"Launching {n_workers} workers over {n_rows} rows in {n_chunks} chunks of {chunk_size}...", flush=True)

    os.makedirs("save_dir", exist_ok=True)
    worker_args = [(chunk, disp_df, lsf_df, i) for i, chunk in enumerate(chunks)]


    results = []
    ctx = mp.get_context("forkserver")  # forkserver is safer than fork with nested pools
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
        futures = [executor.submit(fit_worker, args) for args in worker_args]
    for i, future in enumerate(as_completed(futures)):
        results.append(future.result())
        print(f"[{i+1}/{len(worker_args)} chunks complete]", flush=True)

    # Final merge of all result files (includes previous runs)
    all_files = sorted(glob.glob(save_dir+"/chunk_*.pkl"))
    merged = pd.concat([pd.read_pickle(f) for f in all_files]).sort_index()
    merged.to_pickle(save_dir+"/sims_fitted.pkl")
    print(f"Done. {len(merged)} rows saved to sims_fitted.pkl")

if __name__ == "__main__":
    main()
