"""
prepare_input_files.py

Generates the same input files as input_files_creation.ipynb, but sourced
from red_prog_df.pkl, blue_prog_df.pkl, and line_df_local.pkl instead of
the original .npz fitting-parameter files.

Output files (mirroring the original notebook):
  1. wavelengths_<year>_<v>_<J>.sav  — MATLAB .sav with labwave, v_exc, j_exc, v_ground, j_ground
  2. velocity_<year>_<v>_<J>.pkl     — dict {str(labwave): velocity_array}
  3. flux_<year>_<v>_<J>.pkl         — dict {str(labwave): flux_array}         (raw)
  4. fluxerr_<year>_<v>_<J>.pkl      — dict {str(labwave): fluxerr_array}      (raw)
  5. RU_Lupi_line_parameters_<year>_<v>_<J>.csv
"""

import re
import pickle
import numpy as np
import pandas as pd
from scipy.io import savemat
from pathlib import Path

# ── 0. USER PATHS ─────────────────────────────────────────────────────────────

# Input pickle files
PATH_RED_PROG  = r"/mnt/c/Users/prera/Documents/GitHub/scripts_by_me/research24/data/red_df.pkl"
PATH_BLUE_PROG = r"/mnt/c/Users/prera/Documents/GitHub/scripts_by_me/research24/data/blue_df.pkl"
PATH_LINE_DF   = r"/mnt/c/Users/prera/Documents/GitHub/scripts_by_me/research24/data/line_df_coadded_bootstrapped_run2.pkl"

# Output directories  (will be created if missing)
OUT_WAV  = Path(r"/mnt/c/Users/prera/Documents/GitHub/scripts_by_me/research24/data/wavelengths_sav_files")
OUT_SPEC = Path(r"/mnt/c/Users/prera/Documents/GitHub/scripts_by_me/research24/data/spectra_sav_files/spectra")
OUT_LP   = Path(r"/mnt/c/Users/prera/Documents/GitHub/scripts_by_me/research24/data/line_parameters")

for d in [OUT_WAV, OUT_SPEC, OUT_LP]:
    d.mkdir(parents=True, exist_ok=True)

# ── 1. REFERENCE DICTIONARIES (copied verbatim from the notebook) ─────────────

# h2_lines_dict_total = {
#     '[3,16]': [1418.23, 1513.99, 1593.26, 1621.12],
#     '[4,4]':  [1477.05, 1526.55, 1613.72],
#     '[1,7]':  [1467.08, 1500.45, 1524.65, 1556.87],
#     '[1,4]':  [1431.01, 1446.12, 1489.57, 1504.76],
#     '[0,2]':  [1402.65, 1463.83, 1525.15],
#     '[2,12]': [1434.54, 1453.10, 1555.89, 1588.80],
#     '[0,3]':  [1395.20, 1407.29, 1468.39],
#     '[4,13]': [1415.33, 1509.45, 1613.99],
#     '[3,0]':  [1435.05, 1591.32, 1636.34],
#     '[0,1]':  [1398.95, 1460.17, 1521.59],
# }

h2_lines_dict_total = {
    '[3,16]': [1418.23, 1513.99, 1593.26, 1621.12],
    '[4,4]':  [1613.72],
    '[1,7]':  [1467.08, 1500.45],
    '[1,4]':  [1431.01, 1446.12, 1489.57, 1504.76],
    '[0,2]':  [1463.83],
    '[2,12]': [1555.89, 1588.80],
    '[0,3]':  [1407.29, 1468.39],
    '[4,13]': [1415.33],
    '[3,0]':  [1591.32, 1636.34],
    '[0,1]':  [1398.95, 1460.17, 1521.59],
}

dict_line_Bmn = {
    '1415.33': 0.125, '1509.45': 0.023, '1613.99': 0.092,
    '1418.23': 0.050, '1513.99': 0.057, '1593.26': 0.122, '1621.12': 0.062,
    '1477.05': 0.039, '1526.55': 0.033, '1613.72': 0.150,
    '1467.08': 0.080, '1500.45': 0.101, '1524.65': 0.111, '1556.87': 0.074,
    '1431.01': 0.058, '1446.12': 0.083, '1489.57': 0.094, '1504.76': 0.115,
    '1435.05': 0.118, '1591.32': 0.233, '1636.34': 0.099,
    '1398.95': 0.141, '1460.17': 0.083, '1521.59': 0.032,
    '1402.65': 0.126, '1463.83': 0.074, '1525.15': 0.029,
    '1434.54': 0.066, '1453.10': 0.049, '1555.89': 0.077, '1588.80': 0.0119,
    '1395.20': 0.096, '1407.29': 0.120, '1468.39': 0.070,
}

h2_pumping_dict = {
    '[3,13]': 1213.36, '[4,13]': 1213.68, '[3,16]': 1214.47, '[4,4]': 1214.78,
    '[1,7]':  1215.73, '[1,4]':  1216.07, '[3,0]':  1217.04, '[0,1]': 1217.21,
    '[0,2]':  1217.64, '[2,12]': 1217.90, '[2,15]': 1218.52, '[0,3]': 1219.09,
}

line_ID_dict = {
    '1608.44': '(3-9)P(14)',  '1615.43': '(3-10)R(12)',
    '1415.33': '(4-6)R(12)',  '1509.45': '(4-8)R(12)',  '1613.99': '(4-11)R(12)',
    '1418.23': '(3-5)R(15)',  '1513.99': '(3-7)R(15)',  '1593.26': '(3-9)R(15)',  '1621.12': '(3-10)R(15)',
    '1477.05': '(4-8)P(5)',   '1526.55': '(4-9)P(5)',   '1613.72': '(4-11)P(5)',
    '1467.08': '(1-6)P(8)',   '1500.45': '(1-7)R(6)',   '1524.65': '(1-7)P(8)',   '1556.87': '(1-8)R(6)',
    '1431.01': '(1-6)R(3)',   '1446.12': '(1-6)P(5)',   '1489.57': '(1-7)R(3)',   '1504.76': '(1-7)P(5)',
    '1435.05': '(3-7)P(1)',   '1591.32': '(3-10)P(1)',  '1636.34': '(3-11)P(1)',
    '1398.95': '(0-5)P(2)',   '1460.17': '(0-6)P(2)',   '1521.59': '(0-2)P(2)',
    '1402.65': '(0-5)P(2)',   '1463.83': '(0-6)P(3)',   '1525.15': '(0-7)P(3)',
    '1434.54': '(2-5)P(13)',  '1453.10': '(2-6)R(11)',  '1555.89': '(2-8)R(11)',  '1588.80': '(2-8)P(13)',
    '1612.39': '(2-8)P(16)',  '1617.42': '(2-9)R(14)',
    '1395.20': '(0-5)R(2)',   '1407.29': '(0-5)P(4)',   '1468.39': '(0-6)P(4)',
}

# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────

red_prog_df  = pd.read_pickle(PATH_RED_PROG)
blue_prog_df = pd.read_pickle(PATH_BLUE_PROG)
line_df      = pd.read_pickle(PATH_LINE_DF)

# Combine all three frames; drop exact duplicates
all_df = pd.concat([red_prog_df, blue_prog_df, line_df], ignore_index=True)

# ── 3. HELPER: map obsdate → year string ──────────────────────────────────────

def obsdate_to_year(obsdate_str):
    """Return the 4-digit year string from an obsdate like '2011-07-20'."""
    return str(obsdate_str)[:4]

all_df['year'] = all_df['obsdate'].apply(obsdate_to_year)

# ── 4. HELPER: extract Gaussian centre (mu) from params dict ──────────────────

def get_cen_from_params(params, prefer_double=True):
    """
    Extract the best-fit line centre from a params dict.
    params1/params2 look like {'cen': [val, err]} or
    {'cen1': [val, err], 'cen2': [val, err], ...}.
    Returns a float or NaN.
    """
    if params is None or not isinstance(params, dict):
        return np.nan
    if prefer_double:
        # double-Gaussian: take the primary component centre
        if 'cen1' in params:
            v = params['cen1']
            return float(v[0]) if hasattr(v, '__len__') else float(v)
    # single-Gaussian
    if 'cen' in params:
        v = params['cen']
        return float(v[0]) if hasattr(v, '__len__') else float(v)
    # fall back to first 'cenX' key found
    for k in params:
        if k.startswith('cen'):
            v = params[k]
            return float(v[0]) if hasattr(v, '__len__') else float(v)
    return np.nan

# ── 5. CHOOSE BEST PARAMS PER ROW ────────────────────────────────────────────
# Use BIC to decide single vs. double Gaussian (lower BIC = better fit).
# Fall back gracefully if BIC values are missing.

def best_params(row):
    bic1 = row.get('BIC1', np.nan)
    bic2 = row.get('BIC2', np.nan)
    if pd.isna(bic1) and pd.isna(bic2):
        return row.get('params1', None)
    if pd.isna(bic2) or (not pd.isna(bic1) and bic1 <= bic2):
        return row.get('params1', None)
    return row.get('params2', None)

all_df['best_params'] = all_df.apply(best_params, axis=1)
all_df['mu'] = all_df['best_params'].apply(get_cen_from_params)

# ── 6. NORMALISE labline to float for matching ────────────────────────────────

all_df['labline_float'] = all_df['labline'].apply(
    lambda x: float(x.value) if hasattr(x, 'value') else float(x)
)

# ── 7. BUILD A LOOKUP: (year, progression, labwave_float) → row ──────────────

def transition_to_progression(trans):
    """
    Convert a transition label like '[1,7]' or '[2, 12]' to the canonical
    key used in h2_lines_dict_total, e.g. '[1,7]'.
    """
    if isinstance(trans, str):
        trans = trans.strip()
        # normalise spacing: '[1, 7]' -> '[1,7]'
        trans = re.sub(r'\s+', '', trans)
        return trans
    return None

all_df['progression'] = all_df['transition'].apply(transition_to_progression)

# ── 8. FILE CREATORS ──────────────────────────────────────────────────────────

COLUMNS_WAV = ['labwave', 'v_exc', 'j_exc', 'v_ground', 'j_ground']
CSV_COLS    = ["Line_ID", "['v', J']", 'Wavelength', 'B_mn',
               'Pumping Wavelength', 'Observed Wavelength', 'Velocity Shift']

def parse_line_id(line_str):
    """Parse '(v_exc-v_ground)B(j_ground)' into (v_exc, v_ground, branch, j_ground)."""
    m = re.match(r'\((\d+)-(\d+)\)([PR])\((\d+)\)', line_str)
    if not m:
        return None
    v_exc    = int(m.group(1))
    v_ground = int(m.group(2))
    branch   = m.group(3)
    j_ground = int(m.group(4))
    j_exc    = j_ground + 1 if branch == 'R' else j_ground - 1
    return v_exc, j_exc, v_ground, j_ground


def create_wavelengths_file(year, progression):
    """
    Create wavelengths_<year>_<v>_<J>.sav with labwave, v_exc, j_exc,
    v_ground, j_ground for every line in h2_lines_dict_total[progression].
    This file is purely derived from reference data (no observation needed).
    """
    if progression not in h2_lines_dict_total:
        return
    rows = []
    for lw in h2_lines_dict_total[progression]:
        key = f"{lw:.2f}"
        if key not in line_ID_dict:
            continue
        parsed = parse_line_id(line_ID_dict[key])
        if parsed is None:
            continue
        v_exc, j_exc, v_ground, j_ground = parsed
        rows.append({
            COLUMNS_WAV[0]: float(f"{lw:.2f}"),
            COLUMNS_WAV[1]: v_exc,
            COLUMNS_WAV[2]: j_exc,
            COLUMNS_WAV[3]: v_ground,
            COLUMNS_WAV[4]: j_ground,
        })
    if not rows:
        return
    df_out = pd.DataFrame(rows)
    v_tag  = progression[1]   # e.g. '1' from '[1,7]'
    j_tag  = progression[3]   # e.g. '7'
    fname  = OUT_WAV / f"wavelengths_{year}_{v_tag}_{j_tag}.sav"
    savemat(str(fname), df_out.to_dict('list'), appendmat=False)
    print(f"  Saved {fname.name}")


def create_spectra_files(year, progression, df_year_prog):
    """
    Create velocity_<year>_<v>_<J>.pkl, flux_<year>_<v>_<J>.pkl,
    fluxerr_<year>_<v>_<J>.pkl.

    df_year_prog: sub-DataFrame already filtered to this year+progression.
    For each labline in h2_lines_dict_total[progression], we grab the
    wavelength, flux, and fluxerr arrays from the first matching row
    (rows share the same spectral arrays per line+epoch combination).
    """
    if progression not in h2_lines_dict_total:
        return

    vel_dict      = {}
    flux_dict     = {}
    fluxerr_dict  = {}

    for lw in h2_lines_dict_total[progression]:
        lw_f = float(lw)
        # Match rows whose labline is this wavelength
        mask = np.isclose(df_year_prog['labline_float'], lw_f, atol=0.01)
        sub  = df_year_prog[mask]
        if sub.empty:
            continue

        row = sub.iloc[0]

        # wavelength array → velocity array
        wl_arr = row['wavelength']
        if hasattr(wl_arr, 'value'):          # SpectralAxis → numpy
            wl_arr = wl_arr.value
        wl_arr = np.asarray(wl_arr, dtype=float)
        vel_arr = ((wl_arr - lw_f) / lw_f) * 3e5

        flux_arr    = np.asarray(row['flux'],    dtype=float)
        fluxerr_arr = np.asarray(row['fluxerr'], dtype=float)

        key = str(lw)
        vel_dict[key]     = vel_arr
        flux_dict[key]    = flux_arr
        fluxerr_dict[key] = fluxerr_arr

    if not vel_dict:
        return

    v_tag = progression[1]
    j_tag = progression[3]

    for label, data in [('velocity', vel_dict),
                         ('flux',     flux_dict),
                         ('fluxerr',  fluxerr_dict)]:
        fpath = OUT_SPEC / f"{label}_{year}_{v_tag}_{j_tag}.pkl"
        with open(fpath, 'wb') as f:
            pickle.dump(data, f)
        print(f"  Saved {fpath.name}")


def create_line_params_file(year, progression, df_year_prog):
    """
    Create RU_Lupi_line_parameters_<year>_<v>_<J>.csv.
    Uses the best-fit Gaussian centre (mu) to compute observed wavelength
    and velocity shift, mirroring the notebook's get_info_from_data().
    """
    if progression not in h2_lines_dict_total:
        return
    if progression not in h2_pumping_dict:
        return

    pumping_wl = h2_pumping_dict[progression]
    rows = []

    for lw in h2_lines_dict_total[progression]:
        lw_f   = float(lw)
        lw_str = f"{lw_f:.2f}"
        if lw_str not in line_ID_dict or lw_str not in dict_line_Bmn:
            continue

        mask = np.isclose(df_year_prog['labline_float'], lw_f, atol=0.01)
        sub  = df_year_prog[mask]
        if sub.empty:
            continue

        row = sub.iloc[0]
        mu  = float(row['mu']) if not pd.isna(row['mu']) else np.nan
        if np.isnan(mu):
            continue

        vel_shift = ((mu - lw_f) / lw_f) * 3e5

        rows.append({
            CSV_COLS[0]: line_ID_dict[lw_str],
            CSV_COLS[1]: progression,
            CSV_COLS[2]: lw_str,
            CSV_COLS[3]: dict_line_Bmn[lw_str],
            CSV_COLS[4]: pumping_wl,
            CSV_COLS[5]: mu,
            CSV_COLS[6]: vel_shift,
        })

    if not rows:
        return

    df_csv = pd.DataFrame(rows, columns=CSV_COLS)
    v_tag  = progression[1]
    j_tag  = progression[3]
    fpath  = OUT_LP / f"RU_Lupi_line_parameters_{year}_{v_tag}_{j_tag}.csv"
    df_csv.to_csv(fpath, index=False)
    print(f"  Saved {fpath.name}")


# ── 9. ALSO CREATE A MODELLED-FLUX VERSION OF THE SPECTRA FILES ──────────────
# Your dataframes contain both raw and empirically modelled flux/fluxerr.
# This function writes the Gaussian-model versions with a '_model' suffix.

def create_model_spectra_files(year, progression, df_year_prog):
    """
    Like create_spectra_files but uses gauss_flux / gauss_fluxerr (single
    Gaussian) or double_gauss_flux / double_gauss_fluxerr (double Gaussian),
    selected by BIC, as the flux arrays.  Velocity array is identical to raw.
    """
    if progression not in h2_lines_dict_total:
        return

    flux_dict    = {}
    fluxerr_dict = {}

    for lw in h2_lines_dict_total[progression]:
        lw_f = float(lw)
        mask = np.isclose(df_year_prog['labline_float'], lw_f, atol=0.01)
        sub  = df_year_prog[mask]
        if sub.empty:
            continue

        row = sub.iloc[0]

        bic1 = row.get('BIC1', np.nan)
        bic2 = row.get('BIC2', np.nan)

        if pd.isna(bic2) or (not pd.isna(bic1) and bic1 <= bic2):
            # single-Gaussian model
            fl  = row.get('gauss_flux',    row['flux'])
            fle = row.get('gauss_fluxerr', row['fluxerr'])
        else:
            # double-Gaussian model
            fl  = row.get('double_gauss_flux',    row['flux'])
            fle = row.get('double_gauss_fluxerr', row['fluxerr'])

        key = str(lw)
        flux_dict[key]    = np.asarray(fl,  dtype=float)
        fluxerr_dict[key] = np.asarray(fle, dtype=float)

    if not flux_dict:
        return

    v_tag = progression[1]
    j_tag = progression[3]
    for label, data in [('flux_model',    flux_dict),
                         ('fluxerr_model', fluxerr_dict)]:
        fpath = OUT_SPEC / f"{label}_{year}_{v_tag}_{j_tag}.pkl"
        with open(fpath, 'wb') as f:
            pickle.dump(data, f)
        print(f"  Saved {fpath.name}")


# ── 10. MAIN LOOP ─────────────────────────────────────────────────────────────

progressions_list = list(h2_lines_dict_total.keys())
years_present     = sorted(all_df['year'].unique())

print(f"Years found in data: {years_present}")
print(f"Progressions to process: {progressions_list}\n")

for year in years_present:
    print(f"=== Year {year} ===")
    df_year = all_df[all_df['year'] == year]

    for prog in progressions_list:
        df_yp = df_year[df_year['progression'] == prog]

        # Wavelengths .sav — independent of whether data exists
        create_wavelengths_file(year, prog)

        if df_yp.empty:
            continue

        print(f"  [{prog}] {len(df_yp)} rows found")
        create_spectra_files(year, prog, df_yp)
        create_model_spectra_files(year, prog, df_yp)
        create_line_params_file(year, prog, df_yp)

print("\nDone.")