import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

fnames = ['zero', 'locA_2', 'locA_1', 'locab', 'locB', 'locC', 'locD']
kam_list = {}
kam_x_list = {}  # 🔷 Store horizontal averages

for fname in fnames:
    print(fname)
    # Load CSV data 
    file_path = os.path.join(
        r"C:\Users\oranm\OneDrive - Imperial College London\PHD\HSSCC\KAM Analysis\results",
        fname, f"{fname}_kam_gos_export.csv")
    df = pd.read_csv(file_path, delimiter=',', skip_blank_lines=True)
    
    # combine alpha nd beta phases
    df['kam'] = df['kam'].combine_first(df['kam_beta'])
    
    # Sort and vertically flip KAM/GOS
    df = df.sort_values(['x', 'y'])
    df['kam'] = df['kam'].iloc[::-1].reset_index(drop=True)
    df['gos'] = df['gos'].iloc[::-1].reset_index(drop=True)

    if fname == 'zero':
        # 🔹 Compute the mean of the zero binned KAM
        zero_mean_kam = np.nanmean(df['kam'])
        print(f"Zeropoint mean KAM: {zero_mean_kam:.4f}")


    # Load crack edge data
    crack_edge_path = os.path.join(
        r"C:\Users\oranm\OneDrive - Imperial College London\PHD\HSSCC\KAM Analysis\results",
        fname, "distance_from_crack_edge.npy")
    crack_data = np.load(crack_edge_path, allow_pickle=True).item()
    distance_um_2d = crack_data['distance_um']
    height_px = crack_data['height_px']
    height_um = crack_data['height_um']

    # Pixel size (assuming square pixels)
    pixel_size_y = height_um / height_px
    pixel_size_x = pixel_size_y

    # Shift y to positive domain for indexing
    y_shift = -df['y'].min()
    df['y_shifted'] = df['y'] + y_shift
    
    df = df.dropna(subset=['x', 'y', 'kam'])


    # Map to pixel indices
    
    df['x_px'] = (df['x'] / pixel_size_x).round().astype(int)
    df['y_px'] = (df['y_shifted'] / pixel_size_y).round().astype(int)
    df['y_px'] = df['y_px'].clip(0, distance_um_2d.shape[0] - 1)
    df['x_px'] = df['x_px'].clip(0, distance_um_2d.shape[1] - 1)

    # Assign crack-relative vertical coordinate
    df['y_zeroed'] = -distance_um_2d[df['y_px'], df['x_px']]
    df = df[df['y_zeroed'] > 0]
    
    ''' # --- Plot a vertical profile before binning ---
    if fname == fnames[6]:
        # Pick one x-column (e.g. the middle one)
        unique_x = sorted(df['x'].unique())
        x_val = unique_x[-1]   # choose middle column
        sub = df[df['x'] == x_val]

        # Plot just this one profile
        plt.figure(figsize=(10, 6))
        plt.plot(
            sub['y_zeroed'], sub['kam'],
            color="steelblue",
            marker="*",
            linestyle="None",
            label=f"x={x_val:.2f}"
        )

        plt.xlabel("Distance from crack edge (µm)")
        plt.ylabel("KAM (°)")
        plt.title(f"Single KAM profile vs distance from crack edge ({fname})")
        plt.axvline(0, color="k", linestyle="--", label="Crack Edge")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


    '''



    # Binning parameters
    bin_width = 2 #microns
    min_points = 200
    smooth_window = 2

    y_zeroed_min, y_zeroed_max = df['y_zeroed'].min(), df['y_zeroed'].max()
    bins = np.arange(y_zeroed_min, y_zeroed_max + bin_width, bin_width)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # 🔷 Horizontal binning (averaging along y)
    x_min, x_max = df['x'].min(), df['x'].max()
    x_bins = np.arange(x_min, x_max + bin_width, bin_width)
    x_bin_centers = (x_bins[:-1] + x_bins[1:]) / 2

    

    # Binning and smoothing function
    def bin_and_smooth_with_error(y_vals, quantity, bins, min_points, smooth_window):
        binned_mean = []
        binned_std = []
        binned_counts = []
        for i in range(len(bins) - 1):
            in_bin = (y_vals >= bins[i]) & (y_vals < bins[i + 1])
            values = quantity[in_bin]
            count = np.sum(in_bin)
            binned_counts.append(count)
            if count >= min_points:
                binned_mean.append(np.nanmean(values))
                binned_std.append(np.nanstd(values))
            else:
                binned_mean.append(np.nan)
                binned_std.append(np.nan)
        
        binned_mean = np.array(binned_mean)
        binned_std = np.array(binned_std)
        binned_counts = np.array(binned_counts)
        
        # Calculate raw standard error
        raw_se = binned_std / np.sqrt(binned_counts)
        
        # Smooth mean and smooth se
        smoothed_mean = pd.Series(binned_mean).rolling(window=smooth_window, center=True, min_periods=1).mean().to_numpy()
        smoothed_se = pd.Series(raw_se).rolling(window=smooth_window, center=True, min_periods=1).mean().to_numpy()
        
        return binned_mean, smoothed_mean, binned_std, raw_se, smoothed_se
    

    # Apply binning and smoothing
    binned_kam, smooth_kam, binned_std, raw_se, smooth_se = bin_and_smooth_with_error(
        df['y_zeroed'], df['kam'], bins, min_points, smooth_window)
    kam_list[fname] = (bin_centers, smooth_kam, binned_kam, raw_se)

    binned_kam_x, smooth_kam_x, _, raw_se_x, _ = bin_and_smooth_with_error(
        df['x'], df['kam'], x_bins, min_points, smooth_window)
    kam_x_list[fname] = (x_bin_centers, smooth_kam_x, binned_kam_x, raw_se_x)
    
    # Normalise binned and smoothed KAM curves
    # Normalise by initial value (preserve decay shape and gradient)
    # 🔷 Normalise by zeropoint mean KAM (preserve shape, valid gradients)
    #binned_kam_norm = binned_kam / zero_mean_kam
    #smooth_kam_norm = smooth_kam / zero_mean_kam
    #kam_list[fname + '_norm'] = (bin_centers, smooth_kam_norm, binned_kam_norm, raw_se / zero_mean_kam)


    '''
    # Plot GOS vs Distance from Crack Edge
    plt.figure()
    plt.plot(bin_centers, binned_gos, 'g--', label='Raw GOS')
    plt.plot(bin_centers, smooth_gos, 'm-', label='Smoothed GOS')
    plt.xlabel('Distance from crack edge (µm)')
    plt.ylabel('GOS (°)')
    plt.title('GOS vs Distance from Crack Edge')
    plt.axvline(0, color='k', linestyle='--', label='Crack Edge')
    plt.legend()
    plt.grid(True)

    plt.show()
    '''
    #print(f"Number of points per bin {fname}:", [np.sum((df['y_zeroed'] >= bins[i]) & (df['y_zeroed'] < bins[i+1])) for i in range(len(bins)-1)])

    '''
    # Scatter plot: KAM in original coordinates
    plt.figure(figsize=(12,6))
    sc1 = plt.scatter(df['x'], df['y'], c=df['kam'], cmap='viridis', s=0.1)
    plt.colorbar(sc1, label='KAM (°)')
    plt.xlabel('x (µm)')
    plt.ylabel('y (µm)')
    plt.title(f'KAM in Original Coordinate Space {fname}')
    plt.gca().invert_yaxis()  # Flip y to match image coordinates if needed
    plt.gca().invert_xaxis()
    plt.grid(True)

    # Scatter plot: KAM in crack-referenced coordinates
    plt.figure(figsize=(12, 6))
    sc2 = plt.scatter(df['x'], df['y_zeroed'], c=df['kam'], cmap='viridis', s=0.1)
    plt.colorbar(sc2, label='KAM (°)')
    plt.xlabel('x (µm)')
    plt.gca().invert_xaxis()
    plt.ylabel('Distance from crack edge (µm)')
    plt.title(f'KAM in Crack-Referenced Coordinate Space {fname}')
    plt.axhline(0, color='k', linestyle='--', label='Crack Edge')
    plt.legend()
    plt.grid(True)


    plt.show()
    '''

#%% Comparison Plot
color_map = {
    'locA_1': 'silver',
    'locA_2': 'grey',
    'locab': 'green',
    'locB': 'slateblue',
    'locC': 'darkorange',
    'locD': 'gold',
    'zero': 'mediumseagreen'
}

annotations1 = {
    'locA_2': ("MAP A_2, location = 2.3mm (brittle region)", (33, 0.825)),
    'locA_1': ("MAP A_1, location = 2.5mm (brittle region)", (30, 0.95)),   # (x, y) coords to place text
    'locB': ("MAP B, location = 5.5mm (brittle/plastic transition)", (5, 0.6)),
    'locC': ("MAP C, location = 8.5mm (plastic region)", (38, 0.70)),
    'locD': ("MAP D, location = 9.5mm (plastic region)", (24, 0.75)),
    'zero': ("ZEROPOINT", (10, 0.45))
}

annotations2 = {}
'''
plt.figure(figsize=(15, 10))

for fname in fnames:
    bin_centers, smooth_kam, binned_kam, raw_se = kam_list[fname]
    color = color_map.get(fname, 'gray')
    
    # Compute range and save it
    kam_range = np.nanmax(binned_kam) - np.nanmin(binned_kam)
    #valid_start = next(val for val in binned_kam if not np.isnan(val))
    #valid_end = next(val for val in reversed(binned_kam) if not np.isnan(val))
    #kam_range = valid_start - valid_end


    annotations2[fname] = kam_range
    
    # Plot raw binned mean with raw SE error bars
    plt.errorbar(
        bin_centers, binned_kam, yerr=raw_se, fmt='o', color=color,
        ecolor=color, elinewidth=1, capsize=2, markersize=4
    )
    
    # Add annotation
    if fname in annotations1:
        text, (x_pos, y_pos) = annotations1[fname]
        plt.text(x_pos, y_pos, text, color=color, fontsize=12, weight='bold')

    
    # Plot smoothed mean without error bars (skip for 'zero')
    if fname != 'zero':
        plt.plot(bin_centers, smooth_kam, '--', color=color, linewidth=2, label=f'{fname} smoothed')

    
# Add horizontal line at mean KAM value for the zeropoint
if 'zero' in kam_list:
    _, _, binned_kam_zero, _ = kam_list['zero']
    mean_zero_kam = np.nanmean(binned_kam_zero)
    plt.axhline(mean_zero_kam, color='mediumseagreen', linestyle='--', linewidth=1.5, label=f"Zeropoint Mean KAM = {mean_zero_kam:.2f}°")


plt.text(72, 0.75, f"ΔKAM A_1 = {round(annotations2['locA_1'],2)}°", color=color_map['locA_1'], fontsize=15)
plt.text(72, 0.7, f"ΔKAM A_2 = {round(annotations2['locA_2'],2)}°", color=color_map['locA_2'], fontsize=15)
plt.text(72, 0.65, f"ΔKAM B = {round(annotations2['locB'],2)}°", color=color_map['locB'], fontsize=15)
plt.text(72, 0.60, f"ΔKAM C = {round(annotations2['locC'],2)}°", color=color_map['locC'], fontsize=15)
plt.text(72, 0.55, f"ΔKAM D = {round(annotations2['locD'],2)}°", color=color_map['locD'], fontsize=15)

plt.xlabel('Distance from crack edge (µm)', weight='bold')
plt.ylabel('KAM (°)', weight='bold')
plt.title('KAM (Averaged Along Crack Edge) vs Distance From Crack Edge ', weight='bold')
plt.axvline(0, color='k', linestyle='--', label='Crack Edge')
plt.legend()
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')
plt.grid(False)
plt.show()
'''
#%% Comparison Plot (Normalised)
'''
plt.figure(figsize=(15, 10))

for fname in fnames:
    norm_key = fname + '_norm'
    if norm_key in kam_list:
        bin_centers, smooth_kam_norm, binned_kam_norm, raw_se = kam_list[norm_key]
        color = color_map.get(fname, 'gray')

        # Plot normalised binned mean
        plt.plot(bin_centers, smooth_kam_norm, '-', color=color, linewidth=2, label=f'{fname} (norm)')

        # Optional: plot raw binned with markers
        plt.plot(bin_centers, binned_kam_norm, 'o', color=color, alpha=0.3)

plt.axvline(0, color='k', linestyle='--', label='Crack Edge')
plt.xlabel('Distance from crack edge (µm)', weight='bold')
plt.ylabel('Normalised KAM (0–1)', weight='bold')
plt.title('Normalised KAM Decay Curves', weight='bold')
plt.legend()
plt.grid(True)
plt.show()
'''

#%% Model Fitting
import matplotlib.pyplot as plt


'''

plt.figure(figsize=(12, 6))
for fname in fnames:
    bin_centers, smooth_gos, binned_gos = gos_list[fname]
    color = color_map.get(fname, 'gray')  # fallback color

    # Plot smoothed and binned with same color, different linestyle
    plt.plot(bin_centers, smooth_gos, '-', color=color, linewidth=2, label=f'{fname} (smoothed)')
    plt.plot(bin_centers, binned_gos, '--', color=color, linewidth=1, alpha=0.7, label=f'{fname} (binned)')

plt.xlabel('Distance from crack edge (µm)')
plt.ylabel('GOS (°)')
plt.title('GOS vs Distance from Crack Edge')
plt.axvline(0, color='k', linestyle='--', label='Crack Edge')
plt.legend()
plt.grid(True)
plt.show()
'''

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

def clean_xy(x, y):
    x = np.array(x)
    y = np.array(y)
    mask = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isinf(x)) & (~np.isinf(y))
    return x[mask], y[mask]

def inv_sqrt_model(x, a, y0, c_fixed):
    return y0 + a**2 / np.sqrt(x * c_fixed)

def residuals_inv_sqrt(params, x, y, c_fixed):
    a, y0 = params
    return inv_sqrt_model(x, a, y0, c_fixed) - y

group1 = ['locA_1', 'locA_2', 'locab']


group2 = ['locB', 'locC', 'locD']

# Set bounds to enforce a > 0
lower_bounds = [0, 0]  # a > 0, y0 unbounded
upper_bounds = [np.inf, np.inf]

# Assuming kam_list and color_map are defined elsewhere
# kam_list = {'locA_1': (x_values, y_values), ...}
# color_map = {'locA_1': 'r', 'locB': 'b', ...}

all_bin_centers = []
for group in [group1, group2]:
    for fname in group:
        bc = np.array(kam_list[fname][0])
        all_bin_centers.append(bc)
max_bin_center = max(np.max(bc) for bc in all_bin_centers)
general_bin_centers = np.linspace(0, max_bin_center, 300)

fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)


results = []
c_fixed = 2*np.pi  # to avoid division by zero

for ax, group, title in zip(axes, [group1, group2], ['Corrosion-Dominated Region (i)', 'Fatigue-Dominated Region (ii)']):
    for fname in group:
        bin_centers, smooth_kam, binned_kam, raw_se = kam_list[fname]
        color = color_map.get(fname, 'gray')

        # Plot smooth data with alpha=0.3
        ax.plot(bin_centers, smooth_kam, '-', color=color, linewidth=1, label=f'{fname} smoothed data', alpha=0.5)

        # Clean and restrict to x >= 0
        mask_nonneg = bin_centers >= 0
        xdata = bin_centers[mask_nonneg]
        ydata = binned_kam[mask_nonneg]
        xdata_clean, ydata_clean = clean_xy(xdata, ydata)
        
        if ydata_clean.size == 0:
            print(f"Warning: No valid data points for {fname}, skipping fit.")
            continue  # Skip this file


        # Initial parameter guess
        a_init = max(np.max(ydata_clean) - np.min(ydata_clean), 0.1)
        y0_init = ydata_clean[-1]
        p0_inv = [a_init, y0_init]

        # Inverse sqrt fit
        res_inv = least_squares(
            residuals_inv_sqrt, p0_inv,
            args=(xdata_clean, ydata_clean, c_fixed),
            bounds=(lower_bounds, upper_bounds),
            loss='soft_l1'
        )

        if res_inv.success:
            a_opt, y0_opt = res_inv.x
            y_fit = inv_sqrt_model(general_bin_centers, a_opt, y0_opt, c_fixed)
            
            # NEW: evaluate fit at 80 µm
            kam_at_80 = inv_sqrt_model(80, a_opt, y0_opt, c_fixed)
    
            ss_res = np.sum((ydata_clean - inv_sqrt_model(xdata_clean, a_opt, y0_opt, c_fixed)) ** 2)
            ss_tot = np.sum((ydata_clean - np.mean(ydata_clean)) ** 2)
            r2 = 1 - (ss_res / ss_tot)
            rmse = np.sqrt(np.mean((ydata_clean - inv_sqrt_model(xdata_clean, a_opt, y0_opt, c_fixed)) ** 2))

            try:
                J = res_inv.jac
                cov = np.linalg.pinv(J.T @ J) * np.sum(res_inv.fun**2) / (len(ydata_clean) - len(res_inv.x))
                stderr = np.sqrt(np.diag(cov))
            except np.linalg.LinAlgError:
                stderr = [np.nan, np.nan]

            a_stderr, y0_stderr = stderr

            ax.plot(general_bin_centers, y_fit, '--', color=color, alpha=0.9,
                    label=f'{fname} inv sqrt fit (R²={r2:.2f})')
        else:
            a_opt, y0_opt, r2, rmse, a_stderr = [np.nan] * 5

        dy_dx = np.gradient(ydata_clean, xdata_clean)
        smoothness = np.std(dy_dx)
        decay_range = ydata_clean[0] - ydata_clean[-1]

        results.append({
            'fname': fname,
            'a_fit_inv': a_opt,
            'a_stderr': a_stderr,
            'y0_fit_inv': y0_opt,
            'y0_stderr': y0_stderr,
            'R2_inv': r2,
            'RMSE_inv': rmse,
            'Smoothness': smoothness,
            'DecayRange': decay_range,
            'KAM_80um': kam_at_80
        })

        # Raw binned data with error bars at alpha=0.3
        ax.errorbar(
            bin_centers, binned_kam, yerr=raw_se, fmt='o', color=color,
            ecolor=color, elinewidth=1, capsize=2, markersize=4, alpha=0.3
        )

    # Add zeropoint to both plots for context
    if 'zero' in kam_list:
        z_bin_centers, z_smooth_kam, _, _ = kam_list['zero']
        ax.plot(z_bin_centers, z_smooth_kam, 'o', alpha=0.3, color='mediumseagreen', label='zeropoint data', markersize=4)
        ax.axhline(zero_mean_kam, color='mediumseagreen', label='zeropoint mean', linewidth=2)

    ax.axvline(0, color='k', linestyle='--', label='crack edge')
    ax.set_title(title, weight='bold', fontsize=15)
    ax.set_xlabel('Distance from crack edge (µm)', weight='bold')
    ax.grid(True)
    ax.set_ylabel('KAM (°)', weight = 'bold')
    ax.set_ylim(0.4, 1.5)
    if ax == axes[0]:
        ax.legend(loc='right')   # For left half
    else:
        ax.legend(loc='upper right')  # For right half

    

fig.suptitle('Inverse Square Model Fitting to KAM vs Distance from Crack Edge', fontsize=18, weight='bold')

plt.tight_layout()
plt.show()

print("\nInverse Square Root Fitting Summary:")
print(f"{'Name':<8} {'K_fit':>10} {'y0':>10} {'R2':>10} {'RMSE':>10} "
      f"{'Smooth':>10} {'Decay':>10} {'KAM@80µm':>12}")
for r in results:
    print(f"{r['fname']:<8} "
          f"{r['a_fit_inv']:10.3f} {r['y0_fit_inv']:10.3f} "
          f"{r['R2_inv']:10.3f} {r['RMSE_inv']:10.3f} "
          f"{r['Smoothness']:10.3f} {r['DecayRange']:10.3f} "
          f"{r['KAM_80um']:12.3f}")
    
print(f'zeropoint mean = {round(zero_mean_kam, 3)}')

#%%
import numpy as np
import matplotlib.pyplot as plt

# --- Crack distance mapping ---
crack_distance_map = {
    'locA_2': 2.3,
    'locA_1': 2.5,
    'locab': 4,
    'locB': 5.5,
    'locC': 8.5,
    'locD': 9.5
}

# --- Extract data from results ---
crack_dist, a_fits, rmse, r2_vals, fnames, a_stderr_arr = [], [], [], [], [], []

for res in results:
    fname = res['fname']
    if fname in crack_distance_map:
        crack_dist.append(crack_distance_map[fname])
        a_fits.append(res['a_fit_inv'])
        rmse.append(res['RMSE_inv'])
        r2_vals.append(res['R2_inv'])
        fnames.append(fname)
        a_stderr_arr.append(res.get('a_stderr', np.nan))

# Convert to arrays
crack_dist = np.array(crack_dist)
a_fits = np.array(a_fits)
rmse = np.array(rmse)
r2_vals = np.array(r2_vals)
fnames = np.array(fnames)
a_stderr_arr = np.array(a_stderr_arr)

# Sort by crack distance
sort_idx = np.argsort(crack_dist)
crack_dist = crack_dist[sort_idx]
a_fits = a_fits[sort_idx]
rmse = rmse[sort_idx]
r2_vals = r2_vals[sort_idx]
fnames = fnames[sort_idx]
a_stderr_arr = a_stderr_arr[sort_idx]

# --- Split fatigue vs hydrogen ---
fatigue_mask = ~np.isin(fnames, ['locA_1', 'locA_2'])
sqrt_crack_dist = np.sqrt(crack_dist)

x_fatigue = sqrt_crack_dist[fatigue_mask]
y_fatigue = a_fits[fatigue_mask]
stderr_fatigue = a_stderr_arr[fatigue_mask]

x_hydrogen = sqrt_crack_dist[~fatigue_mask]
y_hydrogen = a_fits[~fatigue_mask]
stderr_hydrogen = a_stderr_arr[~fatigue_mask]

# --- Remove NaNs and zero errors for weights ---
valid_mask = (~np.isnan(y_fatigue)) & (~np.isnan(stderr_fatigue)) & (stderr_fatigue > 0)
x_fit = x_fatigue[valid_mask]
y_fit = y_fatigue[valid_mask]
weights = 1 / (stderr_fatigue[valid_mask] ** 2)
weights /= np.max(weights)  # normalize

# --- Weighted slope-only fit through origin ---
slope_kfit = np.sum(weights * x_fit * y_fit) / np.sum(weights * x_fit**2)
y_pred = slope_kfit * x_fit

# --- Weighted R² for fit through origin ---
mean_y_weighted = np.average(y_fit, weights=weights)
ss_res = np.sum(weights * (y_fit - y_pred) ** 2)
ss_tot = np.sum(weights * (y_fit - mean_y_weighted) ** 2)
r_squared = 1 - ss_res / ss_tot

# --- Prepare smooth curve for plotting ---
x_plot = np.linspace(0, np.max(x_fit), 100)
y_plot = slope_kfit * x_plot

# --- Plot ---
plt.figure(figsize=(8,5))

plt.errorbar(x_fit, y_fit, yerr=stderr_fatigue[valid_mask], fmt='k*', capsize=4, label='Fatigue')
plt.errorbar(x_hydrogen, y_hydrogen, yerr=stderr_hydrogen, fmt='bo', capsize=4, label='Corrosion')

# Updated legend with slope in physical units
plt.plot(x_plot, y_plot, 'r--', 
         label=f'Weighted fit through origin:\nSlope = {slope_kfit:.3f} MPa·√m per √mm\nR² = {r_squared:.4f}')

# Annotate all points
for x, y, r2, c in zip(sqrt_crack_dist, a_fits, r2_vals, fnames):
    plt.annotate(f"R²={r2:.2f}", (x, y-0.2), textcoords="offset points", xytext=(5, 5), ha='left', fontsize=8)
    plt.annotate(f"{c}", (x, y-0.3), textcoords="offset points", xytext=(5, 5), ha='left', fontsize=8)

plt.xlabel('√(Crack length) (mm½)')
plt.ylabel('K_fit from inverse square decay fit (MPa·√m)')
plt.title('Weighted fit of K_fit vs √(Crack length)')
plt.grid(True)
plt.ylim(0,)
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()


print(f'Weighted slope-only fit through origin: slope = {slope_kfit:.3f}, R² = {r_squared:.4f}')


# burgers
'''
# Burgers vector for titanium (in microns)
b = 0.000295  # 0.295 nm -> microns
step_size = 0.085  # µm

plt.figure(figsize=(8,6))

for fname, (bin_centers, smooth_kam, _, _) in kam_list.items():
    # Convert KAM from degrees to radians
    kam_rad = np.deg2rad(smooth_kam)
    
    # Calculate Nye GND density (1/µm^2)
    rho_gnd = kam_rad / (b * step_size)
    
    plt.plot(bin_centers, rho_gnd, label=fname, linewidth=2)

plt.xlabel('Distance from crack edge (µm)')
plt.ylabel('GND Density (1/µm²)')
plt.title('Nye GND Density vs Distance from Crack Edge')
plt.axvline(0, color='k', linestyle='--', label='Crack Edge')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
'''

# GEOM FACTOR IRWIN

import numpy as np
import matplotlib.pyplot as plt

# --- Crack distance mapping ---
crack_distance_map = {
    'locA_2': 2.3,
    'locA_1': 2.5,
    'locab': 4,
    'locB': 5.5,
    'locC': 8.5,
    'locD': 9.5
}

def irwin_rp(a_micron, sigma_y, sigma, w_micron, stress_unit='MPa'):
    a = a_micron * 1e-3 # mm to m
    w = w_micron * 1e-3 # mm to m
    c = a_micron / w_micron

    MG = 0.1*c**2 + 0.29*c + 1.081
    MB = 0.75*c**2 - 0.185*c + 1.019
    MS = 0.9*c**2 - 0.21*c + 1.02
    Y = MG * MB * MS * (2/np.pi)

    if stress_unit == 'MPa':
        sigma_pa = sigma * 1e6
        sigma_y_pa = sigma_y * 1e6
    else:
        sigma_pa = sigma
        sigma_y_pa = sigma_y

    K = Y * sigma_pa * np.sqrt(np.pi * a)  # Pa·√m
    rp_plane_stress = (1/np.pi) * (K / sigma_y_pa)**2
    rp_plane_strain = (1/(6*np.pi)) * (K / sigma_y_pa)**2

    return Y, K / 1e6, rp_plane_stress * 1e6, rp_plane_strain * 1e6

# Parameters for legend
sigma = 700
sigma_y = 1000
w_micron = 1000

# ------------------ PLOT 1: K vs sqrt(a) ------------------
sqrt_a = np.array([np.sqrt(a) for a in crack_distance_map.values()])
K_vals = np.array([irwin_rp(a, sigma_y, sigma, w_micron)[1] for a in crack_distance_map.values()])

# Compute slope through origin for K vs sqrt(a)
slope_K = np.sum(sqrt_a * K_vals) / np.sum(sqrt_a**2)

plt.figure()
plt.plot(sqrt_a, K_vals, 'b^')
plt.plot([0, np.max(sqrt_a)], [0, slope_K * np.max(sqrt_a)], 'r--', label=f'Slope = {slope_K:.3f} MPa·√m per √mm')
plt.title('K vs sqrt a')
plt.xlabel('√(Crack length) [√mm]')
plt.ylabel('K [MPa·√m]')
plt.legend(title=f'σ={sigma} MPa, σ_y={sigma_y} MPa, w={w_micron} μm')
plt.grid()
plt.show()

# ------------------ PLOT 2: Plastic zone vs a ------------------
a_vals = np.array(list(crack_distance_map.values()))
rp_ps_vals = np.array([irwin_rp(a, sigma_y, sigma, w_micron)[2] for a in a_vals])
rp_pe_vals = np.array([irwin_rp(a, sigma_y, sigma, w_micron)[3] for a in a_vals])

plt.figure()
plt.plot(a_vals, rp_ps_vals, 'k^', label='Plane stress')
plt.plot(a_vals, rp_pe_vals, 'r*', label='Plane strain')
plt.title('PZ vs a')
plt.xlabel('Crack length a [mm]')
plt.ylabel('Plastic zone rp [μm]')
plt.legend(title=f'σ={sigma} MPa, σ_y={sigma_y} MPa, w={w_micron} μm')
plt.grid()
plt.show()

print('Factor difference between K_fit slope and literature: m_K_fit / m_K = ' + str(slope_kfit/slope_K))
