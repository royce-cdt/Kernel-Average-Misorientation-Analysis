
# Multi-sample KAM analysis: full + x-axis quarters with labels & point counts



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

kam_list = {}
kam_x_list = {}

#-------------------------#
# Helper: bin + smooth KAM
#-------------------------#
def bin_and_smooth_with_error(y_vals, quantity, bins, min_points, smooth_window):
    binned_mean, binned_std, binned_counts = [], [], []

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
    raw_se = binned_std / np.sqrt(binned_counts)

    smoothed_mean = pd.Series(binned_mean).rolling(
        window=smooth_window, center=True, min_periods=1
    ).mean().to_numpy()

    smoothed_se = pd.Series(raw_se).rolling(
        window=smooth_window, center=True, min_periods=1
    ).mean().to_numpy()

    return binned_mean, smoothed_mean, binned_std, raw_se, smoothed_se

#-------------------------#
# Helper: process dataframe
#-------------------------#
def process_dataframe(df, label, crack_edge_path):
    crack_data = np.load(crack_edge_path, allow_pickle=True).item()
    distance_um_2d = crack_data['distance_um']
    height_px = crack_data['height_px']
    height_um = crack_data['height_um']

    # Pixel size
    pixel_size_y = height_um / height_px
    pixel_size_x = pixel_size_y

    # Map to pixel indices
    df['x_px'] = (df['x'] / pixel_size_x).round().astype(int).clip(0, distance_um_2d.shape[1]-1)
    df['y_px'] = (df['y_zeroed'] / pixel_size_y).round().astype(int).clip(0, distance_um_2d.shape[0]-1)

    # Binning parameters
    bin_width = 3
    min_points = 500
    smooth_window = 2

    # Vertical bins
    y_zeroed_min, y_zeroed_max = df['y_zeroed'].min(), df['y_zeroed'].max()
    bins = np.arange(y_zeroed_min, y_zeroed_max + bin_width, bin_width)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Vertical binning + smoothing
    binned_kam, smooth_kam, binned_std, raw_se, smooth_se = bin_and_smooth_with_error(
        df['y_zeroed'], df['kam'], bins, min_points, smooth_window
    )
    kam_list[label] = (bin_centers, smooth_kam, binned_kam, raw_se)

#-------------------------#
# Main loop
#-------------------------#
fnames = ['zero', 'locA_2', 'locA_1', 'locab', 'locB', 'locC', 'locD']
base_dir = r"C:\Users\oranm\OneDrive - Imperial College London\PHD\HSSCC\KAM Analysis\results"

for fname in fnames:
    print(f"\n--- Processing {fname} ---")

    file_path = os.path.join(base_dir, fname, f"{fname}_kam_gos_export.csv")
    crack_edge_path = os.path.join(base_dir, fname, "distance_from_crack_edge.npy")

    # Load dataframe
    df_full = pd.read_csv(file_path, delimiter=',', skip_blank_lines=True)

    # Combine alpha + beta phases
    if 'kam_beta' in df_full.columns:
        df_full['kam'] = df_full['kam'].combine_first(df_full['kam_beta'])

    # Sort and vertically flip
    df_full = df_full.sort_values(['x', 'y'])
    df_full['kam'] = df_full['kam'].iloc[::-1].reset_index(drop=True)
    df_full['gos'] = df_full['gos'].iloc[::-1].reset_index(drop=True)

    #-------------------------#
    # Compute y_zeroed (distance from crack edge)
    #-------------------------#
    crack_data = np.load(crack_edge_path, allow_pickle=True).item()
    distance_um_2d = crack_data['distance_um']
    pixel_size_y = crack_data['height_um'] / crack_data['height_px']
    pixel_size_x = pixel_size_y
    y_shift = -df_full['y'].min()
    df_full['y_shifted'] = df_full['y'] + y_shift
    df_full['x_px'] = (df_full['x'] / pixel_size_x).round().astype(int).clip(0, distance_um_2d.shape[1]-1)
    df_full['y_px'] = (df_full['y_shifted'] / pixel_size_y).round().astype(int).clip(0, distance_um_2d.shape[0]-1)
    df_full['y_zeroed'] = -distance_um_2d[df_full['y_px'], df_full['x_px']]
    df_full = df_full[df_full['y_zeroed'] > 0]

    # Apply scale only for locC
    if fname == 'locC':
        df_full['y_zeroed'] *= 1.7

    # Process full dataset
    process_dataframe(df_full.copy(), f"{fname}_full", crack_edge_path)

    #-------------------------#
    # Quarter along x-axis
    #-------------------------#
    x_min, x_max = df_full['x'].min(), df_full['x'].max()
    x_edges = np.linspace(x_min, x_max, 5)  # 4 quarters → 5 edges

    quarters = []
    for i in range(4):
        low, high = x_edges[i], x_edges[i + 1]
        label = f"{fname}_x{i+1} ({int(low)}-{int(high)} µm, N={np.sum((df_full['x']>=low)&(df_full['x']<high))})"
        mask = (df_full['x'] >= low) & (df_full['x'] < high)
        quarters.append((label, df_full[mask].copy()))

    # Process each x-quarter
    for label, df_quarter in quarters:
        process_dataframe(df_quarter, label, crack_edge_path)

    #-------------------------#
    # Plot full vs x-quarters
    #-------------------------#
    plt.figure(figsize=(12, 6))
    colors = ["darkorange", "royalblue", "seagreen", "crimson", "purple"]
    labels = [f"{fname}_full"] + [q[0] for q in quarters]

    for i, label in enumerate(labels):
        bin_centers, smooth_kam, binned_kam, raw_se = kam_list[label]
        if label == f"{fname}_full":
            plt.errorbar(
                bin_centers, binned_kam, yerr=raw_se, fmt='o',
                color=colors[i], ecolor=colors[i],
                elinewidth=1, capsize=2, markersize=5, label=f"{label} (raw)"
            )
            plt.plot(bin_centers, smooth_kam, '-', color=colors[i], linewidth=2, label=f"{label} (smoothed)", alpha = 1)
        else:
            plt.errorbar(
                bin_centers, binned_kam, yerr=raw_se, fmt='o',
                color=colors[i], ecolor=colors[i],
                elinewidth=1, capsize=2, markersize=2, label=f"{label} (raw)", alpha = 0.5
            )
            plt.plot(bin_centers, smooth_kam, '--', color=colors[i], linewidth=2, label=f"{label} (smoothed)", alpha = 0.5)

    from itertools import combinations
    
    # Collect smoothed KAM profiles of the x-quarters
    quarter_smooth = [kam_list[q[0]][1] for q in quarters]
    
    # Compute pairwise Pearson correlations
    corrs = []
    for q1, q2 in combinations(quarter_smooth, 2):
        # Align lengths
        min_len = min(len(q1), len(q2))
        q1_trim, q2_trim = q1[:min_len], q2[:min_len]
    
        # Remove NaNs
        mask = ~np.isnan(q1_trim) & ~np.isnan(q2_trim)
        if np.sum(mask) > 2:  # Need at least 3 points to compute correlation
            corrs.append(np.corrcoef(q1_trim[mask], q2_trim[mask])[0,1])
    
    # If no valid correlations, set similarity_score to NaN
    if corrs:
        similarity_score = np.mean(corrs)
    else:
        similarity_score = np.nan

        
        
    plt.axvline(0, color='k', linestyle='--', label='Crack Edge')
    plt.xlabel("Distance from crack edge (µm)", weight='bold')
    plt.ylabel("KAM (°)", weight='bold')
    total_points = len(df_full)
    plt.title(
    f"KAM Profiles: {fname} full vs x-quarters | Total points: {total_points} | Similarity: {similarity_score:.2f}",
    weight='bold'
    )   

    plt.legend(fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.ylim(0,1.5)
    plt.show()
