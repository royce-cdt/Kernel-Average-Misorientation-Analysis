# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 15:11:31 2025

@author: oranm
"""

import numpy as np
from PIL import Image
import os
import matplotlib.pyplot as plt

step_size = 0.085  # microns per pixel


def load_image_rgb(image_path: str) -> np.ndarray:
    """Load an image and verify it is RGB."""
    img = Image.open(image_path)
    plt.imshow(img)
    img_np = np.array(img)
    if img_np.ndim != 3 or img_np.shape[2] < 3:
        raise ValueError("Image must be RGB with at least 3 channels")
    return img_np


def find_pink_mask(img_np: np.ndarray, target_rgb=(193, 0, 81), tolerance=30) -> np.ndarray:
    """Return a boolean mask for pixels matching pink color within tolerance."""
    r, g, b = img_np[:, :, 0], img_np[:, :, 1], img_np[:, :, 2]
    mask = (
        (np.abs(r - target_rgb[0]) < tolerance) &
        (np.abs(g - target_rgb[1]) < tolerance) &
        (np.abs(b - target_rgb[2]) < tolerance)
    )
    return mask


def calculate_crack_edge_per_column(image_path: str) -> dict:
    """
    Calculates the crack edge location as the maximum y pixel with pink color
    for each column (x) in the image.

    Returns:
        dict with:
            - crack_edge_y_px: 1D numpy array (width) of crack edge pixel y position
            - crack_edge_y_um: same array converted to microns
            - image info
    """
    img_np = load_image_rgb(image_path)
    height_px, width_px = img_np.shape[:2]

    pink_mask = find_pink_mask(img_np)
    if not np.any(pink_mask):
        raise ValueError("No pink pixels detected — check color range or image content.")

    crack_edge_y_px = np.full(width_px, np.nan)

    for x in range(width_px):
        pink_y_coords = np.where(pink_mask[:, x])[0]
        if pink_y_coords.size > 0:
            crack_edge_y_px[x] = np.max(pink_y_coords)

    # Interpolate NaNs (columns with no pink pixels)
    valid = ~np.isnan(crack_edge_y_px)
    if np.sum(valid) < 2:
        raise ValueError("Too few pink pixels found across image width.")
    crack_edge_y_px = np.interp(np.arange(width_px), np.where(valid)[0], crack_edge_y_px[valid])

    crack_edge_y_um = crack_edge_y_px * step_size
    total_height_um = height_px * step_size

    return {
        "image": os.path.basename(image_path),
        "crack_edge_y_px": crack_edge_y_px,
        "crack_edge_y_um": crack_edge_y_um,
        "height_px": height_px,
        "height_um": total_height_um
    }


def calculate_distance_from_crack_edge(image_path: str) -> dict:
    """
    Calculates distance from crack edge (in pixels and microns) for every pixel.

    Returns:
        dict with:
            - distance_px: 2D numpy array (height x width) with distance in pixels
            - distance_um: same array in microns
            - crack_edge profile and image info
    """
    img_np = load_image_rgb(image_path)
    height_px, width_px = img_np.shape[:2]


    crack_edge_data = calculate_crack_edge_per_column(image_path)
    crack_edge_y_px = crack_edge_data["crack_edge_y_px"]

    # Create a 2D array where each column contains crack_edge_y_px[x]
    crack_edge_2d = np.tile(crack_edge_y_px, (height_px, 1))

    # Create a 2D array of y pixel indices (rows)
    y_indices = np.arange(height_px).reshape((height_px, 1))

    # Distance from crack edge = pixel y - crack edge y at that x
    distance_px = y_indices - crack_edge_2d
    distance_um = distance_px * step_size

    return {
        "image": crack_edge_data["image"],
        "height_px": height_px,
        "height_um": crack_edge_data["height_um"],
        "crack_edge_y_px": crack_edge_y_px,
        "crack_edge_y_um": crack_edge_data["crack_edge_y_um"],
        "distance_px": distance_px,
        "distance_um": distance_um
    }


if __name__ == "__main__":
    image_file = r"F:\EBSD\Hole12_EBSD\LocAB\edgedetection.tif"  # Replace with your filename

    edge_result = calculate_crack_edge_per_column(image_file)
    dist_result = calculate_distance_from_crack_edge(image_file)
    
    save_path = r"C:\Users\oranm\OneDrive - Imperial College London\PHD\HSSCC\KAM Analysis\results\locab"
    os.makedirs(save_path, exist_ok=True)  # Create directory if it doesn't exist

    file_name = "distance_from_crack_edge.npy"
    full_path = os.path.join(save_path, file_name)

    np.save(full_path, dist_result)
    print(f"Saved distance_um array to: {full_path}")

    print(f"Image: {edge_result['image']}")
    print(f"Image height: {edge_result['height_px']} px = {edge_result['height_um']:.2f} µm")
    print(f"Crack edge (y) range: {np.min(edge_result['crack_edge_y_px']):.2f} px to {np.max(edge_result['crack_edge_y_px']):.2f} px")
    print(f"Distance from crack edge array shape: {dist_result['distance_px'].shape}")
