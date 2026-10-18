# -*- coding: utf-8 -*-
"""

@author: 20213
"""


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd




source_file = (
    "Fig1b_latitude_profile_0.5degree_Source_Data.csv"
)

output_png = "Fig1b_latitude_profile_0.5degree.png"


# Read Source Data
data = pd.read_csv(source_file).dropna(
    subset=[
        "Latitude",
        "Mean_0.5deg",
        "CI_lower_0.5deg",
        "CI_upper_0.5deg",
    ]
)

x = data["Mean_0.5deg"].to_numpy()
y = data["Latitude"].to_numpy()
ci_lower = data["CI_lower_0.5deg"].to_numpy()
ci_upper = data["CI_upper_0.5deg"].to_numpy()


# Prepare gradient line
points = np.column_stack(
    [x, y]
).reshape(-1, 1, 2)

segments = np.concatenate(
    [points[:-1], points[1:]],
    axis=1,
)

segment_values = (
    x[:-1] + x[1:]
) / 2


# Colour scale
colors = [
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
]

cmap = LinearSegmentedColormap.from_list(
    "custom_diverging",
    colors,
    N=256,
)

norm = plt.Normalize(
    -0.005,
    0.005,
)


# Plot
fig, ax = plt.subplots(figsize=(10, 10))

# Confidence interval
ax.fill_betweenx(
    y,
    ci_lower,
    ci_upper,
    color="#9C9C9C",
    alpha=0.30,
    linewidth=0,
    zorder=1,
)

# Gradient mean curve
line = LineCollection(
    segments,
    cmap=cmap,
    norm=norm,
    linewidth=5,
    zorder=3,
)

line.set_array(segment_values)
ax.add_collection(line)


# Colorbar
colorbar = fig.colorbar(
    line,
    ax=ax,
    orientation="vertical",
    extend="both",
    shrink=1.0,
    pad=0.04,
)

colorbar.set_ticks(
    [-0.005, 0, 0.005]
)

colorbar.ax.tick_params(
    direction="in",
    labelsize=40,
    width=4,
    length=10,
)

colorbar.outline.set_linewidth(4)

for label in colorbar.ax.get_yticklabels():
    label.set_fontname("Times New Roman")


# Zero reference line
ax.axvline(
    0,
    color="black",
    linewidth=5,
    linestyle="--",
    alpha=0.8,
    zorder=2,
)


# Axes
ax.set_xlim(-0.01, 0.01)
ax.set_ylim(-45, 35)

ax.set_xticks([-0.01, 0, 0.01])
ax.set_yticks([30, 15, 0, -15, -30])

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_linewidth(5)
ax.spines["left"].set_linewidth(5)

ax.tick_params(
    axis="both",
    direction="in",
    length=15,
    width=5,
    pad=20,
    labelsize=70,
)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontname("Times New Roman")


# Save
fig.savefig(
    output_png,
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
)


plt.show()