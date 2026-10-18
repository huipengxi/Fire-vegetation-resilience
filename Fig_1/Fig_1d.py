# -*- coding: utf-8 -*-
"""

@author: 20213
"""

# Figure 1: plot directly from Source Data

import matplotlib.pyplot as plt
import pandas as pd


source_file = (
    "Fig1d_latitude_profile_0.5degree_Source_Data.csv"
)

output_png = "Fig1_latitude_profile_0.5degree.png"


# Read Source Data
data = pd.read_csv(source_file).dropna(
    subset=[
        "Latitude",
        "Mean_0.5deg",
        "CI_lower_0.5deg",
        "CI_upper_0.5deg",
    ]
)


# Plot
fig, ax = plt.subplots(figsize=(10, 10))

# Mean curve
ax.plot(
    data["Mean_0.5deg"],
    data["Latitude"],
    color="#5673E0",
    linewidth=5,
)

# 95% confidence interval
ax.fill_betweenx(
    data["Latitude"],
    data["CI_lower_0.5deg"],
    data["CI_upper_0.5deg"],
    color="#9C9C9C",
    alpha=0.30,
    linewidth=0,
)

# Zero reference line
ax.axvline(
    0,
    color="black",
    linewidth=5,
    linestyle="--",
    alpha=0.8,
)

# Axis limits and ticks
ax.set_xlim(0, 1)
ax.set_ylim(-45, 35)

ax.set_xticks([0, 0.5, 1])
ax.set_yticks([30, 15, 0, -15, -30])


# Axis style
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