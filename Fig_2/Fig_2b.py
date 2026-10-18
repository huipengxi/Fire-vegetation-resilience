# -*- coding: utf-8 -*-
"""
@author: 20213
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


folder = Path(r"D:\Fig_2b")
input_file = folder / "Fig_2b.csv"
output_png = folder / "Fig_2b.png"
output_pdf = folder / "Fig_2b.pdf"

window_degrees = 0.5

columns = [
    "Latitude_negative",
    "Mean_negative",
    "SE_negative",
    "CI_lower_negative",
    "CI_upper_negative",
    "Latitude_positive",
    "Mean_positive",
    "SE_positive",
    "CI_lower_positive",
    "CI_upper_positive",
    "Latitude_sample",
    "Total_SampleSize",
]


def centred_roll(series, window, operation):
    result = getattr(
        series.rolling(window, center=True, min_periods=1),
        operation,
    )()
    left = (window - 1) // 2
    right = window - left - 1
    result.iloc[:left] = np.nan
    if right:
        result.iloc[-right:] = np.nan
    return result


data = pd.read_csv(
    input_file,
    skiprows=2,
    names=columns,
    na_values="NA",
)
data[columns] = data[columns].apply(pd.to_numeric, errors="coerce")

if not (
    np.allclose(data["Latitude_negative"], data["Latitude_positive"])
    and np.allclose(data["Latitude_negative"], data["Latitude_sample"])
):
    raise ValueError("Latitude columns do not match.")

data = data.rename(columns={"Latitude_negative": "Latitude"})
data = data.sort_values("Latitude", ascending=False).reset_index(drop=True)

latitude_step = np.median(
    np.diff(np.sort(data["Latitude"].dropna().unique()))
)
window_rows = max(1, round(window_degrees / latitude_step))

for column in [
    "Mean_negative",
    "CI_lower_negative",
    "CI_upper_negative",
    "Mean_positive",
    "CI_lower_positive",
    "CI_upper_positive",
]:
    data[column] = centred_roll(data[column], window_rows, "mean")

data["Total_SampleSize"] = centred_roll(
    data["Total_SampleSize"], window_rows, "sum"
)

plot_columns = [
    "Latitude",
    "Mean_negative",
    "CI_lower_negative",
    "CI_upper_negative",
    "Mean_positive",
    "CI_lower_positive",
    "CI_upper_positive",
    "Total_SampleSize",
]
plot_data = data[plot_columns].replace([np.inf, -np.inf], np.nan).dropna()

plt.rcParams.update({"font.family": "Times New Roman"})

fig, (ax, ax_n) = plt.subplots(
    1,
    2,
    figsize=(17, 10),
    sharey=True,
    gridspec_kw={"width_ratios": [4, 2.2], "wspace": 0.22},
)

for prefix, colour in [("negative", "#D16C23"), ("positive", "#1795AA")]:
    ax.plot(
        plot_data[f"Mean_{prefix}"],
        plot_data["Latitude"],
        color=colour,
        linewidth=5,
    )
    ax.fill_betweenx(
        plot_data["Latitude"],
        plot_data[f"CI_lower_{prefix}"],
        plot_data[f"CI_upper_{prefix}"],
        color=colour,
        alpha=0.30,
        linewidth=0,
    )

ax.axvline(0, color="black", linewidth=5, linestyle="--", alpha=0.75)
ax.set(xlim=(-150, 50), ylim=(-45, 35))
ax.set_xticks([-150, -100, -50, 0, 50])
ax.set_yticks([30, 15, 0, -15, -30])

sample_size = plot_data["Total_SampleSize"]
sample_max = max(200000, int(np.ceil(sample_size.max() / 200000) * 200000))

ax_n.fill_betweenx(
    plot_data["Latitude"],
    0,
    sample_size,
    color="#BDBDBD",
    alpha=0.65,
    linewidth=0,
)
ax_n.plot(sample_size, plot_data["Latitude"], color="#707070", linewidth=5)
ax_n.set(xlim=(0, sample_max), ylim=(-45, 35))
ax_n.set_xticks([0, sample_max / 2, sample_max])
ax_n.set_xlabel(r"$n$ ($\times 10^3$)", fontsize=38, labelpad=18)
ax_n.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1000:g}"))

for axis in (ax, ax_n):
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_linewidth(5)
    axis.tick_params(
        axis="x",
        direction="in",
        width=5,
        length=15,
        pad=18,
        labelsize=50,
    )

ax.tick_params(
    axis="y",
    direction="in",
    width=5,
    length=15,
    pad=18,
    labelsize=60,
)
ax_n.tick_params(axis="y", left=False, labelleft=False)

plt.subplots_adjust(left=0.11, right=0.98, bottom=0.17, top=0.95)
fig.savefig(output_png, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(output_pdf, bbox_inches="tight", facecolor="white")
plt.show()
