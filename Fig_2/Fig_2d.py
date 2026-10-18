# -*- coding: utf-8 -*-
"""Plot Figure 2d from the supplied source data."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


folder = Path(__file__).resolve().parent
input_file = folder / "Fig_2d.csv"
output_png = folder / "Fig_2d.png"
output_pdf = folder / "Fig_2d.pdf"

periods = {
    "Post_1": ("Post 1 yr", "#FECC5C"),
    "Post_3": ("Post 3 yr", "#FD8D3C"),
    "Post_5": ("Post 5 yr", "#E6550D"),
    "Pre_1": ("Pre 1 yr", "#28A25A"),
}

plt.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
})


def main():
    data = pd.read_csv(input_file).sort_values("Latitude", ascending=False)

    required = ["Latitude", "Min_sample_size"]
    for code in periods:
        required.extend([
            f"{code}_Mean",
            f"{code}_CI_lower",
            f"{code}_CI_upper",
        ])

    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

    values = data[required].apply(pd.to_numeric, errors="coerce")
    data[required] = values
    mask = np.isfinite(values).all(axis=1)
    data = data.loc[mask].copy()

    latitude = data["Latitude"].to_numpy()
    sample_size = data["Min_sample_size"].to_numpy()

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(17, 10),
        sharey=True,
        gridspec_kw={"width_ratios": [4, 2.2], "wspace": 0.22},
    )

    for code, (label, color) in periods.items():
        ax1.plot(
            data[f"{code}_Mean"],
            latitude,
            color=color,
            linewidth=5,
            label=label,
        )
        ax1.fill_betweenx(
            latitude,
            data[f"{code}_CI_lower"],
            data[f"{code}_CI_upper"],
            color=color,
            alpha=0.20,
            linewidth=0,
        )

    ax1.axvline(0, color="black", linestyle="--", linewidth=5, alpha=0.6)
    ax1.set_xlim(-100, 50)
    ax1.set_xticks([-100, -50, 0, 50])
    ax1.set_ylim(-45, 35)
    ax1.set_yticks([30, 15, 0, -15, -30])

    ax2.fill_betweenx(latitude, 0, sample_size, color="#BDBDBD", alpha=0.65)
    ax2.plot(sample_size, latitude, color="#707070", linewidth=5)

    sample_max = max(200000, int(np.ceil(sample_size.max() / 200000) * 200000))
    ax2.set_xlim(0, sample_max)
    ax2.set_xticks([0, sample_max / 2, sample_max])
    ax2.set_xlabel(r"$n$ ($\times 10^3$)", fontsize=38, labelpad=18)
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1000:g}"))

    for ax in (ax1, ax2):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_linewidth(5)
        ax.tick_params(
            axis="x",
            direction="in",
            length=15,
            width=5,
            pad=18,
            labelsize=50,
        )

    ax1.tick_params(
        axis="y",
        direction="in",
        length=15,
        width=5,
        pad=18,
        labelsize=60,
    )
    ax2.tick_params(axis="y", left=False, labelleft=False)

    plt.subplots_adjust(left=0.11, right=0.98, bottom=0.17, top=0.95)
    fig.savefig(output_png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
