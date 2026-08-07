#!/usr/bin/env python3
"""
Combined recession plot: both 2023 and 2024 Winooski recessions on one axis,
aligned by hours since peak. Designed for wrapfigure (~0.38\textwidth).

Usage:
    uv run --with requests --with pandas --with numpy --with matplotlib \
        analysis/recession_combined.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
from winooski_recession import (
    fetch_event_data, find_peak, fit_exponential_recession,
    SITE_NO, EVENTS, COLOR_DATA, COLOR_FIT,
)

IMAGES_DIR = REPO_ROOT / "images"
IMAGES_DIR.mkdir(parents=True, exist_ok=True)

COLOR_2023 = "#0072B2"
COLOR_2024 = "#CC79A7"

def main():
    events = {}
    for key, cfg in EVENTS.items():
        df, _ = fetch_event_data(SITE_NO, cfg["start"], cfg["end"])
        peak_time, peak_Q, _, _ = find_peak(df)
        primary_end = min(peak_time + pd.Timedelta(hours=60),
                          df.index[df.index > peak_time][-1])
        fit = fit_exponential_recession(df, peak_time, primary_end, peak_Q)
        events[key] = dict(df=df, peak_time=peak_time, peak_Q=peak_Q, fit=fit, label=cfg["label"])

    plt.rcParams.update({
        "font.size": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
    })
    fig, ax = plt.subplots(figsize=(3.2, 2.4))

    for key, color, marker in [("2023", COLOR_2023, None), ("2024", COLOR_2024, None)]:
        ev = events[key]
        df = ev["df"]
        peak_time = ev["peak_time"]
        fit = ev["fit"]

        plot_start = peak_time - pd.Timedelta(hours=12)
        plot_end = peak_time + pd.Timedelta(hours=120)
        window = df.loc[plot_start:plot_end, "discharge_cfs"]
        hours = (window.index - peak_time).total_seconds() / 3600.0

        ax.plot(hours, window.values, color=color, lw=0.8, alpha=0.85,
                label=f"{ev['label']}")

        t_line = np.linspace(0, fit["window_hours"], 80)
        Q_line = np.exp(fit["slope"] * t_line + fit["intercept"])
        ax.plot(t_line, Q_line, color=color, lw=1.2, ls="--", alpha=0.9)

        tau_txt = rf"$\tau$={fit['tau_hours']:.0f}h ($R^2$={fit['r2']:.2f})"
        if key == "2023":
            ax.text(65, 22000, tau_txt, fontsize=5.5, color=color, ha="left")
        else:
            ax.text(65, 3200, tau_txt, fontsize=5.5, color=color, ha="left")

    ax.set_yscale("log")
    ax.set_xlabel("Hours since peak", fontsize=7)
    ax.set_ylabel("Discharge (cfs)", fontsize=7)
    ax.set_xlim(-12, 120)
    ax.legend(loc="lower left", fontsize=6, frameon=False, handlelength=1.5)

    fig.tight_layout(pad=0.4)
    out_pdf = IMAGES_DIR / "recession_combined.pdf"
    out_png = IMAGES_DIR / "recession_combined.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
