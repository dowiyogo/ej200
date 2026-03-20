#!/usr/bin/env python3
import argparse
from collections import defaultdict

import ROOT


def build_time_diff_hist(tree, time_branch="time_ns", event_branch="event_id", plane_branch="plane_id"):
    # Use earliest hit time per event per plane, then compute left-right.
    earliest = defaultdict(lambda: {0: None, 1: None})
    for entry in tree:
        event_id = int(getattr(entry, event_branch))
        plane_id = int(getattr(entry, plane_branch))
        t = float(getattr(entry, time_branch))
        cur = earliest[event_id][plane_id]
        if cur is None or t < cur:
            earliest[event_id][plane_id] = t

    # 10 ps = 0.01 ns bins over [-4, 4] ns => 800 bins
    h = ROOT.TH1D("h_time_diff", "Time diff (left - right);#Deltat_{L-R} [ns];Events", 800, -4.0, 4.0)
    for event_id, planes in earliest.items():
        t_left = planes[0]
        t_right = planes[1]
        if t_left is None or t_right is None:
            continue
        h.Fill(t_left - t_right)
    return h


def fill_yz_hists(tree, y_branch="y_mm", z_branch="z_mm", plane_branch="plane_id"):
    h_left = ROOT.TH2D(
        "h_yz_left",
        "Left plane hits;Y [mm];Z [mm]",
        200, -50.0, 50.0,
        200, -50.0, 50.0,
    )
    h_right = ROOT.TH2D(
        "h_yz_right",
        "Right plane hits;Y [mm];Z [mm]",
        200, -50.0, 50.0,
        200, -50.0, 50.0,
    )

    for entry in tree:
        plane_id = int(getattr(entry, plane_branch))
        y = float(getattr(entry, y_branch))
        z = float(getattr(entry, z_branch))
        if plane_id == 0:
            h_left.Fill(y, z)
        elif plane_id == 1:
            h_right.Fill(y, z)
    return h_left, h_right


def build_photons_per_event_hist(tree, event_branch="event_id", nbins_max=2000):
    counts = defaultdict(int)
    for entry in tree:
        event_id = int(getattr(entry, event_branch))
        counts[event_id] += 1

    max_n = max(counts.values()) if counts else 0
    if max_n <= nbins_max:
        nbins = max(1, max_n + 1)
        xmax = max_n + 0.5
    else:
        nbins = nbins_max
        xmax = max_n + 0.5

    h = ROOT.TH1D(
        "h_photons_per_event",
        "Detected photons per event;N_{photons} (both planes);Events",
        nbins,
        -0.5,
        xmax,
    )
    for n in counts.values():
        h.Fill(n)
    return h


def main():
    parser = argparse.ArgumentParser(description="Plot time diff and YZ hit maps from photon_hits.root")
    parser.add_argument("--input", default="photon_hits.root", help="Input ROOT file")
    parser.add_argument("--tree", default="plane_hits", help="TTree name")
    parser.add_argument("--out", default="plots.root", help="Output ROOT file for histograms")
    parser.add_argument("--png", default="plots.png", help="Output PNG for quick look")
    parser.add_argument("--pdf", default="plots.pdf", help="Output multi-page PDF")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open input file: {args.input}")
    tree = f.Get(args.tree)
    if not tree:
        raise RuntimeError(f"Cannot find TTree: {args.tree}")

    h_time = build_time_diff_hist(tree)
    h_nphot = build_photons_per_event_hist(tree)
    h_left, h_right = fill_yz_hists(tree)

    fout = ROOT.TFile(args.out, "RECREATE")
    h_time.Write()
    h_nphot.Write()
    h_left.Write()
    h_right.Write()
    fout.Close()

    # Quick PNG preview (single canvas with 3 pads)
    c = ROOT.TCanvas("c", "c", 1200, 900)
    c.Divide(2, 2)

    c.cd(1)
    h_time.Draw("hist")

    c.cd(2)
    h_left.Draw("COLZ")

    c.cd(3)
    h_right.Draw("COLZ")

    c.cd(4)
    h_nphot.Draw("hist")

    c.SaveAs(args.png)

    # Multi-page PDF: one plot per page (4 pages)
    cpdf = ROOT.TCanvas("cpdf", "cpdf", 900, 700)
    cpdf.Print(args.pdf + "[")

    h_time.Draw("hist")
    cpdf.Print(args.pdf)

    cpdf.Clear()
    h_nphot.Draw("hist")
    cpdf.Print(args.pdf)

    cpdf.Clear()
    h_left.Draw("COLZ")
    cpdf.Print(args.pdf)

    cpdf.Clear()
    h_right.Draw("COLZ")
    cpdf.Print(args.pdf)

    cpdf.Print(args.pdf + "]")


if __name__ == "__main__":
    main()
