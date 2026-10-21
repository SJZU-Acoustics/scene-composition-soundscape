# Per-photograph image-metric extraction: reference for the metric definitions.
#
# The photographs are not public, so run_all.R does not call this script; the
# respondent-level metrics it produced are deposited in the Photo_metrics sheet
# of the data workbook (built into the pipeline by build/build_photo_metrics.py).
# To run it on a photograph collection:
#   python3 extract_photo_metrics.py <photo_root> <inventory.csv.gz> <out.csv>
# where the gzipped inventory lists one photograph per row with the columns
# folder, filename and stem (stem = the directional views of one panorama).
# Output: one row per readable JPEG with the metrics defined in metrics_for().
# The photographs are the fisheye-corrected segmentation inputs (about 640 px
# wide), so the metrics are comparative, not colorimetric. The respondent-level
# values average the metrics over the record's photographs (hue circularly).
import os, csv, gzip, math, sys
import numpy as np
from PIL import Image
from multiprocessing import Pool

VOL, INV, OUT = sys.argv[1:4]

MAXDIM = 768  # downsample bound; most files are already below this


def metrics_for(args):
    folder, filename, stem = args
    path = os.path.join(VOL, folder, filename) if folder != "." else os.path.join(VOL, filename)
    try:
        im = Image.open(path)
        im = im.convert("RGB")
        w0, h0 = im.size
        if max(im.size) > MAXDIM:
            im.thumbnail((MAXDIM, MAXDIM), Image.BILINEAR)
        a = np.asarray(im, dtype=np.float32) / 255.0
    except Exception as e:
        return {"folder": folder, "filename": filename, "stem": stem, "error": type(e).__name__}
    h, w, _ = a.shape
    R, G, B = a[..., 0], a[..., 1], a[..., 2]
    lum = 0.2126 * R + 0.7152 * G + 0.0722 * B
    V = a.max(axis=2)
    mn = a.min(axis=2)
    C = V - mn
    S = np.where(V > 1e-6, C / np.maximum(V, 1e-6), 0.0)
    # hue in degrees
    H = np.zeros_like(V)
    mask = C > 1e-6
    r_, g_, b_ = R[mask], G[mask], B[mask]
    v_, c_ = V[mask], C[mask]
    hh = np.zeros_like(v_)
    idx = v_ == r_
    hh[idx] = (60 * ((g_ - b_)[idx] / c_[idx])) % 360
    idx = v_ == g_
    hh[idx] = 60 * ((b_ - r_)[idx] / c_[idx]) + 120
    idx = v_ == b_
    hh[idx] = 60 * ((r_ - g_)[idx] / c_[idx]) + 240
    H[mask] = hh % 360

    def circ_mean(deg, wgt=None):
        if deg.size == 0:
            return np.nan
        rad = np.deg2rad(deg)
        if wgt is None:
            s, c = np.sin(rad).mean(), np.cos(rad).mean()
        else:
            tw = wgt.sum()
            if tw <= 0:
                return np.nan
            s, c = (np.sin(rad) * wgt).sum() / tw, (np.cos(rad) * wgt).sum() / tw
        return float(np.rad2deg(math.atan2(s, c)) % 360)

    # vegetation pixels: green hue band with real chroma
    veg = (H >= 60) & (H < 170) & (S >= 0.15) & (V >= 0.08)
    veg_frac = float(veg.mean())
    veg_hue = circ_mean(H[veg])
    veg_sat = float(S[veg].mean()) if veg.any() else np.nan
    veg_val = float(V[veg].mean()) if veg.any() else np.nan
    # excess green on chromaticity
    tot = R + G + B + 1e-6
    exg = 2 * G / tot - R / tot - B / tot
    # sky in top 40% of rows
    top = slice(0, max(1, int(0.4 * h)))
    Ht, St, Vt = H[top], S[top], V[top]
    sky_blue = (Ht >= 185) & (Ht <= 255) & (St >= 0.08) & (St <= 0.75) & (Vt >= 0.5)
    sky_white = (St < 0.12) & (Vt >= 0.75)
    sb, sw = float(sky_blue.mean()), float(sky_white.mean())
    sky_pix = sky_blue | sky_white
    warm = ((H < 60) | (H >= 300)) & (S >= 0.2)
    # colourfulness (Hasler & Suesstrunk)
    rg = R - G
    yb = 0.5 * (R + G) - B
    colf = float(np.sqrt(rg.std() ** 2 + yb.std() ** 2) + 0.3 * np.sqrt(rg.mean() ** 2 + yb.mean() ** 2))
    # edges / blur on luminance
    gy, gx = np.gradient(lum)
    gmag = np.hypot(gx, gy)
    lap = (np.roll(lum, 1, 0) + np.roll(lum, -1, 0) + np.roll(lum, 1, 1) + np.roll(lum, -1, 1) - 4 * lum)[1:-1, 1:-1]
    hist, _ = np.histogram(lum, bins=64, range=(0, 1))
    p = hist / hist.sum()
    p = p[p > 0]
    entropy = float(-(p * np.log2(p)).sum())
    return {
        "folder": folder, "filename": filename, "stem": stem, "error": "",
        "w": w0, "h": h0,
        "mean_R": float(R.mean()), "mean_G": float(G.mean()), "mean_B": float(B.mean()),
        "lum_mean": float(lum.mean()), "lum_sd": float(lum.std()),
        "sat_mean": float(S.mean()), "sat_sd": float(S.std()),
        "val_mean": float(V.mean()), "val_sd": float(V.std()),
        "hue_mean_w": circ_mean(H[mask], (S * V)[mask]),
        "veg_frac": veg_frac, "veg_hue": veg_hue, "veg_sat": veg_sat, "veg_val": veg_val,
        "exg_mean": float(exg.mean()),
        "sky_blue_frac_top": sb, "sky_white_frac_top": sw,
        "sky_frac_top": sb + sw,
        "sky_blue_ratio": (sb / (sb + sw)) if (sb + sw) > 1e-4 else np.nan,
        "sky_sat_top": float(St[sky_pix].mean()) if sky_pix.any() else np.nan,
        "warm_frac": float(warm.mean()),
        "colorfulness": colf,
        "edge_density": float((gmag > 0.08).mean()),
        "blur_lapvar": float(lap.var()),
        "clip_hi": float((V >= 0.98).mean()), "clip_lo": float((V <= 0.02).mean()),
        "entropy": entropy,
    }


def main():
    rows = []
    with gzip.open(INV, "rt") as fh:
        for r in csv.DictReader(fh):
            if r["filename"].lower().endswith((".jpg", ".jpeg")):
                rows.append((r["folder"], r["filename"], r["stem"]))
    print("images:", len(rows))
    fields = None
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with Pool(8) as pool, open(OUT, "w", newline="") as fh:
        wtr = None
        done = 0
        for res in pool.imap_unordered(metrics_for, rows, chunksize=32):
            if wtr is None:
                fields = ["folder", "filename", "stem", "error", "w", "h", "mean_R", "mean_G", "mean_B",
                          "lum_mean", "lum_sd", "sat_mean", "sat_sd", "val_mean", "val_sd", "hue_mean_w",
                          "veg_frac", "veg_hue", "veg_sat", "veg_val", "exg_mean",
                          "sky_blue_frac_top", "sky_white_frac_top", "sky_frac_top", "sky_blue_ratio",
                          "sky_sat_top", "warm_frac", "colorfulness", "edge_density", "blur_lapvar",
                          "clip_hi", "clip_lo", "entropy"]
                wtr = csv.DictWriter(fh, fieldnames=fields, restval="")
                wtr.writeheader()
            wtr.writerow({k: res.get(k, "") for k in fields})
            done += 1
            if done % 1000 == 0:
                print(done, flush=True)
    print("done", done)


if __name__ == "__main__":
    main()
