"""Route A: geometric luminous-axis reconstruction from SiPM counts."""
from __future__ import annotations

import numpy as np
import pandas as pd

from recon.geometry import SIPM_PER_FACE, face_index, sipm_center_mm


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-12:
        return np.array([0.0, 0.0, -1.0])
    return v / n


def face_totals(summary: pd.DataFrame, event_id: int | None = None) -> np.ndarray:
    df = summary if event_id is None else summary[summary["EventID"] == event_id]
    totals = np.zeros(6, dtype=float)
    grouped = df.groupby("SiPMID")["PhotonCount"].sum()
    for sipm_id, count in grouped.items():
        totals[face_index(int(sipm_id))] += float(count)
    return totals


def face_centroid_3d(summary: pd.DataFrame, face_id: int, event_id: int) -> np.ndarray | None:
    df = summary[(summary["EventID"] == event_id)]
    grouped = df.groupby("SiPMID")["PhotonCount"].sum()
    positions, weights = [], []
    lo, hi = face_id * SIPM_PER_FACE, (face_id + 1) * SIPM_PER_FACE
    for sipm_id, count in grouped.items():
        sid = int(sipm_id)
        if lo <= sid < hi:
            positions.append(sipm_center_mm(sid))
            weights.append(float(count))
    if not weights or sum(weights) <= 0:
        return None
    return np.average(np.array(positions), axis=0, weights=np.array(weights))


def reconstruct_route_a(summary: pd.DataFrame, event_id: int) -> dict:
    """Return reconstructed unit direction and diagnostics."""
    df = summary[summary["EventID"] == event_id]
    if df.empty:
        return {"direction": np.array([0.0, 0.0, -1.0]), "vertex_mm": np.zeros(3), "method": "empty"}

    grouped = df.groupby("SiPMID", as_index=False)["PhotonCount"].sum()
    positions = np.array([sipm_center_mm(int(r.SiPMID)) for r in grouped.itertuples()])
    weights = grouped["PhotonCount"].to_numpy(dtype=float)
    if weights.sum() <= 0:
        return {"direction": np.array([0.0, 0.0, -1.0]), "vertex_mm": np.zeros(3), "method": "zero"}

    ftot = face_totals(summary, event_id)
    centroids = [face_centroid_3d(summary, f, event_id) for f in range(6)]
    c_top, c_bot = centroids[0], centroids[1]

    if c_top is not None and c_bot is not None:
        direction = _unit(c_bot - c_top)
        method = "chord"
    else:
        direction = np.array([0.0, 0.0, -1.0])
        method = "default"

    vertex = np.average(positions, axis=0, weights=weights)
    if c_top is not None:
        vertex = 0.85 * vertex + 0.15 * c_top

    return {
        "direction": direction,
        "vertex_mm": vertex,
        "brightest_face": int(np.argmax(ftot)),
        "face_totals": ftot.tolist(),
        "method": method,
    }


def opening_angle_deg(reco: np.ndarray, truth: np.ndarray) -> float:
    reco_u = _unit(reco)
    truth_u = _unit(truth)
    cosang = np.clip(np.dot(reco_u, truth_u), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosang)))


def truth_direction_from_event(row: pd.Series) -> np.ndarray:
    return np.array([row["DirX"], row["DirY"], row["DirZ"]], dtype=float)
