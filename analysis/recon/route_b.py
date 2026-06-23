"""Route B: ML direction regressor using engineered SiPM features."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from recon.io import load_runs
from recon.route_a import face_centroid_3d, face_totals, truth_direction_from_event
from recon.geometry import _FACE_SPECS


def event_features(summary: pd.DataFrame, event_id: int) -> np.ndarray:
    ftot = face_totals(summary, event_id)
    total = ftot.sum()
    if total <= 0:
        return np.zeros(20, dtype=np.float32)
    fn = ftot / total
    c_top = face_centroid_3d(summary, event_id=event_id, face_id=0)
    c_bot = face_centroid_3d(summary, event_id=event_id, face_id=1)
    off = np.zeros(3)
    if c_top is not None and c_bot is not None:
        off = 0.5 * ((c_top - _FACE_SPECS[0][0]) + (c_bot - _FACE_SPECS[1][0]))
    # Pair asymmetries along X,Y,Z.
    asym = np.array([
        ftot[1] - ftot[0],
        ftot[3] - ftot[2],
        ftot[5] - ftot[4],
    ], dtype=np.float32) / total
    return np.concatenate([fn.astype(np.float32), off.astype(np.float32), asym])


def build_feature_dataset(data_dir: Path, tag: str) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    summary, events = load_runs(data_dir, tag)
    xs, ys, meta = [], [], []
    for eid in sorted(summary["EventID"].unique()):
        if eid not in events["EventID"].values:
            continue
        xs.append(event_features(summary, int(eid)))
        row = events[events["EventID"] == eid].iloc[0]
        ys.append(truth_direction_from_event(row))
        meta.append({"EventID": int(eid), "E_primary_GeV": float(row["E_primary_GeV"])})
    return np.asarray(xs, dtype=np.float32), np.asarray(ys, dtype=np.float32), pd.DataFrame(meta)


def train_and_predict(X: np.ndarray, y: np.ndarray, test_fraction: float = 0.2, seed: int = 42):
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.multioutput import MultiOutputRegressor
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_fraction, random_state=seed)
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("gbr", MultiOutputRegressor(
            GradientBoostingRegressor(n_estimators=200, max_depth=4, learning_rate=0.05, random_state=seed)
        )),
    ])
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    norms = np.linalg.norm(y_pred, axis=1, keepdims=True)
    y_pred = y_pred / np.maximum(norms, 1e-9)
    return y_test, y_pred, model


def predict_all(model, X: np.ndarray) -> np.ndarray:
    y = model.predict(X)
    norms = np.linalg.norm(y, axis=1, keepdims=True)
    return y / np.maximum(norms, 1e-9)
