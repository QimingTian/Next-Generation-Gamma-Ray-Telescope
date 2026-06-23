"""SiPM geometry mapping mirroring DetectorConstruction.cc."""
from __future__ import annotations

import numpy as np

GRID = 18
SIPM_PER_FACE = GRID * GRID
HALF_LXE_MM = 500.0
SIPM_SPACING_MM = 50.0
SIPM_OFFSET_MM = 0.5 + 0.01  # half thickness + gap

# Face centers and local axes (u, v) in mm; copy order matches Geant4 placement.
_FACE_SPECS: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = [
    (np.array([0.0, 0.0, HALF_LXE_MM - SIPM_OFFSET_MM]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
    (np.array([0.0, 0.0, -HALF_LXE_MM + SIPM_OFFSET_MM]), np.array([1.0, 0.0, 0.0]), np.array([0.0, -1.0, 0.0])),
    (np.array([HALF_LXE_MM - SIPM_OFFSET_MM, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0])),
    (np.array([-HALF_LXE_MM + SIPM_OFFSET_MM, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    (np.array([0.0, HALF_LXE_MM - SIPM_OFFSET_MM, 0.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, -1.0])),
    (np.array([0.0, -HALF_LXE_MM + SIPM_OFFSET_MM, 0.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])),
]

# Inward normals (toward LXe center).
_FACE_NORMALS: list[np.ndarray] = [
    np.array([0.0, 0.0, -1.0]),
    np.array([0.0, 0.0, 1.0]),
    np.array([-1.0, 0.0, 0.0]),
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, -1.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
]


def face_index(sipm_id: int) -> int:
    return int(sipm_id) // SIPM_PER_FACE


def face_local_id(sipm_id: int) -> int:
    return int(sipm_id) % SIPM_PER_FACE


def face_normal(face_id: int) -> np.ndarray:
    return _FACE_NORMALS[face_id].copy()


def sipm_center_mm(sipm_id: int) -> np.ndarray:
    fid = face_index(sipm_id)
    local = face_local_id(sipm_id)
    row, col = divmod(local, GRID)
    center, u_dir, v_dir = _FACE_SPECS[fid]
    offset_u = (row - (GRID - 1) / 2.0) * SIPM_SPACING_MM
    offset_v = (col - (GRID - 1) / 2.0) * SIPM_SPACING_MM
    return center + offset_u * u_dir + offset_v * v_dir


def all_sipm_centers_mm() -> np.ndarray:
    return np.stack([sipm_center_mm(i) for i in range(SIPM_PER_FACE * 6)], axis=0)
