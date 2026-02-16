from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


def read_king_ids(path: Path) -> list[str]:
    ids = []
    with open(path, "rt") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            ids.append(re.split(r"\s+", s)[0])
    if not ids:
        raise ValueError(f"Empty KING id file: {path}")
    return ids


def read_king_matrix_square(path: Path, n: int) -> np.ndarray:
    rows = []
    with open(path, "rt") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rows.append(np.fromstring(s, sep=" ", dtype=np.float32))
    if len(rows) != n:
        raise ValueError(f"KING rows mismatch: {len(rows)} != {n}")
    mat = np.vstack(rows)
    if mat.shape != (n, n):
        raise ValueError(f"KING shape mismatch: {mat.shape} != {(n, n)}")
    mat = (mat + mat.T) / 2.0
    return mat


def load_meta(meta_path: Path) -> pd.DataFrame:
    m = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("NA")
    need = ["indiv_id", "tf", "cell", "algn_id", "gse", "path"]
    for c in need:
        if c not in m.columns:
            raise ValueError(f"meta must contain column: {c}")
    m = m.drop_duplicates("indiv_id").set_index("indiv_id", drop=False)
    return m


def intersect(ids: list[str], kin: np.ndarray, meta: pd.DataFrame) -> tuple[list[str], np.ndarray, pd.DataFrame]:
    meta_ids = set(meta.index.astype(str).tolist())
    keep = np.array([i in meta_ids for i in ids], dtype=bool)
    keep_n = int(keep.sum())
    if keep_n < 2:
        raise ValueError(f"Too few overlap samples: {keep_n}")
    if keep_n == len(ids):
        meta2 = meta.reindex(ids)
        if int(meta2["indiv_id"].isna().sum()) > 0:
            raise ValueError("Metadata missing for some KING IDs")
        return ids, kin, meta2

    idx = np.where(keep)[0]
    ids2 = [ids[i] for i in idx.tolist()]
    kin2 = kin[np.ix_(idx, idx)]
    meta2 = meta.reindex(ids2)
    if int(meta2["indiv_id"].isna().sum()) > 0:
        raise ValueError("Metadata missing for some kept samples")
    return ids2, kin2, meta2


def apply_floor(kin: np.ndarray, floor: float) -> np.ndarray:
    kin = kin.copy()
    kin[kin < floor] = 0.0
    kin = (kin + kin.T) / 2.0
    return kin


def kinship_to_distance(kin: np.ndarray) -> np.ndarray:
    dist = 1.0 - 2.0 * kin
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    dist[dist < 0] = 0.0
    return dist.astype(np.float32, copy=False)


def labels_to_indiv_ids(ids: list[str], labels: np.ndarray) -> pd.Series:
    df = pd.DataFrame({"old_indiv_id": ids, "lab": labels})
    keys = (
        df.groupby("lab")["old_indiv_id"]
        .min()
        .sort_values(kind="mergesort")
        .index.to_list()
    )
    lab_to_rank = {lab: i + 1 for i, lab in enumerate(keys)}
    ranks = df["lab"].map(lab_to_rank).astype(int).to_numpy()
    indiv = np.array([f"INDIV_{r:04d}" for r in ranks], dtype=object)
    return pd.Series(indiv, index=ids, name="indiv_id")


def split_multicell_clusters(
    out: pd.DataFrame,
    *,
    cluster_col: str = "indiv_id",
    cell_col: str = "cell",
    sep: str = "__CELL_",
) -> pd.DataFrame:
    out2 = out.copy()

    if cluster_col not in out2.columns:
        raise KeyError(f"{cluster_col} is missing in output dataframe")
    if cell_col not in out2.columns:
        raise KeyError(f"{cell_col} is missing in output dataframe")

    n_cells = out2.groupby(cluster_col, dropna = False)[cell_col].nunique()
    multicell = set(n_cells[n_cells > 1].index.astype(str).tolist())

    if not multicell:
        return out2

    out2[cluster_col] = out2[cluster_col].astype(str)

    mask = out2[cluster_col].isin(multicell)
    out2.loc[mask, cluster_col] = out2.loc[mask, cluster_col] + sep + out2.loc[mask, cell_col].astype(str)

    return out2


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--king", type=Path, required=True)
    p.add_argument("--king-id", type=Path, required=True)
    p.add_argument("--meta", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)

    p.add_argument("--floor", type=float, default=0.1)
    p.add_argument("--thr", type=float, default=0.8)
    p.add_argument("--method", type=str, default="average")

    p.add_argument(
        "--with-multicell-clusters",
        action = "store_true",
        help = "Do not split clusters that contain multiple cell lines; keep multicell clusters as-is"
    )

    args = p.parse_args()

    ids = read_king_ids(args.king_id)
    n = len(ids)

    meta = load_meta(args.meta)
    kin = read_king_matrix_square(args.king, n=n)

    ids, kin, meta = intersect(ids, kin, meta)

    kin = apply_floor(kin, floor=float(args.floor))
    dist = kinship_to_distance(kin)

    cond = squareform(dist, checks=False)
    z = hierarchy.linkage(cond, method=str(args.method))

    labels = hierarchy.fcluster(z, t=float(args.thr), criterion="distance")
    indiv = labels_to_indiv_ids(ids, labels)

    out = meta.loc[ids, ["indiv_id", "tf", "cell", "algn_id", "gse", "path"]].copy()
    out = out.rename(columns={"indiv_id": "old_indiv_id"})
    out["indiv_id"] = indiv.loc[ids].to_numpy()

    out = out[["old_indiv_id", "tf", "cell", "algn_id", "gse", "path", "indiv_id"]]
    out = out.reset_index(drop = True)

    if not args.with_multicell_clusters:
        out = split_multicell_clusters(out, cluster_col = "indiv_id", cell_col = "cell")

    out = out.sort_values(['indiv_id', 'algn_id'], kind = 'mergesort').reset_index(drop = True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
