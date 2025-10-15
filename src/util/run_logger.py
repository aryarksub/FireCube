# util/run_logger.py
import os, csv, time, traceback
from datetime import datetime

LOG_FIELDS = [
    "timestamp", "fid", "fire_name", "center_lat", "center_lon",
    "start_time_iso", "hours",
    "era5_success", "pyregence_success", "landfire_success",
    "feds_success", "frp_success", "crop_success", "plots_success",
    "overall_status", "error_stage", "error_message", "elapsed_sec",
]

def _ensure(path: str):
    os.makedirs(os.path.dirname(path or "."), exist_ok=True)
    if not os.path.exists(path):
        with open(path, "w", newline="", encoding="utf-8") as f:
            csv.DictWriter(f, fieldnames=LOG_FIELDS).writeheader()

def _append(path: str, row: dict):
    _ensure(path)
    clean = {k: row.get(k, "") for k in LOG_FIELDS}
    with open(path, "a", newline="", encoding="utf-8") as f:
        csv.DictWriter(f, fieldnames=LOG_FIELDS).writerow(clean)

def safe_fire_name(fire_row, default: str) -> str:
    try:
        if fire_row is not None and not getattr(fire_row, "empty", True) and "Incid_Name" in fire_row.columns:
            return str(fire_row["Incid_Name"].values[0])
    except Exception:
        pass
    return default

def _short_trace(exc: BaseException) -> str:
    tb = traceback.format_exc(limit=2)
    return f"{exc.__class__.__name__}: {exc} | {(tb.strip().splitlines()[-1] if tb else '')}"

class FireRun:
    def __init__(self, path: str, **meta):
        self.path = path
        self.row = {k: "" for k in LOG_FIELDS}
        self.row.update({
            "timestamp": datetime.utcnow().isoformat(timespec="seconds") + "Z",
            **{k: meta.get(k, "") for k in ["fid","fire_name","center_lat","center_lon","start_time_iso","hours"]},
        })
        self._stage = ""
        self._t0 = None

    def at(self, stage: str):
        self._stage = stage

    def ok(self, fieldname: str):
        if fieldname in LOG_FIELDS:
            self.row[fieldname] = True

    def set(self, **kvs):
        for k, v in kvs.items():
            if k in LOG_FIELDS:
                self.row[k] = v

    def __enter__(self):
        self._t0 = time.time()
        return self

    def __exit__(self, et, ev, tb):
        self.row["elapsed_sec"] = round(time.time() - self._t0, 2) if self._t0 else ""
        if et is None:
            self.row["overall_status"] = self.row.get("overall_status") or "ok"
        else:
            self.row["overall_status"] = self.row.get("overall_status") or "error"
            self.row["error_stage"] = self._stage or self.row.get("error_stage", "")
            self.row["error_message"] = _short_trace(ev)
        _append(self.path, self.row)
        return False  # propagate exceptions

def fire_run(path: str, **meta) -> FireRun:
    return FireRun(path, **meta)
