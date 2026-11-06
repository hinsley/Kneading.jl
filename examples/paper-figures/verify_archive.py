"""Verify the exact source snapshot without loading numerical dependencies."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
manifest = json.loads((root / "SOURCE-MANIFEST.json").read_text())
for entry in manifest["files"]:
    data = (root / entry["path"]).read_bytes()
    assert len(data) == entry["bytes"], entry["path"]
    assert hashlib.sha256(data).hexdigest() == entry["sha256"], entry["path"]
print(f"Verified {len(manifest['files'])} archived files against original hashes.")
