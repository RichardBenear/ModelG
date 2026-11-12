"""Resolve PI-retained local raw data; never fetch or upload bulk archives."""
from pathlib import Path
import json
def local_archive(public_base):
    public_base=Path(public_base);record=public_base/'P1_local_archive.json'
    if not record.exists():return public_base
    cfg=json.loads(record.read_text(encoding='utf-8'));checkout=public_base.parents[2];candidate=(checkout.parent/cfg['archive_checkout_sibling']/cfg['path_in_checkout']).resolve()
    if checkout.parent.resolve() not in candidate.parents:raise ValueError('Local archive must remain within the workspace archive directory')
    if not (candidate/'P1_production/ledger.json').is_file():raise FileNotFoundError('Full raw data are local-only. Restore the sibling archive recorded in P1_local_archive.json.')
    return candidate
