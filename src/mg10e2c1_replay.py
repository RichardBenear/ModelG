"""MG075 exact-prefix guard for prospective same-seed replay; no production CLI."""
from pathlib import Path
import hashlib
import numpy as np
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
class PrefixGuard:
    """Checks actual atomic writer arguments before storing a replay prefix."""
    def __init__(self,path,save):
        self.path=Path(path);self.sha256=digest(path);self.save=save;self.verified=0
        with np.load(path) as raw:self.old={key:raw[key].copy() for key in raw.files}
        self.frames=self.old['data'].shape[1]
    def __call__(self,path,**arrays):
        n=arrays['data'].shape[1]
        if n<=self.frames:
            for key in self.old:
                expected=self.old[key]
                if key=='data':expected=expected[:,:n]
                elif key=='times':expected=expected[:n]
                elif key=='stats' and n!=self.frames:continue
                if key not in arrays or not np.array_equal(arrays[key],expected):raise ValueError(f'Exact saved-prefix replay mismatch: {key}, frame{n}')
            self.verified=n
        elif self.verified!=self.frames:raise ValueError('Replay advanced before the saved prefix was verified')
        self.save(path,**arrays)
    def finish(self):
        if self.verified!=self.frames or digest(self.path)!=self.sha256:raise ValueError('Incomplete prefix verification or original changed')
        return dict(original_prefix_sha256=self.sha256,verified_frames=self.verified,exact_replay=True)
