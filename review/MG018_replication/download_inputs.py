"""Download public review inputs at one pinned commit; verify every Git blob."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request, json, hashlib
ROOT=Path(__file__).resolve().parent
m=json.loads((ROOT/'download_manifest.json').read_text())
def get(e):
    p=ROOT/'inputs'/e['path']
    data=p.read_bytes() if p.exists() else urllib.request.urlopen('https://raw.githubusercontent.com/RichardBenear/ModelG/'+m['commit']+'/'+e['path'], timeout=90).read()
    sha=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    assert sha==e['sha'],(e['path'],sha,e['sha'])
    p.parent.mkdir(parents=True,exist_ok=True); p.write_bytes(data)
    return dict(path=e['path'],sha=sha,size=len(data))
rows=list(ThreadPoolExecutor(max_workers=6).map(get,m['files']))
(ROOT/'verified_inputs.json').write_text(json.dumps(rows,indent=2)+'\n')
print('Verified',len(rows),'files;',sum(e['size'] for e in rows),'bytes')
