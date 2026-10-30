import json, urllib.request, hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).parent
manifest=json.loads((ROOT/'download_manifest.json').read_text())
def get(e):
    p=ROOT/'inputs'/e['path'];p.parent.mkdir(parents=True,exist_ok=True)
    data=urllib.request.urlopen('https://raw.githubusercontent.com/RichardBenear/ModelG/'+manifest['commit']+'/'+e['path'],timeout=60).read()
    actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if actual!=e['sha']:raise RuntimeError(e['path'])
    p.write_bytes(data)
    return dict(path=e['path'],sha=actual,size=len(data))
rows=list(ThreadPoolExecutor(max_workers=6).map(get,manifest['files']))
(ROOT/'verified_inputs.json').write_text(json.dumps(rows,indent=2))
print('Verified',len(rows),'files;',sum(x['size'] for x in rows),'bytes; commit',manifest['commit'])
