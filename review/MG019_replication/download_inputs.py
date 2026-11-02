import json, urllib.request, hashlib, argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).parent
parser=argparse.ArgumentParser()
parser.add_argument('--manifest',default='download_manifest.json')
parser.add_argument('--output',default='inputs')
parser.add_argument('--verification',default='verified_inputs.json')
args=parser.parse_args()
manifest=json.loads((ROOT/args.manifest).read_text())
def get(e):
    p=ROOT/args.output/e['path'];p.parent.mkdir(parents=True,exist_ok=True)
    data=urllib.request.urlopen('https://raw.githubusercontent.com/RichardBenear/ModelG/'+manifest['commit']+'/'+e['path'],timeout=60).read()
    actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if actual!=e['sha']:raise RuntimeError(e['path'])
    p.write_bytes(data)
    return dict(path=e['path'],sha=actual,size=len(data))
rows=list(ThreadPoolExecutor(max_workers=6).map(get,manifest['files']))
(ROOT/args.verification).write_text(json.dumps(rows,indent=2))
print('Verified',len(rows),'files;',sum(x['size'] for x in rows),'bytes; commit',manifest['commit'])
