"""MG100b production entry point requiring both qualification layers."""
import argparse
from mg11g5_supplement import hashes,read,EXTRA,BASE,report
from run_mg11g5 import register,production
from mg10d1_pool import atomic_json
from mg11g4_io import digest
from mg11g5 import ROOT
REG=ROOT/'results/MG11G5/registration.json'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--register',action='store_true');p.add_argument('--approval');p.add_argument('--production',action='store_true');a=p.parse_args();q=read(EXTRA/'results.json')
    if not q['passed'] or q['source_hashes']!=hashes():raise ValueError('Supplementary qualification/source block')
    if a.register:
        register(a.approval);r=read(REG);r.update(all_source_hashes=hashes(),supplement_sha256=digest(EXTRA/'results.json'));atomic_json(REG,r)
    elif a.production:
        r=read(REG)
        if r['all_source_hashes']!=hashes() or r['supplement_sha256']!=digest(EXTRA/'results.json'):raise ValueError('Registered supplemental source changed')
        production();report()
    else:p.error('Choose action')
