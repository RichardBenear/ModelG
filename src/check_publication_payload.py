"""Local pre-push size/branch guard requested by the PI's local-raw storage policy."""
import subprocess,sys
ZERO='0'*40
def check(lines):
    objects={}
    for line in lines:
        local_ref,new,remote_ref,old=line.split()
        if new==ZERO:continue
        if local_ref.startswith('refs/heads/archive/'):
            raise ValueError('Local archive branches must not be uploaded. Use the compact publication checkout.')
        args=['git','rev-list','--objects',new]
        if old!=ZERO:args.append('^'+old)
        else:args.append('--not');args.append('--remotes')
        for row in subprocess.check_output(args,text=True).splitlines():
            parts=row.split(' ',1);objects[parts[0]]=parts[1] if len(parts)>1 else ''
    if not objects:return 0
    queries=''.join(key+'\n' for key in objects).encode();result=subprocess.run(['git','cat-file','--batch-check=%(objectname) %(objecttype) %(objectsize)'],input=queries,capture_output=True,check=True).stdout.decode();total=0
    for row in result.splitlines():
        oid,kind,size=row.split();size=int(size)
        if kind!='blob':continue
        path=objects[oid]
        if '/P1_production/' in path or path.endswith('P1_exploratory_mean_fields.npz'):raise ValueError('Bulk MG160 raw data are local-only: '+path)
        if size>5*1024**2:raise ValueError('New file exceeds the5MiB public default; retain locally or obtain explicit size approval: '+path)
        total+=size
    if total>10*1024**2:raise ValueError('New public payload exceeds10MiB. Keep bulk data local; review the publication package.')
    return total
if __name__=='__main__':
    try:total=check(sys.stdin);print(f'Publication guard: {total} bytes of new file content; within compact-publication limits.',file=sys.stderr)
    except Exception as e:print('Publication blocked: '+str(e),file=sys.stderr);sys.exit(1)
