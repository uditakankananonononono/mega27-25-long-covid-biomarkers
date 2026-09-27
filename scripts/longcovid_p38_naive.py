"""P38 registered sorted naive CD4 patient-level pseudobulk six-gene stress test."""
import csv,hashlib,tarfile,io,json,sys
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.stats import ttest_ind,t
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import stats
root=Path(__file__).resolve().parents[1];p=root/'data/geo/p38/GSE334857_RAW.tar';sha='a36f12d7689926e94f25c0703b017f0d58b356f2da3d0038b438b1f1e7cf9b15';assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
r=list(csv.DictReader((root/'results/longcovid_p38_samples.csv').open()));assert len(r)==25 and len({z['gsm'] for z in r})==25
selected=[z for z in r if z['fraction']=='Naive'];assert len(selected)==10 and len({z['person'] for z in selected})==10
fixed=['JAK1','CXCL8','BCL2L1','OSM','MAP3K8','STAT3'];sums={};cells={};gene_index=None
with tarfile.open(p) as tar:
 assert len(tar.getmembers())==25 and {z['tar_member'] for z in r}=={m.name for m in tar.getmembers()}
 for z in selected:
  with h5py.File(io.BytesIO(tar.extractfile(z['tar_member']).read()),'r') as f:
   m=f['matrix'];shape=m['shape'][:];names=pd.Index([b.decode() for b in m['features']['name'][:]]);assert shape[0]==len(names) and shape[1]==len(m['barcodes']) and shape[1]>=100
   if gene_index is None:gene_index=names
   else:assert names.equals(gene_index)
   data=m['data'][:];assert np.isfinite(data).all() and (data>=0).all() and np.array_equal(data,np.floor(data))
   x=sparse.csc_matrix((data,m['indices'][:],m['indptr'][:]),shape=shape)
   nonempty=np.asarray(x.sum(axis=0)).ravel()>0;assert nonempty.sum()>=100
   sums[z['person']]=np.asarray(x.sum(axis=1)).ravel().astype(float);cells[z['person']]=int(nonempty.sum())
counts=pd.DataFrame(sums,index=gene_index);dup=counts.index[counts.index.duplicated(keep=False)];counts=counts.loc[~counts.index.isin(dup)].copy();assert counts.index.is_unique and counts.to_numpy().sum()>0
cpm=counts.div(counts.sum(axis=0),axis=1)*1e6;log=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1)
a=[z['person'] for z in selected if z['status']=='PASC'];b=[z['person'] for z in selected if z['status']=='Recovered'];assert len(a)==6 and len(b)==4
out=[]
coverage={z:dict(raw_symbol_rows=int(sum(gene_index==z)),duplicate_symbol=bool(z in dup),measured=bool(z in log.index),samples_above_one_cpm=int((cpm.loc[z]>1).sum()) if z in cpm.index else None) for z in fixed}
for gene in fixed:
 if gene not in log.index:out.append(dict(gene=gene,status='missing or duplicate or filtered',up=None));continue
 x=log.loc[gene,a].to_numpy();y=log.loc[gene,b].to_numpy();d=float(x.mean()-y.mean());s1=np.var(x,ddof=1)/len(x);s0=np.var(y,ddof=1)/len(y);df=(s1+s0)**2/(s1*s1/(len(x)-1)+s0*s0/(len(y)-1));ci=list(map(float,t.interval(.95,df,loc=d,scale=np.sqrt(s1+s0))));g,v=stats.hedges_g(x.reshape(1,-1),y.reshape(1,-1));pv=float(ttest_ind(x,y,equal_var=False).pvalue)
 out.append(dict(gene=gene,status='measured',n_pasc=6,n_recovered=4,log2cpm_difference=d,welch_ci95=ci,hedges_g=float(g[0]),welch_p=pv,up=bool(d>0),bonferroni_pass=bool(d>0 and pv<.05/6)))
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334857',sha256=sha,n_records=25,n_people=10,n_pasc=6,n_recovered=4,naive_nonempty_cells=cells,duplicate_symbol_rows=len(dup),n_measurable_genes=len(log),target_coverage=coverage,genes=out,registered_descriptive_support=bool(len(out)==6 and all(z.get('bonferroni_pass') for z in out)),caveat='Six genes previously selected from GSE226260 paper; ten patients only, sorted naive CD4 single-cell pseudobulk versus whole blood; no novel discovery or comparable published benchmark.')
(root/'results/longcovid_p38_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
