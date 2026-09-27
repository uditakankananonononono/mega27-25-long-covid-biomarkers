"""Frozen P37 six-gene cross-assay long-COVID stress test."""
import csv,hashlib,json
from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import ttest_ind,t
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import stats
root=Path(__file__).resolve().parents[1];p=root/'data/geo/p37/GSE275334_File_1_Normalised.xlsx';sha='9db33a05d0becbc51693cd392b03da5a1def06b0b36d3160e28a04ccfea7dfc5';assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
r=list(csv.DictReader((root/'results/longcovid_p37_samples.csv').open()));assert len(r)==47 and len({a['gsm'] for a in r})==47 and len({a['column'] for a in r})==47
assert all(a['column'].startswith({'Healthy control':'HC','Long COVID':'LC','ME/CFS':'ME'}[a['disease']]) for a in r)
d=pd.read_excel(p,sheet_name='Normalised counts');assert d.Symbol.is_unique and d.Symbol.notna().all() and set(d.columns[4:])=={a['column'] for a in r}
fixed=['JAK1','CXCL8','BCL2L1','OSM','MAP3K8','STAT3'];case=[a['column'] for a in r if a['disease']=='Long COVID'];ctrl=[a['column'] for a in r if a['disease']=='Healthy control'];assert len(case)==15 and len(ctrl)==18
out=[]
for gene in fixed:
 z=d.loc[d.Symbol==gene];assert len(z)==1
 x=z[case].iloc[0].to_numpy(dtype=float);y=z[ctrl].iloc[0].to_numpy(dtype=float)
 assert np.isfinite(x).all() and np.isfinite(y).all()
 diff=float(x.mean()-y.mean());s1=np.var(x,ddof=1)/len(x);s0=np.var(y,ddof=1)/len(y);df=(s1+s0)**2/(s1*s1/(len(x)-1)+s0*s0/(len(y)-1));ci=list(map(float,t.interval(.95,df,loc=diff,scale=np.sqrt(s1+s0))));g,v=stats.hedges_g(x.reshape(1,-1),y.reshape(1,-1));pv=float(ttest_ind(x,y,equal_var=False).pvalue)
 out.append(dict(gene=gene,n_lc=len(x),n_healthy=len(y),lc_mean=float(x.mean()),healthy_mean=float(y.mean()),normalized_difference=diff,welch_ci95=ci,hedges_g=float(g[0]),welch_p=pv,up=bool(diff>0),bonferroni_pass=bool(diff>0 and pv<.05/6)))
record=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275334',published_study='https://doi.org/10.1172/jci.insight.183810',workbook_sha256=sha,n_records=47,n_lc=15,n_healthy=18,n_mecfs_context=14,selected_genes=out,registered_cross_assay_support=bool(all(a['bonferroni_pass'] for a in out)),caveat='Source genes post-selected from GSE226260; PBMC NanoString LC versus healthy is not whole-blood RNA-seq PASC versus infected recovered. No novel biomarker or benchmark comparison.')
(root/'results/longcovid_p37_result.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
