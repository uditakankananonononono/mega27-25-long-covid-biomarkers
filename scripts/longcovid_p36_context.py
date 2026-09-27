"""P36 locked published-gene source-context check, one deposited title-token per person."""
from pathlib import Path
import sys,hashlib,json,csv
import numpy as np,pandas as pd
from scipy.stats import ttest_ind,t
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import stats,geo
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p36'
sha={'GSE226260_AdditionalSamples.rawCounts.csv.gz':'7edfd3e993143fdf439ab3ed59f14708b2813fe5afa7e842ec855eddf122761a','GSE226260-GPL34284_series_matrix.txt.gz':'816dc4b1cc4ee6ad17d2d59b4d0d5d9911f82c5210efbbc303be2c32a06c8eb6'}
for n,digest in sha.items():assert hashlib.sha256((folder/n).read_bytes()).hexdigest()==digest
rows=list(csv.DictReader((root/'results/longcovid_p36_samples.csv').open()));assert len(rows)==103 and len({r['title'] for r in rows})==103 and len({r['gsm'] for r in rows})==103
selected=[r for r in rows if r['selected']=='1'];assert len(selected)==37 and {s:sum(r['status']==s for r in selected) for s in ('PASC','Recovered')}=={'PASC':17,'Recovered':20}
assert len({r['person_token'] for r in selected})==37
for p in {r['person_token'] for r in rows}:
 group=[r for r in rows if r['person_token']==p];assert len({r['status'] for r in group})==1
 if p not in {r['person_token'] for r in selected}:assert all(r['window']=='Acute' for r in group)
x=pd.read_csv(folder/'GSE226260_AdditionalSamples.rawCounts.csv.gz',index_col=0);assert x.index.is_unique and set(x.columns)=={r['title'] for r in rows};assert x.index.str.fullmatch(r'ENSG\d+').all();assert np.isfinite(x.values).all() and (x.values>=0).all() and np.array_equal(x.values,np.floor(x.values))
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)];mapping=dict(zip(h.ensembl_gene_id,h.symbol))
x=x[[r['title'] for r in selected]];x['gene']=x.index.map(mapping);expr=x.dropna(subset=['gene']).groupby('gene').sum();assert expr.index.is_unique and (expr.sum(axis=0)>0).all()
cpm=expr.div(expr.sum(axis=0),axis=1)*1e6;log=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1)
a=[r['title'] for r in selected if r['status']=='PASC'];b=[r['title'] for r in selected if r['status']=='Recovered'];g,v=stats.hedges_g(log[a].to_numpy(),log[b].to_numpy());eff=pd.Series(g,index=log.index)
fixed=['JAK1','CXCL8','BCL2L1','OSM','MAP3K8','STAT3'];obs=[z for z in fixed if z in log.index];rng=np.random.default_rng(20260925);pool=np.array(log.index.difference(fixed));assert len(pool)>1000
null=[]
for i in range(10000):null.append(int((eff.loc[rng.choice(pool,len(obs),replace=False)]>0).sum()))
results=[]
for z in fixed:
 if z not in obs:results.append(dict(gene=z,status='unmeasurable',up=None,welch_p=None));continue
 x1=log.loc[z,a].to_numpy();x0=log.loc[z,b].to_numpy();delta=float(x1.mean()-x0.mean());s1=np.var(x1,ddof=1)/len(x1);s0=np.var(x0,ddof=1)/len(x0);df=(s1+s0)**2/(s1*s1/(len(x1)-1)+s0*s0/(len(x0)-1));ci=t.interval(.95,df,loc=delta,scale=np.sqrt(s1+s0));p=float(ttest_ind(x1,x0,equal_var=False).pvalue)
 results.append(dict(gene=z,status='measured',n_pasc=17,n_recovered=20,pasc_mean=float(x1.mean()),recovered_mean=float(x0.mean()),log2cpm_difference=delta,welch_ci95=list(map(float,ci)),hedges_g=float(eff[z]),welch_p=p,up=bool(delta>0),bonferroni_pass=bool(delta>0 and p<.05/6)))
count=sum(z.get('up')==True for z in results);pv=float((1+sum(q>=count for q in null))/10001);rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226260',source_paper='https://www.nature.com/articles/s41590-025-02353-x',sha256=sha,n_source_samples=103,n_selected_people=37,n_pasc=17,n_recovered=20,n_mapped_symbols=len(expr),n_measurable_symbols=len(log),ambiguous_ensembl_ids=len(amb),n_fixed_measured=len(obs),up_signs=count,null_mean=float(np.mean(null)),empirical_p=pv,published_gene_context_support=bool(len(obs)==6 and count==6 and pv<.01 and any(z.get('bonferroni_pass') for z in results)),genes=results,limitation='Source paper selected these genes from related data; one sample per inferred title token, no independently supplied person crosswalk, windows mixed. Not independent validation or discovery.')
(root/'results/longcovid_p36_result.json').write_text(json.dumps(rec,indent=2)+'\n');pd.DataFrame({'random_up_genes':null}).to_csv(root/'results/longcovid_p36_null.csv.gz',index=False);print(json.dumps(rec,indent=2))
