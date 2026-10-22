"""Shared styling and read-only provenance; no analysis or figure-generation logic."""
from pathlib import Path
import argparse, hashlib, json, platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT = Path(__file__).resolve().parents[3]
P='04_results/processed/figures2_6_curated_v13/'
V='04_results/processed/figures3_5_v18/'
W='04_results/processed/figures3_5_refined_v29/'
H='outputs/hydrology_three_layers_20261011_v2_final/'
SEASONS=['Spring','Summer','Autumn','Winter']
MONTH=dict(zip(SEASONS,[2,5,8,11]))
MONTH_NAME={2:'February',5:'May',8:'August',11:'November'}
MONTH_LABEL={s:MONTH_NAME[m] for s,m in MONTH.items()}
SC=dict(zip(SEASONS,['#0072B2','#D98D00','#009E73','#CC79A7']))
EC={'C':'#8971AA','N':'#487AA1','P':'#50946C','S':'#C97153'}
SHORT=['C fixation','Asn deamidation','Peptide hydrolysis','Sulfolactate desulfonation','Ammonium assimilation','Assim. nitrate reduction','DNRA terminal reduction','Denitrification steps','Gln deamidation','Nitrate reduction','N fixation','Urea hydrolysis','Glycerophosphodiester hydrolysis','Organic P hydrolysis','Phosphate uptake','Phosphonate catabolism','PolyP metabolism','Dissimilatory sulfur redox','Sulfite oxidation (Soe)','Sulfur oxidation','Taurine desulfonation']
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
class Delivery:
    def __init__(self,number,script):
        self.n=number;self.inputs={};self.mapping=[];self.script=Path(script)
        parser=argparse.ArgumentParser(description=f'Replot Figure {number} from existing results only.')
        parser.add_argument('--out',type=Path,default=ROOT/f'05_figures/science_questions_20261013_final/Figure{number}')
        parser.add_argument('--replace',action='store_true',help='Replace generated outputs only, after archiving the previous edition.')
        args=parser.parse_args();self.out=args.out.resolve();self.replace=args.replace
        if self.out.exists() and not self.replace:raise FileExistsError(f'Refusing to overwrite: {self.out}. Archive old outputs and use --replace, or use a new --out.')
        self.evidence(self.script,'Drawing code','Independent single-figure script')
        self.evidence(Path(__file__),'Utility','Style, file reading, export and provenance only')
    def evidence(self,path,panel,purpose):
        path=Path(path)
        if not path.is_absolute():path=ROOT/path
        self.inputs[str(path)]=sha(path)
        self.mapping.append({'Figure':self.n,'Panel':panel,'Input_file':str(path),'Use':purpose})
        return path
    def read(self,path,panel,purpose,**kwargs):
        return pd.read_csv(self.evidence(path,panel,purpose),sep='\t',**kwargs)
    def finish(self,fig,legend_en,legend_cn,changes,checks=None,tables=None):
        assert all(sha(Path(p))==h for p,h in self.inputs.items()),'A source file changed during plotting.'
        self.out.mkdir(parents=True,exist_ok=self.replace)
        for ext in ['png','pdf']:fig.savefig(self.out/f'Figure_{self.n}.{ext}',dpi=400,facecolor='white')
        fig.savefig(self.out/f'Figure_{self.n}_preview.png',dpi=130,facecolor='white');plt.close(fig)
        pd.DataFrame(self.mapping).to_csv(self.out/'INPUT_FILES.tsv',sep='\t',index=False,encoding='utf-8-sig')
        for name,df in (tables or {}).items():df.to_csv(self.out/name,sep='\t',index=False)
        (self.out/'LEGEND_EN_CN.md').write_text(f'# Figure {self.n}\n\n## English\n\n{legend_en}\n\n## 中文\n\n{legend_cn}\n',encoding='utf-8')
        (self.out/'PANEL_CHANGES_CN.md').write_text(changes,encoding='utf-8')
        report={'figure':self.n,'source_SHA256':self.inputs,'sources_unchanged':all(sha(Path(p))==h for p,h in self.inputs.items()),'original_analysis_recomputed':False,'new_hypothesis_tests':False,'checks':checks or {},'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}}
        (self.out/'PROVENANCE_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps({'figure':self.n,'output':str(self.out),'sources_unchanged':report['sources_unchanged']}))
def clean(ax):
    for s in ax.spines.values():s.set_visible(False)
    ax.tick_params(length=0)
def season_handles():
    return [Line2D([],[],ls='',marker='o',color=SC[s],label=MONTH_LABEL[s],markersize=6) for s in SEASONS]
