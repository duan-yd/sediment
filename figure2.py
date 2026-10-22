"""Figure 2: existing cross-level results, month labels, no new inference."""
from _support import *
from matplotlib.patches import FancyBboxPatch
def main():
    d=Delivery(2,__file__)
    scores=d.read(P+'layer_PCoA_scores.tsv','A–C','Existing 36-layer Bray–Curtis PCoA')
    per=d.read(P+'Fig3B_PERMANOVA.tsv','A–C / supplement','Original partial R2 and nine-test BH q; annotate only q<0.05')
    corr=d.read(P+'Fig3C_cross_layer_associations.tsv','D','Existing 12-core Mantel and PROTEST results')
    env=d.read(V+'Fig3D_environment_source.tsv','E','Existing bulk-chemistry db-RDA adjusted R2')
    d.evidence('03_workflows/R/46_figures2_6_curated_v13.R','Methods','Original statistics; not executed')
    assert scores.groupby('Profile').size().eq(36).all()
    fig=plt.figure(figsize=(14.2,10.2))
    fig.suptitle('Cross-level organization of microbial communities and functions',y=.972,fontsize=17,weight='bold')
    h=season_handles()+[Line2D([],[],ls='',marker=m,color='#62717D',label=l) for m,l in [('o','0–30 cm'),('^','30–60 cm'),('s','60–90 cm')]]
    fig.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,.935),ncol=7,frameon=False)
    gs=fig.add_gridspec(2,6,left=.065,right=.97,top=.832,bottom=.12,hspace=.56,wspace=.72,height_ratios=[1,1])
    annotations=[]
    for i,p in enumerate(['Taxonomy','KEGG','Nutrient']):
        ax=fig.add_subplot(gs[0,2*i:2*i+2]);z=scores[scores.Profile.eq(p)]
        for s in SEASONS:
            for dep,m in [('Surface','o'),('Middle','^'),('Deep','s')]:
                q=z[z.OriginalSeasonLabel.eq(s)&z.Depth.eq(dep)]
                ax.scatter(q.PCoA1,q.PCoA2,s=43,c=SC[s],marker=m,edgecolors='white',lw=.5,zorder=3)
        ax.axhline(0,c='#E0E7EC',lw=.65);ax.axvline(0,c='#E0E7EC',lw=.65)
        ax.set(xlabel=f'PCoA1 ({z.Axis1_pct.iloc[0]:.2f}%)',ylabel=f'PCoA2 ({z.Axis2_pct.iloc[0]:.2f}%)')
        ax.set_aspect('equal',adjustable='datalim');ax.margins(.12)
        ax.set_title(chr(65+i)+'  '+['Taxonomic PCoA','KEGG functional PCoA','Nutrient-marker PCoA'][i],loc='left',weight='bold',pad=15)
        supported=per[per.Profile.eq(p)&per.Term.isin(['SamplingPeriod','Depth'])&per.q_BH.lt(.05)]
        if len(supported):
            lines=['PERMANOVA']+[('Depth' if r.Term=='Depth' else 'Sampling month')+f' partial R² = {r.partial_R2:.3f}' for r in supported.itertuples()]
            ax.text(.02,.98,'\n'.join(lines),transform=ax.transAxes,va='top',fontsize=8)
            annotations+=lines
    ax=fig.add_subplot(gs[1,:3])
    ax.set_title('D  Cross-level correspondence',loc='left',weight='bold',pad=16)
    pairs=['Tax–KEGG','Tax–Nutrient','KEGG–Nutrient']
    lookup={(r.Layer1,r.Layer2):r for r in corr.itertuples()}
    rows=[lookup[('Taxonomy','KEGG')],lookup[('Taxonomy','Nutrient')],lookup[('KEGG','Nutrient')]]
    y=np.arange(3)[::-1]
    for i,r in enumerate(rows):
        strong=(r.Layer1,r.Layer2)==('KEGG','Nutrient')
        ax.plot([r.Mantel_r,r.Procrustes_r],[y[i],y[i]],color='#28786D' if strong else '#A9BAC4',lw=3 if strong else 1.5,zorder=1)
        ax.scatter(r.Mantel_r,y[i],s=58,c='#28786D' if strong else '#527589',marker='o',label='Mantel r' if i==0 else None,zorder=3)
        ax.scatter(r.Procrustes_r,y[i],s=58,c='white',edgecolors='#28786D' if strong else '#527589',lw=1.4,marker='o',label='Procrustes r' if i==0 else None,zorder=3)
        ax.text(max(r.Mantel_r,r.Procrustes_r)+.025,y[i]+.10,f'Mantel q={r.Mantel_q:.4f}\nPROTEST q={r.Procrustes_q:.4f}',fontsize=7.1,va='center',color='#405563')
        ax.text(-.015,y[i],pairs[i],ha='right',va='center',fontsize=9)
        ax.text(.02,y[i]-.22,f'within-month P: Mantel {r.Within_period_Mantel_P:.4f}; PROTEST {r.Within_period_Procrustes_P:.4f}',fontsize=7.2,color='#52616D')
    ax.axvline(0,c='#DCE5EA',lw=.7);ax.set(xlim=(-.25,1.02),ylim=(-.55,2.55),xlabel='Correspondence coefficient (r)',yticks=[])
    ax.grid(axis='x',c='#E2E9ED',lw=.7);ax.set_axisbelow(True)
    ax.legend(loc='upper left',frameon=False,ncol=2,fontsize=8)
    ax.text(.99,-.50,'12 equal-weight core profiles; q values are BH-adjusted.',ha='right',fontsize=7.5,color='#52616D')
    ax=fig.add_subplot(gs[1,3:]);ax.set_title('E  Explanation by measured bulk chemistry',loc='left',weight='bold',pad=16)
    for i,r in env.iterrows():
        ax.scatter(r.Adjusted_R2,2-i+.12,s=58,c='#356F8C')
        ax.scatter(r.Partial_period_adjusted_R2,2-i-.12,s=58,facecolors='white',edgecolors='#356F8C',lw=1.3)
        if r.q_BH<.05:ax.text(r.Adjusted_R2+.015,2-i+.12,'*',va='center',fontsize=16)
    ax.set(yticks=[2,1,0],yticklabels=['Taxonomy','KEGG','Nutrient markers'],xlim=(0,.46),ylim=(-.55,2.7),xlabel='Adjusted R²',xticks=[0,.1,.2,.3,.4])
    ax.grid(axis='x',c='#E2E9ED',lw=.7);ax.set_axisbelow(True);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    ax.legend(handles=[Line2D([],[],ls='',marker='o',c='#356F8C',label='Overall'),Line2D([],[],ls='',marker='o',mfc='white',mec='#356F8C',label='Month-conditioned')],loc='upper left',ncol=2,frameon=False,fontsize=9)
    fig.text(.065,.066,'A–C: 36 layers nested in 12 cores. PCoA displays structure; the original PERMANOVA results are reported in the legend and tables.',fontsize=9,color='#52616D')
    fig.text(.065,.035,'E: OM + TN + TP + pH, 12 cores. Overall KEGG q = 0.0159; month-conditioned P = 0.0694. Water level is not included.',fontsize=9,color='#52616D')
    en='Cross-level organization of microbial communities and functions. (A–C) Existing Bray–Curtis PCoA coordinates for 36 sediment layers nested within 12 cores, with February, May, August and November indicated by colour and depth by shape. Relative profiles, Lingoes correction and axis percentages are unchanged. The point clouds describe ordination structure; they are not interpreted as evidence of group separation without the corresponding PERMANOVA tests. No depth or sampling-month term passed the original nine-test BH correction. (D) Existing correspondence statistics for the three pairwise comparisons among 12 equal-weight core profiles. Filled points are Spearman Mantel r and open points are two-dimensional PROTEST r. Labels report BH-adjusted q values for the primary tests and the original within-month permutation P values. The KEGG–nutrient-marker pair has the largest coefficient in both analyses. These coefficients describe profile correspondence, not ecological interactions or activity. (E) Existing adjusted db-RDA R² for OM, TN, TP and pH, overall and conditioned on sampling month. Overall tests retain 9,999 unrestricted permutations and BH correction across three profiles; conditional tests retain 1,295 within-month permutations. Overall KEGG q = 0.0159, whereas the month-conditioned result is P = 0.0694. No hydrological model or new test is introduced.'
    cn='跨层级群落与功能组织。A–C为原有36个分层样品的Bray–Curtis PCoA，颜色为February、May、August和November，形状为深度。点云用于展示排序结构，组间差异仍以原PERMANOVA结果判断；原九项BH校正后没有深度或采样月份项达到显著。D完整展示三组层级比较的Mantel r和PROTEST r，实心点为Mantel、空心点为PROTEST，旁边列出主检验BH q和采样月份内限制置换P；KEGG与营养标记的对应系数最高。E为OM、TN、TP和pH的整体及控制采样月份后的db-RDA调整R²；总体KEGG q=0.0159，控制月份后P=0.0694，二者不可混用。未加入水位检验。'
    d.finish(fig,en,cn,'Fig.2：D由三角网络改为完整统计点图，展示三组Mantel/PROTEST系数、BH q及采样月份内限制置换P；E明确区分总体KEGG q=0.0159与控制月份P=0.0694。月份标签统一为英文月份名。',{'PCoA_points_per_layer':36,'supported_depth_month_terms':len(annotations),'new_tests':False},{'correspondence_existing_statistics.tsv':corr,'environment_existing_statistics.tsv':env,'PERMANOVA_existing_statistics.tsv':per})
if __name__=='__main__':main()
