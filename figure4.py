"""Figure 4: weighted bipartite attribution network, breadth, and classification coverage."""
from _support import *
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Patch
def main():
    d=Delivery(4,__file__)
    defs=d.read(P+'group_definitions.tsv','A–C','Existing process labels').set_index('Code')
    codes=['G07','G08','G14','G16'];gids=defs.loc[codes,'Group_ID'].tolist()
    source=d.read(V+'Fig4A_all_genus_mean_contributions.tsv','A/C','Existing full-pool mean genus contributions')
    div=d.read(V+'Fig4B_mean_layer_effective_genus_per_core.tsv','B','Existing per-core effective resolved genus numbers')
    stats=d.read(P+'Fig5_contributor_PERMANOVA.tsv','C','Existing contributor-composition PERMANOVA effect sizes and BH q values')
    for g in gids:assert np.isclose(source[source.Group_ID.eq(g)].Mean_contribution.sum(),1)
    display=[];named=set();coverage=[]
    for c,g in zip(codes,gids):
        a=source[source.Group_ID.eq(g)];un=a.loc[a.Genus.eq('Unassigned'),'Mean_contribution'].sum()
        top=a[~a.Genus.isin(['Unassigned','Other resolved'])].sort_values(['Mean_contribution','Genus'],ascending=[False,True]).head(5)
        named.update(top.Genus)
        for r in top.itertuples():display.append({'Code':c,'Group_ID':g,'Genus':r.Genus,'Mean_fraction':r.Mean_contribution})
        display.extend([{'Code':c,'Group_ID':g,'Genus':'Other resolved','Mean_fraction':1-un-top.Mean_contribution.sum()},{'Code':c,'Group_ID':g,'Genus':'Unassigned','Mean_fraction':un}])
        coverage.append({'Code':c,'Resolved_pct':100*(1-un),'Unassigned_pct':100*un})
    edges=pd.DataFrame(display);cov=pd.DataFrame(coverage)
    assert np.allclose(edges.groupby('Code').Mean_fraction.sum(),1)
    taxa=sorted(named,key=lambda x:-edges.loc[edges.Genus.eq(x),'Mean_fraction'].sum())+['Other resolved','Unassigned']
    labels={'G07':'DNRA terminal\nreduction','G08':'Denitrification\nsteps','G14':'Organic P\nhydrolysis','G16':'Phosphonate\ncatabolism'}
    colors=dict(zip(codes,['#487AA1','#7E6B9D','#50946C','#BF8550']))
    fig=plt.figure(figsize=(15,10.8));fig.suptitle('Taxonomic architecture of selected nutrient transformation processes',y=.977,fontsize=17,weight='bold')
    ax=fig.add_axes([.045,.14,.53,.70]);ax.set(xlim=(-.28,1.55),ylim=(-.04,1.05));ax.axis('off')
    ax.set_title('A  Weighted taxon–marker attribution network',loc='left',weight='bold',pad=22)
    yp=dict(zip(codes,[.88,.64,.40,.16]));yt=dict(zip(taxa,np.linspace(.97,.03,len(taxa))))
    for r in edges.itertuples():
        if r.Mean_fraction<=0:continue
        verts=[(.08,yp[r.Code]),(.39,yp[r.Code]),(.63,yt[r.Genus]),(.97,yt[r.Genus])]
        p=MPath(verts,[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4])
        ax.add_patch(PathPatch(p,fc='none',ec=colors[r.Code],lw=12*r.Mean_fraction,alpha=.60 if r.Genus not in ['Other resolved','Unassigned'] else .30))
    for c in codes:
        ax.scatter(.08,yp[c],s=95,c=colors[c],zorder=4)
        ax.text(.015,yp[c],c+'\n'+labels[c],ha='right',va='center',fontsize=10,weight='bold',color=colors[c])
    for t in taxa:
        ax.scatter(.97,yt[t],s=28,fc='white',ec='#5E7687',lw=.9,zorder=4)
        ax.text(1.02,yt[t],t,va='center',fontsize=9,fontstyle='normal' if t in ['Other resolved','Unassigned'] else 'italic',color='#253B4A')
    ax.text(.05,1.025,'Marker groups',ha='right',weight='bold',fontsize=9);ax.text(1.02,1.025,'Genus attribution',weight='bold',fontsize=9)
    ax.legend(handles=[Line2D([],[],c='#658393',lw=12*v,label=f'{100*v:g}%') for v in [.05,.15,.30]],title='Mean fraction of full marker pool',loc='upper center',bbox_to_anchor=(.52,-.03),ncol=3,frameon=False,fontsize=9,title_fontsize=9)
    bx=fig.add_axes([.70,.52,.26,.31]);bx.set_title('B  Effective contributor number',loc='left',weight='bold',pad=22)
    for i,g in enumerate(gids):
        q=div[div.Group_ID.eq(g)];assert len(q)==12
        for k,s in enumerate(SEASONS):
            v=q[q.OriginalSeasonLabel.eq(s)].Effective_genus_number.to_numpy()
            bx.scatter(v,i+np.linspace(-.21,.21,4)[k]+np.array([-.035,0,.035]),s=31,c=SC[s],edgecolors='white',lw=.4,zorder=3)
        med=q.Effective_genus_number.median();bx.plot([med,med],[i-.30,i+.30],c='black',lw=1.3,zorder=4)
    bx.set(yticks=range(4),yticklabels=codes,ylim=(3.6,-.6),xlim=(0,22),xlabel="Resolved effective genus number, exp(H′)")
    bx.grid(axis='x',c='#E5EBEF',lw=.7);bx.set_axisbelow(True);bx.spines['left'].set_visible(False);bx.tick_params(axis='y',length=0)
    bx.legend(handles=season_handles(),loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2,frameon=False,fontsize=8)
    cx=fig.add_axes([.70,.18,.26,.19]);cx.set_title('C  Sampling-period effect on contributor composition',loc='left',weight='bold',pad=18,fontsize=10)
    st=stats[stats.Analysis.eq('Genus_with_Unassigned') & stats.Term.eq('SamplingPeriod') & stats.Group_ID.isin(gids)].copy()
    st['Code']=st.Group_ID.map(dict(zip(gids,codes)));st=st.set_index('Code').loc[codes].reset_index()
    for i,r in st.iterrows():
        cx.scatter(r.partial_R2,i,s=62,c='#356F8C',zorder=3)
        qlabel=f'q={r.q_BH:.3f}' if r.q_BH>=.001 else 'q<0.001'
        cx.text(r.partial_R2+.018,i,qlabel,va='center',fontsize=8,color='#344F60')
    cx.axvline(.05,c='#B7C4CB',ls='--',lw=.9)
    cx.set(yticks=range(4),yticklabels=codes,ylim=(3.5,-.5),xlim=(0,.78),xlabel='Sampling-period partial R²')
    cx.grid(axis='x',c='#E5EBEF',lw=.7);cx.set_axisbelow(True);cx.spines['left'].set_visible(False);cx.tick_params(axis='y',length=0)
    cx.text(.02,-.38,'Dashed line: q threshold is shown for reference; q is printed beside each estimate.',transform=cx.transAxes,fontsize=7.2,color='#52616D')
    fig.text(.045,.040,'A: top five named genera per pool; remaining resolved signal and Unassigned retained. Edge widths encode attribution, not ecological interactions.',fontsize=9,color='#52616D')
    fig.text(.045,.017,'B: 12 core means per group; black ticks indicate medians. Effective genus number describes resolved breadth and evenness, not functional redundancy.',fontsize=9,color='#52616D')
    en='Taxonomic architecture of selected nutrient transformation processes. (A) Weighted bipartite graph for four prespecified representative marker pools: G07 DNRA terminal reduction, G08 denitrification steps, G14 organic-P hydrolysis and G16 phosphonate catabolism. These pools were selected to represent nitrogen, phosphorus and contrasting organic-P functions; they do not represent all 21 marker groups. The graph shows each pool’s five highest mean named-genus fractions. Other resolved and Unassigned are retained separately, and fractions sum to one within each pool. Edge width is proportional to the mean full-pool fraction. Edges represent taxonomic attribution of marker genes, not ecological interactions, co-occurrence or measured activity. (B) Existing core means of layer-level resolved exp(H′), 12 cores per pool. Unassigned was excluded and resolved fractions renormalized in the original diversity calculation. Month colours indicate sampling context and black ticks are descriptive medians; no new hypothesis test is fitted. (C) Existing genus-with-Unassigned contributor-composition PERMANOVA results for the four selected pools. Points show sampling-period partial R² and labels show BH-adjusted q values. The displayed statistics are inherited from the existing analysis and do not establish activity or causation. The previous resolved/Unassigned fractions remain available in the supplementary table.'
    cn='所选营养转化过程的分类架构。A展示G07、G08、G14和G16四个代表性标记组的加权二分网络，覆盖氮、磷及两类有机磷过程，不代表全部21组标记。每组突出平均份额最高的五个命名属，Other resolved和Unassigned单独保留；边宽表示基因标记的分类归属比例，不表示互作、共现或实际活性。B展示已解析贡献属的有效数量exp(H′)，颜色为英文月份，黑线为中位数；它反映分类贡献广度和均匀度，不是功能冗余的直接证据。C沿用已有Genus_with_Unassigned贡献者组成PERMANOVA结果，点为采样时期partial R²，旁边为BH校正q值；该统计不证明活性或因果。解析比例仍保留在补充表。'
    d.finish(fig,en,cn,'Fig.4：C由解析/未解析比例改为已有贡献者组成PERMANOVA的采样时期效应量；A明确四组为代表性过程并保留Other resolved和Unassigned，B保留有效贡献者数。',{'displayed_named_genera':len(named),'top_named_per_pool':5,'fractions_sum_to_one':True,'C_analysis':'Genus_with_Unassigned SamplingPeriod rows from existing Fig5_contributor_PERMANOVA.tsv'},{'network_display_edges.tsv':edges,'mean_resolution_coverage.tsv':cov,'existing_core_diversity.tsv':div,'existing_contributor_period_effects.tsv':st})
if __name__=='__main__':main()
