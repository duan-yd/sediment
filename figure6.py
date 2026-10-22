"""Figure 6 only: reuse stored paired contrasts, three-layer profiles and water trajectories."""
from _support import *
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Rectangle
def main():
    d=Delivery(6,__file__)
    raw=d.read(H+'three_layer_profiles_252.tsv','A–C','Existing paired deep–surface and middle-curvature TPM, and three observed layers')
    means=d.read(H+'period_marker_summary_84.tsv','C','Existing three-core means at each depth')
    water=d.read(H+'water_windows.tsv','D','Existing 30-day water trajectories including sampling day')
    defs=d.read(P+'group_definitions.tsv','A–C','Fixed 21 marker identities')
    d.evidence('03_workflows/python/hydrology_three_layers_20261011_v2.py','Lineage','Source of stored contrasts/profiles; not executed')
    d.evidence('03_workflows/python/89_figures2_6_integrated_v22.py','Lineage','Predecessor conceptual Figure 6 replaced with observations')
    assert len(raw)==252 and not raw.duplicated(['Code','CoreID']).any()
    # Arithmetic verification only: the figure reads the saved contrasts, not new statistical estimates.
    assert np.allclose(raw.Delta_TPM,raw.Deep_TPM-raw.Surface_TPM)
    assert np.allclose(raw.Middle_curvature_TPM,raw.Middle_TPM-(raw.Surface_TPM+raw.Deep_TPM)/2)
    cores=[]
    for s in SEASONS:
        for pos in ['Near','Middle','Far']:
            a=raw[raw.Season.eq(s)&raw.DistanceClass.eq(pos)].CoreID.unique();assert len(a)==1;cores.append(a[0])
    codes=defs.Code.tolist()
    mats=[raw.pivot(index='Code',columns='CoreID',values=f).loc[codes,cores] for f in ['Delta_TPM','Middle_curvature_TPM']]
    maxval=max(np.abs(m.to_numpy()).max() for m in mats)
    vmax=10**np.ceil(np.log10(maxval));norm=SymLogNorm(linthresh=10,linscale=1,vmin=-vmax,vmax=vmax,base=10)
    fig=plt.figure(figsize=(14,14.4));fig.suptitle('Vertical organization of nutrient-related functional potential',y=.979,fontsize=17,weight='bold')
    gs=fig.add_gridspec(1,2,left=.265,right=.962,bottom=.515,top=.89,wspace=.19)
    for i,(matrix,title) in enumerate(zip(mats,['A  Deep − surface','B  Middle − endpoint mean'])):
        ax=fig.add_subplot(gs[i]);im=ax.imshow(matrix,aspect='auto',interpolation='nearest',cmap='RdBu_r',norm=norm)
        ax.set_title(title,loc='left',weight='bold',pad=37)
        ax.set(xticks=range(12),xticklabels=['N','M','F']*4,yticks=range(21),yticklabels=[f'{c}  {lab}' for c,lab in zip(codes,SHORT)] if i==0 else codes)
        ax.tick_params(labelsize=8.5,length=0,pad=5)
        for j,t in enumerate(ax.get_yticklabels()):t.set_color(EC[defs.iloc[j].Element.split(';')[0]])
        for k,s in enumerate(SEASONS):
            ax.add_patch(Rectangle((3*k-.5,-1.08),3,.24,fc=SC[s],clip_on=False,lw=0));ax.text(3*k+1,-1.40,s,ha='center',fontsize=9,color=SC[s],weight='bold')
        for x in [2.5,5.5,8.5]:ax.axvline(x,color='white',lw=1.7)
        for y in [3.5,11.5,16.5]:ax.axhline(y,color='white',lw=1.7)
        clean(ax)
    cax=fig.add_axes([.36,.479,.48,.012]);cb=fig.colorbar(im,cax=cax,orientation='horizontal')
    ticks=[-vmax,-100,-10,0,10,100,vmax];ticks=sorted(set(x for x in ticks if abs(x)<=vmax));cb.set_ticks(ticks);cb.set_ticklabels([f'{x:g}' for x in ticks])
    cb.set_label('TPM difference (symmetric-log colour scale; linear within ±10 TPM)',fontsize=9)
    fig.text(.04,.436,'C  Three-layer profiles of representative markers',fontsize=12,weight='bold')
    gs2=fig.add_gridspec(1,3,left=.075,right=.962,bottom=.24,top=.392,wspace=.36)
    for j,(code,title) in enumerate([('G14','Organic P hydrolysis'),('G18','Dissimilatory sulfur redox'),('G07','DNRA terminal reduction')]):
        ax=fig.add_subplot(gs2[j])
        for s in SEASONS:
            subset=raw[raw.Code.eq(code)&raw.Season.eq(s)]
            for r in subset.itertuples():ax.plot(range(3),[r.Surface_TPM,r.Middle_TPM,r.Deep_TPM],c=SC[s],lw=.75,alpha=.32)
            m=means[means.Code.eq(code)&means.Season.eq(s)].iloc[0]
            ax.plot(range(3),[m.Surface_mean,m.Middle_mean,m.Deep_mean],c=SC[s],lw=2,marker='o',ms=4)
        ax.set(xticks=[0,1,2],xticklabels=['0–30','30–60','60–90'],xlabel='Sediment depth (cm)',ylabel='Marker-pool TPM',ylim=(0,None))
        ax.set_title(code+'  '+title,fontsize=10,loc='left',weight='bold');ax.grid(axis='y',c='#E5EBEF',lw=.6);ax.set_axisbelow(True)
    fig.legend(handles=season_handles(),loc='upper center',bbox_to_anchor=(.71,.44),ncol=4,frameon=False,fontsize=9)
    fig.text(.075,.198,'Thin lines: individual cores. Thick lines: existing three-core means. Each marker has its own TPM axis.',fontsize=9,color='#53616C')
    fig.text(.04,.171,'D  Antecedent 30-day hydrological context',fontsize=12,weight='bold')
    gs3=fig.add_gridspec(1,4,left=.075,right=.962,bottom=.065,top=.144,wspace=.32)
    lo=np.floor(water.Water_level_m.min())-.15;hi=np.ceil(water.Water_level_m.max())+.15
    for i,s in enumerate(SEASONS):
        ax=fig.add_subplot(gs3[i]);q=water[water.Season.eq(s)].sort_values('Day_relative')
        assert q.Day_relative.tolist()==list(range(-30,1))
        ax.plot(q.Day_relative,q.Water_level_m,color=SC[s],lw=1.8)
        ax.scatter(q.Day_relative.iloc[-1],q.Water_level_m.iloc[-1],s=24,c=SC[s],zorder=3)
        ax.set(xlim=(-30,1),ylim=(lo,hi),xticks=[-30,-15,0],xlabel='Days before sampling')
        ax.set_title(s,color=SC[s],weight='bold',fontsize=10,pad=3)
        if i==0:ax.set_ylabel('Station water level (m)')
        ax.grid(axis='y',color='#E5EBEF',lw=.6)
    fig.text(.075,.015,'Chenglingji station; daily 08:00 observations. Four sampling-period contexts only; no water-level–marker association test.',fontsize=9,color='#53616C')
    en='Vertical organization of nutrient-related functional potential. (A) Existing paired deep-minus-surface TPM differences for all 21 groups and 12 cores: TPM(60–90 cm) − TPM(0–30 cm). (B) Existing middle-depth deviations: TPM(30–60 cm) − [TPM(0–30 cm) + TPM(60–90 cm)]/2. Each period contains separate near, middle and far cores (N/M/F), retaining among-core variation. Both panels share one symmetric-log colour scale in original TPM units, linear within ±10 TPM; values were not row standardized or altered. Red indicates a positive contrast and blue a negative contrast. A positive middle deviation need not indicate a middle-layer maximum; actual shapes are resolved by the three-layer profiles. (C) Existing profiles for G14, G18 and G07, selected as specified in the study figure plan. Thin lines are individual cores and thick lines are stored three-core period means. TPM axes differ among marker groups. (D) Existing Chenglingji-station water-level observations from day −30 to sampling day 0, using identical vertical scales. Daily values are 08:00 observations. Water level provides four campaign-level contexts, not 12 or 36 independent hydrological replicates; no hydrological correlation or causal model is fitted. Spring/Summer/Autumn/Winter retain the project-defined campaign labels. This figure is descriptive and adds no significance tests.'
    cn='营养相关功能潜力的垂向组织。（A）直接使用已存深浅差值：60–90 cm TPM−0–30 cm TPM。（B）使用已存中层偏离：30–60 cm TPM−两端TPM均值。21组×12根柱，每期保留近、中、远三根柱，展示柱间变异。两图共用原始TPM单位的对称对数色标，±10 TPM内线性；未做逐行标准化。红色为正、蓝色为负。中层偏离为正不一定等于中层最高，须结合三层实际值判断。（C）G14/G18/G07三层TPM剖面；细线是每根柱，粗线为原表已有三柱均值，各组纵轴不同。（D）城陵矶站采样前第30日至采样当天的08时逐日水位，共用纵轴。四条轨迹仅为四个采样期水文背景，不将它们扩展为12或36个水文重复，不开展水位—标记关联检验。春夏秋冬保持项目命名。'
    changes='Fig.6：移除原概念模型，改为已有观测结果。A展示深浅差值，B补充中层偏离，避免遗漏非单调结构；每期三个柱分别保留。C按用户指定展示三层真实剖面，D并排展示前30日水位。仅读取已有差值、均值和轨迹；公式检查用于核对，不运行原分析或新模型。色标使用共同对称对数尺度以容纳TPM动态范围，单位仍是TPM，未改动数值。'
    d.finish(fig,en,cn,changes,{'paired_contrasts':252,'water_points_per_period':31,'stored_contrasts_pass_arithmetic_check':True,'colour_limit_TPM':float(vmax),'linear_colour_interval_TPM':[-10,10]}, {'existing_three_layer_profiles.tsv':raw,'existing_period_means.tsv':means,'existing_water_trajectories.tsv':water})
if __name__=='__main__':main()
