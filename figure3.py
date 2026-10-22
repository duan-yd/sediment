"""Figure 3: process-specific and depth-dependent organization of nutrient markers.

The figure uses existing core summaries and paired depth values.  It changes only
the display: the deep-surface and middle-layer contrasts are combined into one
panel; no new statistical test is run here.
"""
from _support import *
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Rectangle


def main():
    d = Delivery(3, __file__)
    defs = d.read(P + 'group_definitions.tsv', 'A', 'Fixed 21 marker groups')
    z = d.read(V + 'Fig3C_row_zscore_core_mean_TPM.tsv', 'A', 'Existing core-mean row-zscore matrix').set_index('Group_ID').loc[defs.Group_ID]
    cm = d.read(V + 'core_display_annotations.tsv', 'A', 'Existing core/month/order').set_index('CoreID').loc[z.columns]
    raw = d.read(H + 'three_layer_profiles_252.tsv', 'B-D', 'Existing paired contrasts and three-depth TPM')
    means = d.read(H + 'period_marker_summary_84.tsv', 'D', 'Existing three-core means for depth profiles')
    assert z.shape == (21, 12) and len(raw) == 252
    assert np.allclose(raw.Delta_TPM, raw.Deep_TPM - raw.Surface_TPM)
    assert np.allclose(raw.Middle_curvature_TPM, raw.Middle_TPM - (raw.Surface_TPM + raw.Deep_TPM) / 2)

    # Three-core aggregation is for display only; all inferential results remain unchanged.
    m = raw.groupby(['Code', 'Season'], sort=False)[['Delta_TPM', 'Middle_curvature_TPM']].mean().reset_index()
    m['Month'] = m.Season.map(MONTH)
    matrices = [m.pivot(index='Code', columns='Season', values=v).loc[defs.Code, SEASONS] for v in ['Delta_TPM', 'Middle_curvature_TPM']]
    # Interleave contrasts within each month: February D-S, February M-dev, ...
    contrast = np.column_stack([matrices[k].iloc[:, i].to_numpy() for i in range(4) for k in range(2)])
    vmax = 10 ** np.ceil(np.log10(max(abs(a.to_numpy()).max() for a in matrices)))
    norm = SymLogNorm(linthresh=10, vmin=-vmax, vmax=vmax, base=10)

    fig = plt.figure(figsize=(16.4, 13.8))
    fig.suptitle('Process-specific and depth-dependent organization\nof nutrient-related functional potential', y=.982, fontsize=18, weight='bold')
    gs = fig.add_gridspec(1, 2, left=.20, right=.972, top=.875, bottom=.38, width_ratios=[1.35, 1.0], wspace=.28)
    axa, axb = [fig.add_subplot(gs[i]) for i in range(2)]

    lim = np.ceil(abs(z.to_numpy()).max() * 10) / 10
    im0 = axa.imshow(z, aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim, interpolation='nearest')
    axa.set(xticks=range(12), xticklabels=['N', 'M', 'F'] * 4, yticks=range(21), yticklabels=[f'{c}  {l}' for c, l in zip(defs.Code, SHORT)])
    axa.tick_params(axis='y', labelsize=8)
    axa.set_title('A  Marker profiles across 12 sediment cores', loc='left', weight='bold', pad=44)
    for i, t in enumerate(axa.get_yticklabels()):
        t.set_color(EC[defs.iloc[i].Element.split(';')[0]])
    for k, s in enumerate(SEASONS):
        axa.text(3 * k + 1, -1.34, MONTH_LABEL[s], ha='center', color=SC[s], weight='bold', fontsize=10)
        axa.add_patch(Rectangle((3 * k - .5, -1.03), 3, .21, color=SC[s], clip_on=False, lw=0))
    for x in [2.5, 5.5, 8.5]:
        axa.axvline(x, color='white', lw=1.8)
    for lo, hi, e in [(0, 4, 'C'), (4, 12, 'N'), (12, 17, 'P'), (17, 21, 'S')]:
        axa.add_patch(Rectangle((-.65, lo - .5), .10, hi - lo, fc=EC[e], clip_on=False, lw=0))

    imb = axb.imshow(contrast, aspect='auto', cmap='RdBu_r', norm=norm, interpolation='nearest')
    axb.set(xticks=range(8), xticklabels=['D−S', 'M dev'] * 4, yticks=range(21), yticklabels=[])
    axb.tick_params(axis='x', labelrotation=35, labelsize=8)
    axb.set_title('B  Depth contrasts across nutrient markers', loc='left', weight='bold', pad=44)
    axb.set_xlabel('')
    for k, s in enumerate(SEASONS):
        axb.add_patch(Rectangle((2 * k - .5, -1.03), 2, .21, color=SC[s], clip_on=False, lw=0))
        axb.text(2 * k + .5, -1.34, MONTH_LABEL[s], ha='center', color=SC[s], weight='bold', fontsize=9)
    for ax in [axa, axb]:
        clean(ax)
        for y in [3.5, 11.5, 16.5]:
            ax.axhline(y, color='white', lw=1.8)

    cax = fig.add_axes([.29, .335, .28, .014])
    cb = fig.colorbar(im0, cax=cax, orientation='horizontal')
    cb.set_label('Row z-score of core-mean TPM', fontsize=9)
    cax = fig.add_axes([.68, .335, .25, .014])
    cb = fig.colorbar(imb, cax=cax, orientation='horizontal')
    ticks = sorted(set([-vmax, -100, -10, 0, 10, 100, vmax])); ticks = [v for v in ticks if abs(v) <= vmax]
    cb.set_ticks(ticks); cb.set_ticklabels([f'{x:g}' for x in ticks]); cb.set_label('TPM difference; symmetric-log colour scale', fontsize=9)
    fig.text(.20, .296, 'A: N / M / F = Near / Middle / Far.  B: D−S = 60–90 cm − 0–30 cm; M dev = 30–60 cm − mean of the two end layers. Positive values indicate higher representation in the named contrast. B values are three-core mean TPM differences.', fontsize=9, color='#52616D')

    fig.text(.075, .255, 'C  Representative three-depth profiles', fontsize=13, weight='bold')
    fig.legend(handles=season_handles(), loc='center', bbox_to_anchor=(.73, .257), ncol=4, frameon=False, fontsize=10)
    gs2 = fig.add_gridspec(1, 3, left=.075, right=.967, bottom=.084, top=.215, wspace=.35)
    markers = ['o', '^', 's', 'D']
    for j, (code, title) in enumerate([('G07', 'DNRA terminal reduction'), ('G14', 'Organic P hydrolysis'), ('G18', 'Dissimilatory sulfur redox')]):
        ax = fig.add_subplot(gs2[j])
        for k, s in enumerate(SEASONS):
            a = raw[raw.Code.eq(code) & raw.Season.eq(s)]
            for r in a.itertuples():
                ax.plot([0, 1, 2], [r.Surface_TPM, r.Middle_TPM, r.Deep_TPM], color=SC[s], lw=.7, alpha=.27)
            r = means[means.Code.eq(code) & means.Season.eq(s)].iloc[0]
            ax.plot([0, 1, 2], [r.Surface_mean, r.Middle_mean, r.Deep_mean], color=SC[s], marker=markers[k], ms=4.5, lw=2)
        ax.set(xticks=[0, 1, 2], xticklabels=['0–30', '30–60', '60–90'], xlabel='Sediment depth (cm)', ylabel='Marker-pool TPM', ylim=(0, None))
        ax.set_title(code + '  ' + title, loc='left', fontsize=11, weight='bold'); ax.grid(axis='y', c='#E5EBEF', lw=.6); ax.set_axisbelow(True)
    fig.text(.075, .027, 'Thin lines: individual cores; thick lines: three-core means. Profiles use marker-specific TPM axes. G07, G14 and G18 were prespecified to represent nitrogen, phosphorus and sulfur processes, respectively; their selection was not based on curve appearance.', fontsize=9, color='#52616D')

    en = 'Process-specific and depth-dependent organization of nutrient-related functional potential. (A) Existing row-standardized TPM profiles for 21 marker groups and 12 cores. Columns represent the three depth layers averaged within each core; colours represent within-group deviations, not cross-group absolute abundance. (B) A combined display of two saved depth contrasts: deep minus surface and the middle-layer deviation from the mean of the two end layers. Values are three-core means for each sampling month; individual-core values are supplied separately. The two contrasts share a symmetric-log TPM scale and are not row-normalized. (C) G07, G14 and G18 were prespecified to span nitrogen, phosphorus and sulfur processes. Thin lines represent individual cores and thick lines the existing three-core means. Marker-specific y-axis ranges mean line heights are not comparable across panels. No new significance test or hydrological association was calculated.'
    cn = '营养相关功能潜力的过程特异性与垂向组织。A为21组标记在12根柱中的行标准化TPM图，颜色表示每组自身的相对偏离，不能跨组比较绝对丰度。B将两个已保存的垂向差值合并展示：深层减表层，以及中层相对于浅、深层均值的偏离；数值为每月三根柱的均值，逐柱值另行提供。两个差值共用对称对数TPM色标，未进行行标准化。C预先选择G07、G14和G18分别代表氮、磷和硫过程；细线为单柱，粗线为三柱均值。各标记组使用独立纵轴，不能比较不同面板的曲线高度。未新增显著性检验或水文关联检验。'
    d.finish(fig, en, cn, 'Fig.3：将深层—表层和中层偏离合并为同一对比面板；说明A的行标准化颜色与B的TPM差值含义不同；统一使用英文月份标签。', {'cores': 12, 'groups': 21, 'cores_per_month': 3, 'new_tests': False}, {'marker_dictionary.tsv': defs, 'vertical_contrast_display.tsv': pd.DataFrame(contrast, index=defs.Code, columns=[f'{MONTH_LABEL[s]}_{k}' for s in SEASONS for k in ['Deep_minus_surface', 'Middle_deviation']]).reset_index().rename(columns={'index': 'Code'}), 'vertical_core_source_values.tsv': raw, 'existing_depth_profile_means.tsv': means})


if __name__ == '__main__':
    main()
