"""Figure 5: rectangular genome-resolved architecture of nutrient markers."""
from _support import *
import ast
from matplotlib.patches import Rectangle, Patch


def main():
    d = Delivery(5, __file__)
    defs = d.read(P + 'group_definitions.tsv', 'A-B', 'Fixed marker groups; mixed-element definitions retained')
    binary = d.read(W + 'Figure5A_binary_131x21.tsv', 'A', 'Existing binary accepted-marker matrix').set_index('MAG')
    order = d.read(W + 'Figure5A_row_order_taxonomy.tsv', 'A', 'Existing taxonomy and tree-tip order').set_index('MAG').loc[binary.index]
    pc = d.read(W + 'Figure5A_phylum_colors.tsv', 'A', 'Existing phylum palette').set_index('Phylum').Color.to_dict()
    link = d.read('04_results/processed/figure5_abundance_v30/Figure5B_all255_linkage.tsv', 'C', 'Existing full-catalogue linkage fractions')
    genes = d.read(P + 'accepted_21_gene_taxonomy_MAG_TPM.tsv.gz', 'C', 'Existing linked accepted gene count', usecols=['Gene_ID', 'MAG'])
    cat = d.read(P + 'MAG_catalogue.tsv', 'A', '255-MAG catalogue')

    utility = d.evidence('03_workflows/python/98_figures3_5_refined_v29.py', 'A tree geometry', 'Only read existing parse/prune/walk utility definitions')
    module = ast.parse(utility.read_text(encoding='utf-8'))
    fs = [n for n in module.body if isinstance(n, ast.FunctionDef) and n.name in {'parse', 'prune', 'walk'}]
    ns = {'np': np, 'inputs': set()}
    exec(compile(ast.Module(body=fs, type_ignores=[]), str(utility), 'exec'), ns)
    parse, prune, walk = [ns[k] for k in ['parse', 'prune', 'walk']]

    trees, names = [], []
    for dom, fn, expected in [('Bacteria', 'MAG.unrooted.tree', 227), ('Archaea', 'MAG.unrooted.Archaea.tree', 28)]:
        path = d.evidence('02_audit/server_audit_package/result/gtdbtk/tree/' + fn, 'A', 'Existing local GTDB tree; no new inference')
        root = parse(path)
        assert sum(not n['children'] for n in walk(root)) == expected
        keep = set(order.index[order.Domain.eq(dom)])
        root = prune(root, keep)
        leaves = [n['name'] for n in walk(root) if not n['children']]
        assert set(leaves) == keep
        trees.append((dom, root))
        names.extend(order.loc[order.Domain.eq(dom)].sort_values('Plot_y').index.tolist())
    assert names == order.sort_values('Plot_y').index.tolist()
    assert len(cat) == 255 and genes.loc[genes.MAG.isin(cat.MAG), 'Gene_ID'].nunique() == 476
    assert binary.shape == (131, 21) and binary.sum(axis=1).gt(0).all()

    # Mutually exclusive C/N/P/S combinations from the fixed marker catalogue.
    e = pd.DataFrame(index=binary.index)
    for element in ['C', 'N', 'P', 'S']:
        cols = defs.loc[defs.Element.str.split(';').map(lambda a: element in a), 'Code']
        e[element] = binary[cols].any(axis=1).astype(int)
    combo = e.apply(lambda row: tuple(row.astype(int)), axis=1)
    counts = combo.value_counts(); combinations = counts.index.tolist()
    combo_table = pd.DataFrame([
        dict(Combination='+'.join(el for el, v in zip('CNPS', c) if v), MAG_count=int(counts[c]), **dict(zip('CNPS', c)))
        for c in combinations
    ])
    assert combo_table.MAG_count.sum() == 131

    fig = plt.figure(figsize=(17.0, 13.6))
    fig.suptitle('Genome-resolved architecture of nutrient-related functional potential', y=.978, fontsize=18, weight='bold')
    fig.text(.965, .948, '255 recovered MAGs  →  131 marker-bearing MAGs  →  476 linked accepted genes', ha='right', fontsize=11, color='#536675')
    fig.text(.025, .915, 'A  GTDB phylogeny, taxonomy and G01–G21 marker presence', weight='bold', fontsize=13)

    # Main rectangular tree and aligned matrix.
    ax = fig.add_axes([.025, .245, .645, .65])
    ax.set_xlim(0, 1.0); ax.set_ylim(133, -4); ax.axis('off')
    y_map = order['Plot_y'].astype(float).to_dict()
    maxx = max(n['x'] for _, t in trees for n in walk(t))
    tree_left, tree_right = .015, .285
    id_x, phylum_x, matrix_left, matrix_right = .295, .405, .44, .995

    def assign_y(n):
        if n['children']:
            for ch in n['children']: assign_y(ch)
            n['y'] = float(np.mean([ch['y'] for ch in n['children']]))
        else: n['y'] = float(y_map[n['name']])
        n['xplot'] = tree_left + (n['x'] / maxx) * (tree_right - tree_left)

    def draw_tree(n):
        if not n['children']: return
        ys = [ch['y'] for ch in n['children']]
        ax.plot([n['xplot'], n['xplot']], [min(ys), max(ys)], color='#788D9A', lw=.45, zorder=1)
        for ch in n['children']:
            ax.plot([n['xplot'], ch['xplot']], [ch['y'], ch['y']], color='#788D9A', lw=.45, zorder=1)
            draw_tree(ch)

    for _, root in trees: assign_y(root); draw_tree(root)

    cellw = (matrix_right - matrix_left) / 21
    for mag in names:
        row_y = y_map[mag]; ph = order.loc[mag, 'Phylum']
        ax.add_patch(Rectangle((phylum_x, row_y - .42), .025, .84, facecolor=pc.get(ph, '#A0ADB5'), edgecolor='none'))
        ax.text(id_x, row_y, mag, va='center', ha='left', fontsize=4.25, color='#334E5E')
        for j, code in enumerate(defs.Code):
            elem = defs.iloc[j].Element.split(';')[0]
            col = EC[elem] if int(binary.loc[mag, code]) else '#EEF2F4'
            ax.add_patch(Rectangle((matrix_left + j * cellw, row_y - .42), cellw * .93, .84, facecolor=col, edgecolor='white', lw=.18))

    for start, end, elem in [(0, 4, 'C'), (4, 12, 'N'), (12, 17, 'P'), (17, 21, 'S')]:
        x0 = matrix_left + start * cellw; x1 = matrix_left + end * cellw
        ax.add_patch(Rectangle((x0, -2.95), x1 - x0, .36, facecolor=EC[elem], alpha=.85, clip_on=False, edgecolor='none'))
        ax.text((x0 + x1) / 2, -2.77, {'C': 'Carbon', 'N': 'Nitrogen', 'P': 'Phosphorus', 'S': 'Sulfur'}[elem], ha='center', va='center', fontsize=8, color='white', weight='bold')
    for j, code in enumerate(defs.Code):
        ax.text(matrix_left + (j + .46) * cellw, -1.95, code, rotation=55, ha='left', va='bottom', fontsize=6.1, color='#2E4858')
    for x in [matrix_left + 4 * cellw, matrix_left + 12 * cellw, matrix_left + 17 * cellw]:
        ax.plot([x, x], [-3.0, 132.2], color='#8A9AA4', lw=.75, zorder=2)
    ax.text(tree_left, -2.0, 'GTDB tree', ha='left', va='bottom', fontsize=8.5, color='#52616D')
    ax.text(phylum_x + .012, -2.0, 'Phylum', ha='center', va='bottom', fontsize=8.5, color='#52616D')
    ax.text(.0, 61, 'Bacteria', rotation=90, ha='center', va='center', fontsize=8.5, color='#536675')
    ax.text(.0, 127.5, 'Archaea', rotation=90, ha='center', va='center', fontsize=8.5, color='#536675')
    ax.axhline(123.3, color='#BAC6CD', lw=.8, ls=':')

    # B: element combinations among all marker-bearing MAGs.
    fig.text(.715, .915, 'B  C/N/P/S marker combinations among\n    marker-bearing MAGs', fontsize=11.5, weight='bold', va='top')
    bars = fig.add_axes([.735, .765, .23, .11]); dots = fig.add_axes([.735, .645, .23, .095])
    xs = np.arange(len(combinations))
    bars.bar(xs, combo_table.MAG_count, color='#577B8B', width=.7)
    for i, v in enumerate(combo_table.MAG_count): bars.text(i, v + .4, str(v), ha='center', fontsize=8)
    bars.set(xlim=(-.6, len(combinations) - .4), xticks=[], ylabel='MAGs (n)', ylim=(0, combo_table.MAG_count.max() * 1.22))
    bars.grid(axis='y', c='#E5EBEF', lw=.5); bars.set_axisbelow(True)
    dots.set(xlim=(-.6, len(combinations) - .4), ylim=(3.5, -.5), xticks=[], yticks=range(4), yticklabels=list('CNPS')); clean(dots)
    for i, c in enumerate(combinations):
        dots.scatter([i] * 4, range(4), s=17, c='#DEE5EA')
        on = [j for j, v in enumerate(c) if v]
        if len(on) > 1: dots.plot([i, i], [min(on), max(on)], c='#536F7E', lw=1.4)
        for j in on: dots.scatter(i, j, s=37, c=EC['CNPS'[j]], zorder=3)
    fig.text(.715, .618, 'Each bar is one disjoint element combination; counts sum to 131.\nMarker combinations do not establish complete pathways.', fontsize=8.5, color='#52616D')
    fig.legend(handles=[Patch(fc=c, label=p) for p, c in pc.items()], loc='upper left', bbox_to_anchor=(.715, .596), ncol=2, frameon=False, fontsize=7.4, title='Phylum (GTDB)', title_fontsize=8.5, columnspacing=.8, handlelength=1.2)

    # C: focused linkage plot; percentages are unchanged source results.
    fig.text(.715, .405, 'C  Community marker-pool TPM linked to recovered MAGs', fontsize=12.5, weight='bold')
    cx = fig.add_axes([.735, .255, .23, .13])
    labels = ['DNRA terminal', 'Organic P hydrolysis', 'Denitrification steps']; colours = ['#487AA1', '#50946C', '#BE865E']
    for i, r in link.iterrows():
        cx.barh(i, r.Linked_TPM_pct, color=colours[i], height=.55)
        cx.text(r.Linked_TPM_pct + .18, i, f'{r.Linked_TPM_pct:.2f}%  |  {r.MAG_count} MAGs', va='center', fontsize=8.3, color='#344F60')
    cx.set(xlim=(0, 12.5), ylim=(2.6, -.6), yticks=range(3), yticklabels=labels, xticks=[0, 3, 6, 9, 12], xlabel='Linked community marker TPM (%)')
    cx.grid(axis='x', c='#E5EBEF', lw=.5); cx.set_axisbelow(True); cx.spines['left'].set_visible(False); cx.tick_params(axis='y', length=0)
    fig.text(.715, .192, 'Denominator = complete community marker-pool TPM.\nUnlinked signal remains unresolved at the MAG level.', fontsize=8.5, color='#52616D')
    fig.text(.025, .045, 'Marker presence indicates genomic potential. Non-detection is not confirmed biological absence; recovered MAGs represent only part of the community gene pool.', fontsize=9, color='#52616D')

    en = 'Genome-resolved architecture of nutrient-related functional potential. (A) The 131 accepted-marker-bearing MAGs are arranged by the existing GTDB tree order, with a phylum strip and an aligned binary G01-G21 marker matrix. The tree is displayed as a rectangular phylogram using supplied branch lengths, without new phylogenetic inference. (B) Each bar represents one mutually exclusive C/N/P/S marker combination among the 131 MAGs; the dot rows identify its elements. A detected combination does not demonstrate a complete metabolic pathway. (C) The focused 0-12.5% plot reports the unchanged fraction of community marker-gene TPM linked to recovered MAGs: DNRA terminal reduction 9.82% (6 MAGs), organic-P hydrolysis 3.39% (7 MAGs), and denitrification steps 2.90% (17 MAGs). These are linkage fractions, not activity, process contribution or pathway completeness. Non-detection is not confirmed biological absence.'
    cn = '基因组层面的营养相关功能架构。A按照已有GTDB系统树顺序展示131个携带接受标记的MAG，并对齐门水平色带和G01–G21二元标记矩阵；树使用已提供的分支长度进行矩形展示，未进行新的系统发育推断。B每根柱代表131个MAG中的一种互斥C/N/P/S标记组合，点阵表示该组合包含的元素；检出组合不等于完整代谢通路。C将横轴聚焦于0–12.5%，展示未改变的群落标记基因TPM中可链接到恢复MAG的比例：DNRA末端还原9.82%（6个MAG）、有机磷水解3.39%（7个MAG）和反硝化步骤2.90%（17个MAG）。这些数值是链接比例，不是活性、过程贡献率或通路完成率；未检出不等于基因真实缺失。'
    d.finish(fig, en, cn, 'Fig.5：主图恢复为矩形GTDB树与G01–G21对齐矩阵；B明确元素组合及MAG数量；C将链接覆盖率聚焦到0–12.5%，并明确其不是活性、贡献率或通路完成率。未检出不等于真实缺失。', {'shown_MAGs': 131, 'displayed_bacteria': 123, 'displayed_archaea': 8, 'linked_genes': 476, 'combination_counts_sum': int(combo_table.MAG_count.sum())}, {'element_combinations_display.tsv': combo_table, 'MAG_element_membership.tsv': e.reset_index(), 'MAG_row_order_taxonomy.tsv': order.reset_index(), 'MAG_trait_matrix_existing.tsv': binary.reset_index(), 'linkage_existing.tsv': link})


if __name__ == '__main__': main()
