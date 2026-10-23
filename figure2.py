"""Figure 2: cross-level organisation of microbial communities and functions.

This script is a display-only replot. It reads saved ordination scores,
PERMANOVA rows, cross-layer correspondence statistics and db-RDA estimates;
no statistic is recomputed here.
"""

from _support import *
import os
import tempfile
from shutil import copy2

PAIR_ORDER = [
    ("Taxonomy", "KEGG", "Taxonomy–KEGG"),
    ("Taxonomy", "Nutrient", "Taxonomy–Nutrient markers"),
    ("KEGG", "Nutrient", "KEGG–Nutrient markers"),
]


def _column(frame, candidates, label):
    for name in candidates:
        if name in frame.columns:
            return name
    raise KeyError(f"{label}: none of {candidates!r} is present; columns={list(frame.columns)!r}")


def _profile_name(value):
    text = str(value).strip().lower().replace("_", " ").replace("-", " ")
    if text.startswith("tax"):
        return "Taxonomy"
    if text.startswith("kegg"):
        return "KEGG"
    if "nutrient" in text or "marker" in text:
        return "Nutrient"
    return str(value)


def _stat(value):
    return f"{float(value):.4f}"


def _style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.spines["left"].set_linewidth(0.65)
    ax.spines["bottom"].set_linewidth(0.65)
    ax.tick_params(axis="both", length=3, width=0.55, color="#222222", pad=2)


def _read_conditional_p(env_row):
    """Read stored conditional P when present, otherwise use the source note."""
    candidates = [
        "Partial_period_P", "Month_conditioned_P", "Month_conditioned_p",
        "P_period", "Partial_period_p", "Partial_period_P_value",
    ]
    for name in candidates:
        if name in env_row.index and pd.notna(env_row[name]):
            value = float(env_row[name])
            if not np.isclose(value, 0.0694, atol=5e-7):
                raise ValueError(
                    "Existing KEGG month-conditioned P differs from documented "
                    f"0.0694: {name}={value:g}. Resolve the source discrepancy first."
                )
            return value, name
    return 0.0694, "original figure note (no P column in curated table)"


def main():
    d = Delivery(2, __file__)
    scores = d.read(P + "layer_PCoA_scores.tsv", "A–C", "Existing 36-layer Bray–Curtis PCoA scores")
    per = d.read(P + "Fig3B_PERMANOVA.tsv", "A–C", "Existing partial R² and BH q values; no new test")
    corr = d.read(P + "Fig3C_cross_layer_associations.tsv", "D", "Existing 12-core Mantel and PROTEST results")
    env = d.read(V + "Fig3D_environment_source.tsv", "E", "Existing bulk-chemistry db-RDA adjusted R²")
    d.evidence("03_workflows/R/46_figures2_6_curated_v13.R", "Methods", "Original statistics and permutation design; not executed")

    assert scores.groupby("Profile").size().eq(36).all(), "Each PCoA profile must have 36 layers."
    expected_pct = {"Taxonomy": (23.80, 10.50), "KEGG": (28.99, 18.69), "Nutrient": (33.02, 12.27)}
    scores["_display_profile"] = scores["Profile"].map(_profile_name)
    if set(expected_pct) - set(scores["_display_profile"]):
        raise ValueError(f"PCoA profile labels do not contain the expected layers: {sorted(scores['_display_profile'].unique())!r}")
    axis_pct_differences = {}
    for profile, (axis1, axis2) in expected_pct.items():
        z = scores[scores["_display_profile"].eq(profile)]
        # The documented percentages are printed to two decimals. Compare at
        # that reporting precision so harmless floating-point/rounding noise
        # does not abort plotting, while still rejecting a real source mismatch.
        axis1_diff = float(np.max(np.abs(z["Axis1_pct"].to_numpy(dtype=float) - axis1)))
        axis2_diff = float(np.max(np.abs(z["Axis2_pct"].to_numpy(dtype=float) - axis2)))
        axis_pct_differences[profile] = {"PCoA1_max_abs_diff": axis1_diff, "PCoA2_max_abs_diff": axis2_diff}
        if axis1_diff > 0.01 or axis2_diff > 0.01:
            raise ValueError(f"{profile} axis percentages differ from documented values: {z[['Axis1_pct', 'Axis2_pct']].drop_duplicates().to_dict('records')}")

    for name in ["Mantel_r", "Procrustes_r", "Mantel_q", "Procrustes_q"]:
        if name not in corr.columns:
            raise KeyError(f"D correspondence table lacks required column {name!r}.")
    rows = []
    for layer1, layer2, label in PAIR_ORDER:
        hit = corr[corr.Layer1.eq(layer1) & corr.Layer2.eq(layer2)]
        if len(hit) != 1:
            raise ValueError(f"Expected one D row for {layer1}–{layer2}, found {len(hit)}.")
        row = hit.iloc[0]
        if min(float(row.Mantel_r), float(row.Procrustes_r)) < 0:
            raise ValueError("D contains a negative coefficient although the documented range is 0–1.")
        rows.append(row)

    env = env.copy()
    env["_display_profile"] = env["Profile"].map(_profile_name)
    env_order = ["Taxonomy", "KEGG", "Nutrient"]
    if set(env_order) - set(env["_display_profile"]):
        raise ValueError(f"E table lacks one of the three response layers: {list(env.Profile)!r}")
    env = env.set_index("_display_profile").loc[env_order].reset_index()
    overall_col = _column(env, ["Adjusted_R2", "Overall_adjusted_R2"], "E overall adjusted R²")
    conditional_col = _column(env, ["Partial_period_adjusted_R2", "Month_conditioned_adjusted_R2"], "E month-conditioned adjusted R²")
    q_col = _column(env, ["q_BH", "Overall_q_BH"], "E overall BH q")
    kegg_row = env.loc[env["_display_profile"].eq("KEGG")].iloc[0]
    overall_kegg_q = float(kegg_row[q_col])
    if not np.isclose(overall_kegg_q, 0.0159, atol=5e-7):
        raise ValueError(f"Existing KEGG overall q differs from documented 0.0159: {q_col}={overall_kegg_q:g}.")
    conditional_kegg_p, conditional_p_source = _read_conditional_p(kegg_row)

    plt.rcParams.update({"font.family": "Arial", "font.size": 8.0, "axes.titlesize": 10.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(14.2, 10.5), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=0.060, right=0.985, top=0.865, bottom=0.205, height_ratios=[1.00, 1.16], hspace=0.40, wspace=0.28)

    panel_titles = {"Taxonomy": "A  Taxonomic PCoA", "KEGG": "B  KEGG functional PCoA", "Nutrient": "C  Nutrient-marker PCoA"}
    supported_terms = []
    for i, profile in enumerate(["Taxonomy", "KEGG", "Nutrient"]):
        ax = fig.add_subplot(gs[0, i])
        z = scores[scores["_display_profile"].eq(profile)]
        for season in SEASONS:
            for depth, marker in [("Surface", "o"), ("Middle", "^"), ("Deep", "s")]:
                q = z[z.OriginalSeasonLabel.eq(season) & z.Depth.eq(depth)]
                if not q.empty:
                    ax.scatter(q.PCoA1, q.PCoA2, s=29, c=SC[season], marker=marker, edgecolors="white", linewidths=0.42, alpha=0.88, zorder=3)
        ax.axhline(0, color="#D9E0E4", linewidth=0.60, zorder=1)
        ax.axvline(0, color="#D9E0E4", linewidth=0.60, zorder=1)
        ax.set(xlabel=f"PCoA1 ({float(z.Axis1_pct.iloc[0]):.2f}%)", ylabel=f"PCoA2 ({float(z.Axis2_pct.iloc[0]):.2f}%)")
        ax.set_title(panel_titles[profile], loc="left", weight="bold", pad=8)
        ax.set_aspect("equal", adjustable="datalim")
        ax.margins(0.12)
        _style_axes(ax)
        if {"Profile", "Term", "q_BH", "partial_R2"}.issubset(per.columns):
            source_profile = z.Profile.iloc[0]
            supported = per[per.Profile.eq(source_profile) & per.Term.isin(["SamplingPeriod", "Depth"]) & per.q_BH.lt(0.05)]
            if len(supported):
                lines = [("Depth" if r.Term == "Depth" else "Sampling month") + f" partial R² = {float(r.partial_R2):.3f}" for r in supported.itertuples()]
                ax.text(0.03, 0.97, "PERMANOVA\n" + "\n".join(lines), transform=ax.transAxes, va="top", fontsize=7.0, color="#465B67")
                supported_terms.extend(f"{profile}: {line}" for line in lines)

    month_handles = season_handles()
    depth_handles = [Line2D([], [], linestyle="", marker=m, color="#5E707C", label=label, markersize=5.4) for m, label in [("o", "0–30 cm"), ("^", "30–60 cm"), ("s", "60–90 cm")]]
    fig.legend(handles=month_handles + depth_handles, loc="upper center", bbox_to_anchor=(0.52, 0.965), ncol=7, frameon=False, handletextpad=0.35, columnspacing=1.05)

    bottom = gs[1, :].subgridspec(1, 2, width_ratios=[1.35, 1.00], wspace=0.28)
    douter = bottom[0].subgridspec(2, 1, height_ratios=[3.15, 1.12], hspace=0.25)
    dtop = douter[0].subgridspec(1, 2, width_ratios=[3.15, 1.35], wspace=0.025)
    dcoef = fig.add_subplot(dtop[0])
    dqax = fig.add_subplot(dtop[1], sharey=dcoef)
    dtable = fig.add_subplot(douter[1])
    dcoef.set_title("D  Cross-level correspondence", loc="left", weight="bold", pad=9)
    y = np.arange(len(rows))[::-1]
    labels = [p[2] for p in PAIR_ORDER]
    # The third comparison is the strongest in both statistics; a restrained
    # row tint and bold label make that pattern discoverable without changing
    # point sizes or encoding the conclusion by colour alone.
    dcoef.axhspan(-0.45, 0.45, color="#E8F2F1", alpha=0.72, zorder=0)
    for idx, (yi, row) in enumerate(zip(y, rows)):
        mantel, protest = float(row.Mantel_r), float(row.Procrustes_r)
        edge = "#28786D" if idx == 2 else "#2F7189"
        dcoef.plot([mantel, protest], [yi, yi], color="#AABEC5", linewidth=1.10, zorder=1)
        dcoef.scatter(mantel, yi, s=44, color=edge, edgecolor="white", linewidth=0.45, zorder=3)
        dcoef.scatter(protest, yi, s=44, facecolor="white", edgecolor=edge, linewidth=1.25, zorder=3)
        # Values are printed at the points, with opposing offsets so close
        # Mantel/PROTEST estimates remain legible at journal size.
        dcoef.annotate(f"{mantel:.3f}", (mantel, yi), xytext=(-2, 7), textcoords="offset points", ha="right", va="bottom", fontsize=6.8, color="#405763")
        dcoef.annotate(f"{protest:.3f}", (protest, yi), xytext=(3, -8), textcoords="offset points", ha="left", va="top", fontsize=6.8, color="#405763")
    dcoef.set(xlim=(0, 1), ylim=(-0.55, len(rows) - 0.45), yticks=y, yticklabels=labels, xlabel="Correspondence coefficient (r)", xticks=np.linspace(0, 1, 5))
    for tick, label in zip(dcoef.get_yticklabels(), labels):
        if label == "KEGG–Nutrient markers":
            tick.set_weight("bold")
            tick.set_color("#28786D")
    dcoef.grid(axis="x", color="#E2E8EB", linewidth=0.60)
    dcoef.set_axisbelow(True)
    _style_axes(dcoef)
    dcoef.legend(handles=[Line2D([], [], linestyle="", marker="o", color="#2F7189", label="Mantel r", markersize=5.5), Line2D([], [], linestyle="", marker="o", markerfacecolor="white", markeredgecolor="#2F7189", label="PROTEST r", markersize=5.5)], loc="lower right", frameon=False, ncol=2, handletextpad=0.25, columnspacing=0.75, borderaxespad=0.1)

    dqax.set(xlim=(0, 1), ylim=dcoef.get_ylim(), xticks=[], yticks=y)
    dqax.tick_params(axis="y", left=False, labelleft=False)
    dqax.set_frame_on(False)
    dqax.text(0.00, 1.04, "Mantel q", transform=dqax.transAxes, ha="left", va="bottom", weight="bold", color="#405763")
    dqax.text(0.57, 1.04, "PROTEST q", transform=dqax.transAxes, ha="left", va="bottom", weight="bold", color="#405763")
    for yi, row in zip(y, rows):
        dqax.text(0.00, yi, _stat(row.Mantel_q), transform=dqax.get_yaxis_transform(), ha="left", va="center", color="#405763")
        dqax.text(0.57, yi, _stat(row.Procrustes_q), transform=dqax.get_yaxis_transform(), ha="left", va="center", color="#405763")

    dtable.axis("off")
    dtable.text(0.00, 1.04, "Within-month permutation P", transform=dtable.transAxes, ha="left", va="bottom", weight="bold", color="#405763", fontsize=7.6)
    table_data = [[p[2], _stat(r.Within_period_Mantel_P), _stat(r.Within_period_Procrustes_P)] for p, r in zip(PAIR_ORDER, rows)]
    tab = dtable.table(cellText=table_data, colLabels=["Comparison", "Mantel P", "PROTEST P"], colWidths=[0.53, 0.22, 0.25], cellLoc="left", colLoc="left", bbox=[0.0, 0.0, 1.0, 0.88])
    tab.auto_set_font_size(False)
    tab.set_fontsize(7.0)
    for (ri, ci), cell in tab.get_celld().items():
        cell.set_edgecolor("#D9E2E6")
        cell.set_linewidth(0.45)
        cell.set_facecolor("#F7F9FA" if ri == 0 else "white")
        if ri == 0:
            cell.set_text_props(weight="bold", color="#405763")

    eax = fig.add_subplot(bottom[1])
    eax.set_title("E  Explanation by measured bulk chemistry", loc="left", weight="bold", pad=9)
    ey = np.arange(len(env_order))[::-1]
    x_values = []
    for yi, profile in zip(ey, env_order):
        row = env[env["_display_profile"].eq(profile)].iloc[0]
        overall, conditional = float(row[overall_col]), float(row[conditional_col])
        x_values.extend([overall, conditional])
        eax.plot([overall, conditional], [yi + 0.115, yi - 0.115], color="#B8C8CF", linewidth=0.85, zorder=1)
        eax.scatter(overall, yi + 0.115, s=44, color="#356F8C", edgecolor="white", linewidth=0.45, zorder=3)
        eax.scatter(conditional, yi - 0.115, s=44, facecolor="white", edgecolor="#356F8C", linewidth=1.2, zorder=3)
        eax.annotate(f"{overall:.3f}", (overall, yi + 0.115), xytext=(4, 5), textcoords="offset points", fontsize=6.7, color="#405763", va="bottom")
        eax.annotate(f"{conditional:.3f}", (conditional, yi - 0.115), xytext=(4, -6), textcoords="offset points", fontsize=6.7, color="#405763", va="top")
        q_value = float(row[q_col])
        if q_value < 0.05 or profile == "KEGG":
            eax.annotate(f"q={q_value:.4f}", (overall, yi + 0.115), xytext=(5, 4), textcoords="offset points", fontsize=6.8, color="#405763", va="bottom")
        if profile == "KEGG":
            eax.annotate(f"P={conditional_kegg_p:.4f}", (conditional, yi - 0.115), xytext=(5, -5), textcoords="offset points", fontsize=6.8, color="#405763", va="top")
    xmax = max(0.40, float(np.nanmax(x_values)) * 1.22)
    xmax = min(1.0, np.ceil(xmax * 20) / 20)
    eax.set(xlim=(0, xmax), ylim=(-0.55, len(env_order) - 0.45), yticks=ey, yticklabels=["Taxonomy", "KEGG", "Nutrient markers"], xlabel="Adjusted R²", xticks=np.linspace(0, xmax, 5))
    eax.grid(axis="x", color="#E2E8EB", linewidth=0.60)
    eax.set_axisbelow(True)
    eax.tick_params(axis="y", length=0)
    eax.spines["left"].set_visible(False)
    _style_axes(eax)
    eax.legend(handles=[Line2D([], [], linestyle="", marker="o", color="#356F8C", label="Overall", markersize=5.5), Line2D([], [], linestyle="", marker="o", markerfacecolor="white", markeredgecolor="#356F8C", label="Month-conditioned", markersize=5.5)], loc="lower right", frameon=False, ncol=1, handletextpad=0.3, borderaxespad=0.1)

    fig.text(0.060, 0.136, "A–C: 36 stratified samples nested in 12 sediment cores.", fontsize=7.4, color="#52616D", va="bottom")
    fig.text(0.060, 0.104, "D: Mantel and Procrustes correspondence across three layer pairs; BH q values and month-conditioned permutation P values are reported separately for 12 equal-weight core profiles.", fontsize=7.4, color="#52616D", va="bottom")
    fig.text(0.060, 0.072, "E: OM, TN, TP and pH explain taxonomy, KEGG and nutrient-marker profiles; water level is not included. Overall KEGG q = 0.0159; month-conditioned P = 0.0694.", fontsize=7.4, color="#52616D", va="bottom")

    legend_en = ("Cross-level organization of microbial communities and functions. (A–C) Existing Bray–Curtis PCoA coordinates for 36 stratified samples nested within 12 sediment cores. Colours indicate February, May, August and November; point shapes indicate 0–30, 30–60 and 60–90 cm. Axis percentages and coordinates are read unchanged from the curated score table. The point clouds are descriptive; no new ellipse, cluster boundary or significance claim is added. (D) Existing correspondence statistics for Taxonomy–KEGG, Taxonomy–Nutrient markers and KEGG–Nutrient markers using 12 equal-weight core profiles. Filled circles are Mantel r and open circles are two-dimensional PROTEST r. BH-adjusted q values are in a separate aligned column, and within-month permutation Mantel and PROTEST P values are in the lower table. (E) Existing adjusted R² for models using OM, TN, TP and pH, shown overall and conditioned on sampling month. The connector is a visual aid; the two endpoints remain separate model estimates. Overall KEGG q = 0.0159; month-conditioned KEGG P = 0.0694. Water level is not included in the current model.")
    legend_cn = ("跨层级微生物群落与功能组织。A–C直接使用36个分层样品的Bray–Curtis PCoA坐标，样品嵌套于12根沉积柱；颜色为February、May、August和November，形状为0–30、30–60和60–90 cm。坐标解释率与点位来自原始结果表，未添加椭圆、聚类边界或新的显著性推断。D使用12根柱的等权核心剖面展示Taxonomy–KEGG、Taxonomy–Nutrient markers和KEGG–Nutrient markers对应关系，实心点为Mantel r，空心点为二维PROTEST r；BH校正q值独立列出，采样月份内限制置换的Mantel与PROTEST P值置于下方表格。E展示使用OM、TN、TP和pH的整体及控制采样月份后的调整R²；连接线只用于视觉比较，两个端点仍是分别的模型估计。总体KEGG q=0.0159，控制月份后KEGG P=0.0694，水位未纳入当前模型。")
    changes = ("Fig.2 revised/v2：A–C统一为等尺寸PCoA面板，顶部只保留一份月份与深度图例；移除过大的整图标题。D改为0–1横向效应量点图，Mantel r使用实心圆、PROTEST r使用空心圆，Mantel/PROTEST BH q值置于右侧独立对齐列，within-month置换P值置于下方小表。E用连接线展示整体与月份条件化调整R²，取消漂浮星号，并在对应KEGG点旁保留q=0.0159与P=0.0694。全部数值直接读取已有结果表；未运行新检验。")
    checks = {"PCoA_points_per_profile": {k: int((scores["_display_profile"] == k).sum()) for k in expected_pct}, "PCoA_axis_percentages_verified": True, "PCoA_axis_percentage_max_abs_differences": axis_pct_differences, "D_pairs": len(rows), "D_q_values_from_source": True, "D_within_month_P_values_from_source": True, "E_overall_KEGG_q": overall_kegg_q, "E_month_conditioned_KEGG_P": conditional_kegg_p, "E_conditional_P_source": conditional_p_source, "supported_PERMANOVA_annotations": supported_terms, "new_tests": False}

    fd, tmp_svg_name = tempfile.mkstemp(prefix="figure2_", suffix=".svg")
    os.close(fd)
    tmp_svg = Path(tmp_svg_name)
    try:
        fig.savefig(tmp_svg, format="svg", facecolor="white", bbox_inches="tight")
        d.finish(fig, legend_en, legend_cn, changes, checks, {"correspondence_existing_statistics.tsv": corr, "environment_existing_statistics.tsv": env.drop(columns=["_display_profile"]), "PERMANOVA_existing_statistics.tsv": per, "within_month_permutation_statistics.tsv": pd.DataFrame([{"Comparison": label, "Mantel_P_within_month": float(row.Within_period_Mantel_P), "PROTEST_P_within_month": float(row.Within_period_Procrustes_P)} for (_, _, label), row in zip(PAIR_ORDER, rows)])})
        os.replace(tmp_svg, d.out / "Figure_2.svg")
        # Keep explicit revised filenames alongside the standard delivery names.
        for ext in ("pdf", "png", "svg"):
            copy2(d.out / f"Figure_2.{ext}", d.out / f"Fig2_revised.{ext}")
    finally:
        if tmp_svg.exists():
            tmp_svg.unlink()


if __name__ == "__main__":
    main()
