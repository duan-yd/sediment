"""Package already-rendered figures and verify source hashes. Does not generate any figure."""
from pathlib import Path
import hashlib,json,zipfile,argparse,re
import pandas as pd
from PIL import Image

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'05_figures/science_questions_20261013_final'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replace',action='store_true',help='Replace delivery indexes after backing up the previous edition.')
    args=parser.parse_args()
    maps=[];legends=[];changes=[];checks=[]
    for n in range(2,6):
        f=OUT/f'Figure{n}'
        record=json.loads((f/'PROVENANCE_QA.json').read_text(encoding='utf-8'))
        assert record['sources_unchanged'] and not record['new_hypothesis_tests']
        assert all(sha(Path(p))==v for p,v in record['source_SHA256'].items())
        maps.append(pd.read_csv(f/'INPUT_FILES.tsv',sep='\t'))
        legends.append((f/'LEGEND_EN_CN.md').read_text(encoding='utf-8'))
        changes.append(f'## Figure {n}\n\n'+(f/'PANEL_CHANGES_CN.md').read_text(encoding='utf-8'))
        with Image.open(f/f'Figure_{n}.png') as im:
            assert min(im.size)>3000 and abs(im.info['dpi'][0]-400)<1
            size=im.size
        assert (f/f'Figure_{n}.pdf').read_bytes().startswith(b'%PDF')
        checks.append({'Figure':n,'Sources_unchanged':True,'PNG_pixels':str(size),'PNG_dpi':400,'PDF_header_valid':True})
    for name in ['ALL_INPUT_FILES.tsv','ALL_LEGENDS_EN_CN.md','ALL_PANEL_CHANGES_CN.md','DELIVERY_QA.tsv','FIGURES2_6_SCRIPTS_LEGENDS.zip']:
        if (OUT/name).exists() and not args.replace:raise FileExistsError(f'Refuse overwrite: {OUT/name}')
    pd.concat(maps,ignore_index=True).to_csv(OUT/'ALL_INPUT_FILES.tsv',sep='\t',index=False,encoding='utf-8-sig')
    (OUT/'ALL_LEGENDS_EN_CN.md').write_text('\n\n---\n\n'.join(legends),encoding='utf-8')
    (OUT/'ALL_PANEL_CHANGES_CN.md').write_text('\n\n'.join(changes),encoding='utf-8')
    pd.DataFrame(checks).to_csv(OUT/'DELIVERY_QA.tsv',sep='\t',index=False)
    # Legacy archive filename is retained; content is current Figures 2–5 only.
    paths=[OUT/name for name in ['README_CN.md','ALL_INPUT_FILES.tsv','ALL_LEGENDS_EN_CN.md','ALL_PANEL_CHANGES_CN.md','DELIVERY_QA.tsv']]
    table_names={2:['correspondence_existing_statistics.tsv','environment_existing_statistics.tsv','PERMANOVA_existing_statistics.tsv'],3:['marker_dictionary.tsv','vertical_contrast_display.tsv','vertical_core_source_values.tsv','existing_depth_profile_means.tsv'],4:['network_display_edges.tsv','mean_resolution_coverage.tsv','existing_core_diversity.tsv','existing_contributor_period_effects.tsv'],5:['element_combinations_display.tsv','MAG_element_membership.tsv','MAG_row_order_taxonomy.tsv','MAG_trait_matrix_existing.tsv','linkage_existing.tsv']}
    for n in range(2,6):
        f=OUT/f'Figure{n}'
        paths.extend(f/name for name in [f'Figure_{n}.png',f'Figure_{n}.pdf',f'Figure_{n}_preview.png','INPUT_FILES.tsv','LEGEND_EN_CN.md','PANEL_CHANGES_CN.md','PROVENANCE_QA.json'])
        for name in table_names[n]:
            assert (f/name).exists(),f'Missing current display table: {f/name}'
            paths.append(f/name)
    scripts=[Path(__file__).parent/name for name in ['_support.py','package_delivery.py','figure2.py','figure3.py','figure4.py','figure5.py']]
    with zipfile.ZipFile(OUT/'FIGURES2_6_SCRIPTS_LEGENDS.zip','w' if args.replace else 'x',zipfile.ZIP_DEFLATED) as z:
        for p in paths:z.write(p,p.relative_to(OUT))
        for p in scripts:z.write(p,'scripts/'+p.name)
    print(json.dumps({'figures':checks,'zip':str(OUT/'FIGURES2_6_SCRIPTS_LEGENDS.zip'),'direct_input_hashes_verified':True}))
if __name__=='__main__':main()
