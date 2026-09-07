"""Read-only numerical/source consistency checks plus a derived audit JSON."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
import numpy as np
from PIL import Image
import xml.etree.ElementTree as ET
from plot_fault_event_alignment import read_csv


def main():
    p=argparse.ArgumentParser();p.add_argument('result',type=Path)
    root=p.parse_args().result.resolve()
    manifest=json.loads((root/'run_manifest.json').read_text(encoding='utf-8'))
    assert manifest['complete']
    for name, expected in manifest['code_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    rows=[]
    for condition in ['BL','PF50']:
        with np.load(root/condition/'evaluated_envelopes.npz') as data:
            assert np.all(np.diff(data['t'])>0)
            assert np.max(np.abs(data['path'].sum(axis=1)))<1e-12
            test=np.zeros(len(data['t']),bool)
            for a,b in manifest['split_seconds']['test']:
                test|=(data['t']>=a)&(data['t']<b)
            for method in ['none','waveform_xcorr','v1_constraint','envelope_xcorr','joint_envelope']:
                assert np.isfinite(data[method]).all()
                assert np.max(np.abs(data[method+'_shift']))<=.0080001
            assignments=data['joint_assignment'][test]
            valid=assignments>=0
            rows.append(dict(condition=condition,main_comparison_events=int(np.sum(test&data['main_mask'])),
                five_method_comparison_events=int(np.sum(test&data['common_mask'])),
                switched_supported_channel_fraction=float(np.mean(assignments[valid]>0))))
    metrics=read_csv(root/'metrics.csv')
    assert all(float(r['time_map_min_derivative'])>=.199 for r in metrics)
    exports=[]
    for svg in sorted((root/'figures').glob('*.svg')):
        count=len(ET.parse(svg).findall('.//{http://www.w3.org/2000/svg}text'))
        assert count>10
        with Image.open(svg.with_suffix('.tiff')) as raster:
            assert min(raster.info['dpi'])>=599
            size=raster.size
        for ext in ['png','pdf','tiff']:
            assert svg.with_suffix('.'+ext).stat().st_size>1000
        exports.append(dict(figure=svg.stem,editable_svg_text_nodes=count,tiff_pixels=size,tiff_dpi=600))
    tests=json.loads((ROOT/'results'/'event_arrival_v2_tests.json').read_text(encoding='utf-8'))
    assert tests['passed']
    qa=json.loads((root/'figure_qa.json').read_text(encoding='utf-8'))
    assert all(not item['text_outside_canvas'] for item in qa['figures'])
    audit=dict(passed=True,core_code_matches_executed_manifest=True,time_maps_monotonic=True,
               local_shifts_bounded_8ms=True,relative_path_zero_mean=True,event_arrays_finite=True,
               support_and_candidate_switches=rows,exports=exports,controlled_tests=tests,
               limitation='Numerical and export checks do not establish physical fault identity or unique speed/path separation.')
    (root/'numerical_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
