"""Verify directory navigation and provenance without rerunning research."""
import ast
import hashlib
import importlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'scripts'))
from project_paths import project_path, source_matches, pre_layout_source


def verify(snapshot=False):
    manifest = json.loads((ROOT / 'maintenance/layout-migration.json').read_text())
    counts = {}
    for row in manifest['files']:
        old, new = ROOT / row['old_path'], ROOT / row['new_path']
        assert not old.exists(), f'File recreated in old location: {old}'
        assert new.is_file(), f'Missing migrated file: {new}'
        actual = hashlib.sha256(new.read_bytes()).hexdigest()
        if snapshot:
            assert actual == row['new_sha256'], f'Changed since migration: {new}'
        if snapshot and new.suffix in ('.json', '.npz', '.png', '.pdf'):
            assert actual == row['old_sha256'], f'Research artifact altered: {new}'
        assert project_path(old.name) == new
        assert project_path(old) == new
        counts[new.parent.name] = counts.get(new.parent.name, 0) + 1

    assert project_path('README.md') == ROOT / 'README.md'
    assert project_path('TASK_STATE.md') == ROOT / 'TASK_STATE.md'
    assert project_path('results/custom.json') == ROOT / 'results/custom.json'
    assert project_path(ROOT.parent/'external.json') == ROOT.parent/'external.json'
    assert not source_matches('spinful_return.py', 'invalid-hash')
    old_source = next(r for r in manifest['files'] if r['old_path']=='spinful_return.py')
    assert not source_matches('spinful_return.py', old_source['old_sha256'], 'changed-source')

    sources = sorted((ROOT / 'scripts').glob('*.py'))
    for path in sources:
        ast.parse(path.read_bytes(), filename=str(path))
        importlib.import_module(path.stem)

    pairs = [('incoming_modes', 'incoming_modes'),
             ('dipole_collar_modes', 'dipole_collar_modes'),
             ('balanced_dihole', 'balanced_dihole'),
             ('threshold_return_bound', 'threshold_return_bound'),
             ('complete_angular_source', 'complete_angular_source'),
             ('long_throat_source', 'long_throat_source'),
             ('magnetic_mouth_matching', 'magnetic_mouth_matching'),
             ('smooth_magnetic_collar', 'smooth_magnetic_collar')]
    for source, result in pairs:
        saved = json.loads(project_path(result+'_results.json').read_text())
        assert source_matches(source+'.py', saved['source_sha256']), source
    pilot = json.loads(project_path('complete_angular_source_results.json').read_text())
    for name, digest in pilot['input_sha256'].items():
        assert source_matches(name, digest), f'Frozen pilot input: {name}'

    from smooth_collar_mode import OLD_COLLAR_SHA
    old = pre_layout_source('dirac_return_closure.py')
    reverted = old.replace(b', *, return_profile=False):', b'):').replace(
        b"    record={'clock':clock", b"    return {'clock':clock").replace(
        b"    if return_profile:\n        return record, {'solution':sol,'central_radius_offset':offset,'alpha':alpha}\n    return record\n", b'')
    assert hashlib.sha256(reverted).hexdigest() == OLD_COLLAR_SHA

    from spinful_return import operator_identity
    assert operator_identity() == ['0', '0']
    broken, links = [], 0
    documents = [ROOT/n for n in ('README.md', 'DIRECTORY_GUIDE.md', 'TASK_STATE.md', 'AGENTS.md')]
    documents += list((ROOT/'docs').glob('*.md'))
    for doc in documents:
        text = re.sub(r'```.*?```', '', doc.read_text(encoding='utf-8'), flags=re.S)
        for match in re.finditer(r'\]\(([^)]+)\)', text):
            target = match[1].strip('<>').split('#', 1)[0]
            if not target or target.startswith(('http:', 'https:', 'mailto:', 'codex:')):
                continue
            if not re.search(r'\.(py|md|json|npz|png|pdf)$|/$', target):
                continue
            links += 1
            if not (doc.parent/target).exists():
                broken.append(f'{doc.relative_to(ROOT)}: {target}')
    assert not broken, '\n'.join(broken)
    assert not any(p.is_file() and p.suffix in ('.py','.json','.npz','.png','.pdf') for p in ROOT.iterdir())
    print(json.dumps({'migrated_files': counts, 'modules_imported': len(sources),
                      'historical_source_guards': len(pairs),
                      'frozen_pilot_hashes': len(pilot['input_sha256']),
                      'markdown_links': links, 'operator_identity': 'pass'}, indent=2))


if __name__ == '__main__':
    verify(snapshot='--snapshot' in sys.argv)
