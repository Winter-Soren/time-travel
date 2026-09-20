"""Paths and explicit source-provenance translation for the directory migration."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'maintenance/layout-migration.json'
FOLDERS = {'.py': 'scripts', '.json': 'results', '.npz': 'data',
           '.png': 'figures', '.pdf': 'figures', '.md': 'docs'}
ROOT_DOCUMENTS = {'README.md', 'DIRECTORY_GUIDE.md', 'TASK_STATE.md', 'AGENTS.md'}


def project_path(name):
    """Resolve project filenames; preserve explicit external/custom paths."""
    path = Path(name)
    if path.is_absolute():
        if path.parent != ROOT:
            return path
        path = Path(path.name)  # A legacy absolute filename at the project root.
    if path.parent != Path('.'):
        return ROOT / path
    folder = '' if path.name in ROOT_DOCUMENTS else FOLDERS.get(path.suffix, '')
    return ROOT / folder / path


def source_matches(name, recorded_sha256, actual_sha256=None):
    """Accept exact bytes or an explicitly recorded path-only source migration.

    New results still record their actual source hashes. Historical result hashes
    are never rewritten. Any later source edit invalidates this translation.
    """
    path = project_path(name)
    actual = actual_sha256 or hashlib.sha256(path.read_bytes()).hexdigest()
    if actual == recorded_sha256:
        return True
    if path.suffix != '.py':
        return False
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != manifest['path_helper_sha256']:
        return False
    relative = path.relative_to(ROOT).as_posix()
    return any(row['new_path'] == relative and row['old_sha256'] == recorded_sha256
               and row['new_sha256'] == actual for row in manifest['files'])


def pre_layout_source(name):
    """Read the authenticated historical bytes needed by the collar audit."""
    path = ROOT / 'provenance/pre-layout' / Path(name).name
    data = path.read_bytes()
    if not source_matches(name, hashlib.sha256(data).hexdigest()):
        raise ValueError(f'Historical source does not match the migration: {name}')
    return data
