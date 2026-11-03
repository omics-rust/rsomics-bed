"""Verify source packaging and the installed fourteen-operation BED candidate."""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile

PRODUCT = 'rsomics-bed'
ROOT = Path(__file__).resolve().parents[1]
OPERATIONS = ['sort', 'merge', 'intersect', 'subtract', 'complement', 'cluster', 'window', 'closest', 'coverage', 'annotate', 'map', 'genomecov', 'slop', 'shift']
HOSTS = {'linux-x86_64': 'x86_64-unknown-linux-gnu', 'linux-aarch64': 'aarch64-unknown-linux-gnu',
         'macos-x86_64': 'x86_64-apple-darwin', 'macos-aarch64': 'aarch64-apple-darwin'}


def identity(path):
    body = Path(path).read_bytes()
    return dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def package_contract(package):
    parent = Path(package['manifest_path']).parent
    result = {k: v for k, v in package.items() if k not in ('id', 'manifest_path', 'targets')}
    result['targets'] = [{**t, 'src_path': str(Path(t['src_path']).relative_to(parent))} for t in package['targets']]
    return result


def extract_package(source, destination, expected, head):
    destination = Path(destination)
    destination.mkdir(parents=True)
    with tarfile.open(source, 'r:gz') as archive:
        members = archive.getmembers()
        require(len(members) == len({m.name for m in members}), 'duplicate source package member')
        require(len(members) == len({m.name.casefold() for m in members}), 'case-colliding source package member')
        bodies = {}
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and '..' not in path.parts and path.parts[0] == destination.name,
                    'unsafe source package path')
            require(member.isfile() or member.isdir(), 'unsupported source package member')
            if member.isfile():
                name = '/'.join(path.parts[1:])
                mode = expected[name]['mode'] if name in expected else 0o644
                require(name and member.mode == mode, 'source package file mode differs: ' + name)
                bodies[name] = archive.extractfile(member).read()
        require(set(bodies) == set(expected) | {'Cargo.toml.orig', '.cargo_vcs_info.json'},
                'complete source package membership differs')
        require(bodies['Cargo.toml.orig'] == expected['Cargo.toml']['body'], 'original Cargo manifest differs')
        vcs = json.loads(bodies['.cargo_vcs_info.json'])
        require(vcs['git']['sha1'] == head and not vcs['git'].get('dirty', False) and vcs['path_in_vcs'] == '',
                'source package is not the exact clean CI head')
        for name, row in expected.items():
            if name != 'Cargo.toml':
                require(bodies[name] == row['body'], 'source package body differs: ' + name)
        for name, body in bodies.items():
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('xb') as stream:
                stream.write(body)
            path.chmod(expected[name]['mode'] if name in expected else 0o644)
    return dict(**identity(source), members=len(bodies), original_members=len(expected),
                original_bodies_modes_match=True, source_commit=head, clean_vcs=True)


def source():
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    metadata = json.loads((delivery / 'package-metadata.json').read_bytes())
    package, = [p for p in metadata['packages'] if p['name'] == PRODUCT]
    prefix = PRODUCT + '-' + package['version']
    archive = Path(metadata['target_directory']) / 'package' / (prefix + '.crate')
    shutil.copy2(archive, delivery / archive.name)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    require(head == os.environ['GITHUB_SHA'], 'source head differs from actual CI head')
    require(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT) == b'', 'source worktree is dirty')
    expected = {}
    for row in subprocess.check_output(['git', 'ls-tree', '-rz', head], cwd=ROOT).rstrip(b'\0').split(b'\0'):
        header, name = row.split(b'\t')
        mode, kind, blob = header.decode().split()
        name = name.decode()
        require(kind == 'blob' and mode in ('100644', '100755'), 'unexpected source Git type/mode: ' + name)
        body = (ROOT / name).read_bytes()
        require(hashlib.sha1(b'blob ' + str(len(body)).encode() + b'\0' + body).hexdigest() == blob,
                'source differs from exact Git body: ' + name)
        expected[name] = dict(mode=int(mode, 8) & 0o777, body=body)
    destination = temporary / 'packaged-source' / prefix
    verification = extract_package(archive, destination, expected, head)
    unpacked = json.loads(subprocess.check_output(['cargo', 'metadata', '--manifest-path', str(destination / 'Cargo.toml'),
                                                 '--locked', '--offline', '--no-deps', '--format-version', '1']))
    installed, = [p for p in unpacked['packages'] if p['name'] == PRODUCT]
    require(package_contract(installed) == package_contract(package), 'normalized source package metadata differs')
    verification['normalized_metadata_match'] = True
    save(delivery / 'package-verification.json', verification)
    metadata_body = subprocess.check_output(['cargo', 'metadata', '--locked', '--format-version', '1'], cwd=ROOT)
    tree_body = subprocess.check_output(['cargo', 'tree', '--locked', '--offline', '--target', os.environ['EXPECTED_HOST'],
                                        '--edges', 'normal,build', '--prefix', 'none', '--format', '{p}'], cwd=ROOT)
    (delivery / 'native-dependencies.metadata.json').write_bytes(metadata_body)
    (delivery / 'native-dependencies.tree.txt').write_bytes(tree_body)
    by_name = {(p['name'], p['version']): p for p in json.loads(metadata_body)['packages']}
    enabled = set()
    for line in tree_body.decode().splitlines():
        name, version = re.match(r'([^ ]+) v([^ ]+)', line).groups()
        if name != PRODUCT:
            enabled.add(by_name[name, version]['id'])
    inventory = json.loads((ROOT / 'licenses/provenance.json').read_bytes())
    require(enabled == set(inventory['target_graphs'][os.environ['EXPECTED_HOST']]), 'enabled target dependency inventory differs')
    require(identity(ROOT / 'Cargo.lock')['sha256'] == inventory['lock_sha256'], 'notice inventory lock differs')
    for row in [*inventory['copied_texts'].values(), inventory['upstream_source']]:
        require(identity(ROOT / row['destination']) == {k: row[k] for k in ('bytes', 'sha256')}, 'bundled notice body differs')
    save(delivery / 'native-dependencies.json', dict(target=os.environ['EXPECTED_HOST'], enabled_packages=sorted(enabled),
         inventory_sha256=identity(ROOT / 'licenses/provenance.json')['sha256'],
         metadata_sha256=hashlib.sha256(metadata_body).hexdigest(), tree_sha256=hashlib.sha256(tree_body).hexdigest(),
         enabled_normal_build_graph_verified=True, linked_binary_inventory=False))
    with Path(os.environ['GITHUB_ENV']).open('a') as stream:
        stream.write('PACKAGE_DIRECTORY=' + str(destination) + '\n')
        stream.write('PACKAGE_VERSION=' + package['version'] + '\n')
        stream.write('SOURCE_PACKAGE=' + str(delivery / archive.name) + '\n')


def smoke(binary):
    binary = Path(binary).resolve(strict=True)
    directory = Path(os.environ['RUNNER_TEMP']) / 'delivery/lifecycle'
    directory.mkdir(parents=True)
    shutil.copy2(binary, directory / 'installed-binary')
    fixture = ROOT / 'tests/golden'
    records = []

    def run(arguments, *, stdin=None, status=0):
        command = [str(binary), *map(str, arguments)]
        index = len(records)
        row = dict(argv=command, stdin=None if stdin is None else dict(path=str(stdin), **identity(stdin)))
        records.append(row)
        save(directory / 'commands.json', records)
        with (directory / f'{index}.stdout').open('xb') as out, (directory / f'{index}.stderr').open('xb') as err:
            result = subprocess.run(command, cwd=directory, input=None if stdin is None else Path(stdin).read_bytes(),
                                    stdout=out, stderr=err, env=dict(os.environ, NO_COLOR='1'), timeout=120)
        row['returncode'] = result.returncode
        save(directory / 'commands.json', records)
        require(result.returncode == status, 'installed command exit differs: ' + str(command))
        return directory / f'{index}.stdout', directory / f'{index}.stderr'

    def golden(arguments, expected):
        out, err = run(arguments)
        require(err.stat().st_size == 0 and out.read_bytes() == (fixture / expected).read_bytes(),
                'installed golden output differs: ' + expected)
        return out

    out, _ = run(['--version'])
    require(out.read_text().strip() == PRODUCT + ' ' + os.environ['PACKAGE_VERSION'], 'installed version differs')
    out, _ = run(['--help'])
    commands = out.read_text().split('Commands:\n', 1)[1].split('\n\n', 1)[0]
    require([line.split()[0] for line in commands.splitlines() if line.split()[0] != 'help'] == OPERATIONS,
            'installed stable operation scope differs')
    sorted_bed, err = run(['sort', '-'], stdin=fixture / 'sort.input.bed')
    require(err.stat().st_size == 0 and sorted_bed.read_bytes() == (fixture / 'sort.expected.bed').read_bytes(),
            'installed stdin sort differs')
    merged = directory / 'merged.bed'
    out, err = run(['--json', 'merge', sorted_bed, '--output', merged])
    result = json.loads(out.read_bytes())
    require(err.stat().st_size == 0 and result['tool'] == PRODUCT and result['status'] == 'ok' and
            result['schema_version'] == '1.0' and result['tool_version'] == os.environ['PACKAGE_VERSION'] and
            result['result'] is None, 'installed common success envelope differs')
    require(merged.read_bytes() == b'chr1\t50\t200\nchr1\t300\t400\nchr2\t50\t200\n',
            'installed sort-to-merge reload differs')
    golden(['intersect', '-a', fixture / 'a.bed', '-b', fixture / 'b.bed'], 'intersect.expected.bed')
    golden(['subtract', '-a', fixture / 'a.bed', '-b', fixture / 'b.bed'], 'subtract.expected.bed')
    golden(['complement', fixture / 'complement.input.bed', '-g', fixture / 'complement.genome.tsv'],
           'complement.expected.bed')
    golden(['cluster', fixture / 'cluster.input.bed'], 'cluster.default.expected.bed')
    golden(['window', '-w', '5', '-a', fixture / 'window.a.bed', '-b', fixture / 'window.b.bed'],
           'window.pairs.expected.bed')
    golden(['closest', '-a', fixture / 'closest.a.bed', '-b', fixture / 'closest.b.bed', '-D', 'b'],
           'closest.b-oriented.expected.bed')
    coverage_a = directory / 'coverage-a.bed'
    coverage_b1 = directory / 'coverage-b1.bed'
    coverage_b2 = directory / 'coverage-b2.bed'
    covered = directory / 'covered.bed'
    coverage_a.write_bytes(b'chr1\t0\t10\tA\t0\t+\n')
    coverage_b1.write_bytes(b'chr1\t0\t5\tfirst\n')
    coverage_b2.write_bytes(b'chr1\t3\t8\tsecond\nchr1\t3\t8\tduplicate\n')
    out, err = run(['--json', 'coverage', '-a', '-', '-b', coverage_b1, coverage_b2, '--output', covered],
                   stdin=coverage_a)
    result = json.loads(out.read_bytes())
    require(err.stat().st_size == 0 and result['tool'] == PRODUCT and result['status'] == 'ok' and
            result['schema_version'] == '1.0' and result['tool_version'] == os.environ['PACKAGE_VERSION'] and
            result['result'] is None and
            covered.read_bytes() == b'chr1\t0\t10\tA\t0\t+\t3\t8\t10\t0.8000000\n',
            'installed stdin/multi-B coverage JSON or complete output differs')
    annotated = directory / 'annotated.bed'
    out, err = run(['--json', 'annotate', '-i', '-', '--files', coverage_b1, coverage_b2,
                    '--both', '--names', 'first', 'second', '--output', annotated], stdin=coverage_a)
    result = json.loads(out.read_bytes())
    require(err.stat().st_size == 0 and result['tool'] == PRODUCT and result['status'] == 'ok' and
            result['schema_version'] == '1.0' and result['tool_version'] == os.environ['PACKAGE_VERSION'] and
            result['result'] is None and
            annotated.read_bytes() == b'#\t\t\t\t\t\tfirst_cnt\tfirst_pct\tsecond_cnt\tsecond_pct\nchr1\t0\t10\tA\t0\t+\t1\t0.500000\t2\t0.500000\n',
            'installed stdin/multi-file/names annotate JSON or complete output differs')
    mapped = directory / 'mapped.bed'
    out, err = run(['--json', 'map', '-a', '-', '-b', fixture / 'map/b.bed',
                    '-c', '4,5', '--operations', 'distinct,mean', '--delimiter', '|',
                    '--null', 'NA', '--precision', '4', '--output', mapped],
                   stdin=fixture / 'map/a.bed')
    result = json.loads(out.read_bytes())
    require(err.stat().st_size == 0 and result['tool'] == PRODUCT and result['status'] == 'ok' and
            result['schema_version'] == '1.0' and result['tool_version'] == os.environ['PACKAGE_VERSION'] and
            result['result'] is None and mapped.read_bytes() == (fixture / 'map/custom-controls.stdout').read_bytes(),
            'installed stdin/multireducer map JSON or complete output differs')
    genome_coverage = directory / 'genome-coverage.bedgraph'
    out, err = run(['--json', 'genomecov', '-i', '-', '-g', fixture / 'genomecov/genome.txt',
                    '--bedgraph-all', '--output', genome_coverage], stdin=fixture / 'genomecov/ordinary.bed')
    result = json.loads(out.read_bytes())
    require(err.stat().st_size == 0 and result['tool'] == PRODUCT and result['status'] == 'ok' and
            result['tool_version'] == os.environ['PACKAGE_VERSION'] and result['result'] is None and
            genome_coverage.read_bytes() == (fixture / 'genomecov/bedgraph-all.stdout').read_bytes(),
            'installed stdin/genomecov JSON or complete output differs')
    golden(['slop', '--input', fixture / 'coordinate/input.bed', '--genome', fixture / 'coordinate/genome.tsv',
            '--header', '--left', '3', '--right', '5', '--strand-relative'],
           'coordinate/slop-absolute.stdout')
    golden(['shift', '--input', fixture / 'coordinate/input.bed', '--genome', fixture / 'coordinate/genome.tsv',
            '--header', '--plus', '0.5', '--minus', '-2.5', '--fraction'],
           'coordinate/shift-fraction.stdout')
    protected = directory / 'protected.bed'
    protected.write_bytes(b'previous result\n')
    out, err = run(['--json', 'intersect', '-a', fixture / 'intersect.malformed-a.bed',
                    '-b', fixture / 'b.bed', '--output', protected], status=1)
    diagnostic = json.loads(err.read_bytes())
    require(out.stat().st_size == 0 and diagnostic['tool'] == PRODUCT and diagnostic['status'] == 'error' and
            diagnostic['error']['kind'] == 'InvalidInput' and diagnostic['exit_code'] == 1 and
            'greater than end' in err.read_text() and protected.read_bytes() == b'previous result\n',
            'installed malformed-input failure lost input-error diagnostics or output atomicity')
    require(len(records) == 17, 'installed lifecycle call count differs')
    save(directory / 'VERIFIED.json', dict(installed_binary=identity(binary), calls=records,
         stable_operations=OPERATIONS, stdin_sort=True, sort_output_reloaded_by_merge=True,
         coverage_stdin_pooled_json_byte_exact=True, release_acceptance=False,
         annotate_stdin_files_names_both_byte_exact=True,
         map_stdin_multireducer_json_byte_exact=True,
         genomecov_stdin_zero_runs_json_byte_exact=True,
         common_json_named_output=True, golden_outputs_byte_exact=True, failure_preserves_prior_output=True,
         files={p.name: identity(p) for p in directory.iterdir() if p.is_file()},
         scope='installed fourteen-operation candidate lifecycle; not upstream, performance or public release acceptance'))


def bundle(binary):
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    binary = Path(binary).resolve(strict=True)
    native_target = os.environ['NATIVE_TARGET']
    rustc = subprocess.check_output(['rustc', '-vV'], text=True)
    require(HOSTS[native_target] == os.environ['EXPECTED_HOST'] and
            '\nhost: ' + os.environ['EXPECTED_HOST'] + '\n' in rustc, 'bundle is not from the declared native host')
    stem = PRODUCT + '-' + os.environ['PACKAGE_VERSION'] + '-' + native_target
    folder = delivery / stem
    folder.mkdir()
    shutil.copy2(binary, folder / PRODUCT)
    for name in ('README.md', 'MIGRATION.md', 'PERFORMANCE.md', 'LICENSE-MIT', 'LICENSE-APACHE-2.0',
                 'THIRD_PARTY_LICENSES.md', 'THIRD_PARTY_NOTICES.md'):
        shutil.copy2(ROOT / name, folder / name)
    shutil.copytree(ROOT / 'licenses', folder / 'licenses')
    provenance = dict(product=PRODUCT, version=os.environ['PACKAGE_VERSION'], target=native_target,
         head_sha=os.environ['GITHUB_SHA'], run_id=int(os.environ['GITHUB_RUN_ID']),
         run_attempt=int(os.environ['GITHUB_RUN_ATTEMPT']), native_host=os.environ['EXPECTED_HOST'], rustc=rustc,
         machine=subprocess.check_output(['uname', '-a'], text=True).strip(), binary=identity(binary),
         lock_sha256=identity(ROOT / 'Cargo.lock')['sha256'],
         source_package=json.loads((delivery / 'package-verification.json').read_bytes()),
         native_dependencies=json.loads((delivery / 'native-dependencies.json').read_bytes()),
         installed_lifecycle=json.loads((delivery / 'lifecycle/VERIFIED.json').read_bytes()),
         accepted_implementation_base='02b85a1a348c271e485cca629dc4e71fa075388a',
         historical_performance_source='d7b1507b178053a087862255d84a244e4921f192',
         historical_performance_binary_sha256='9a55972f96e11bf087515b61799575505291d12d67902b4c965b8b9e32b51ff8',
         this_binary_remeasured=False, registry_publication=False, signed_release=False)
    require(provenance['source_package']['clean_vcs'] and provenance['source_package']['source_commit'] ==
            os.environ['GITHUB_SHA'], 'native bundle source package is not the clean CI head')
    require(provenance['installed_lifecycle']['installed_binary'] == identity(binary), 'bundle/lifecycle binary differs')
    save(folder / 'BUILD.json', provenance)
    paths = [folder, *folder.rglob('*')]
    require(len(paths) == len({str(p.relative_to(folder)).casefold() for p in paths}), 'case-colliding bundle path')
    for path in paths:
        require(not path.is_symlink() and (path.is_dir() or path.is_file()), 'unsupported bundle member')
        path.chmod(0o755 if path.is_dir() or path.name == PRODUCT else 0o644)
    archive = delivery / (stem + '.tar.gz')
    with tarfile.open(archive, 'w:gz') as output:
        output.add(folder, arcname=stem)
    source_archive = Path(os.environ['SOURCE_PACKAGE'])
    (delivery / 'SHA256SUMS').write_text(''.join(identity(path)['sha256'] + '  ' + path.name + '\n'
                                              for path in (archive, source_archive)))
    print(json.dumps(provenance, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('source', 'smoke', 'bundle'))
    parser.add_argument('--binary', type=Path)
    args = parser.parse_args()
    if args.action == 'source':
        source()
    elif args.action == 'smoke':
        smoke(args.binary)
    else:
        bundle(args.binary)


if __name__ == '__main__':
    main()
