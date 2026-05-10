import os
import subprocess
import sys
from pathlib import Path, PurePosixPath

from drakkar.cli_context import ERROR, INFO, RESET
from drakkar.output import print

def list_files_recursive(base_dir, exclude_snakemake=False):
    for root, dirs, files in os.walk(base_dir):
        if exclude_snakemake:
            dirs[:] = [d for d in dirs if d != ".snakemake"]
        for filename in files:
            yield Path(root) / filename

def collect_transfer_files(base_dir, args):
    selected_files = set()
    missing_paths = []
    warnings = []

    def add_annotations():
        add_file(Path("annotating/cluster_annotations.tsv.xz"))
        add_file(Path("annotating/gene_annotations.tsv.xz"))
        add_file(Path("annotating/genome_taxonomy.tsv"))

    def add_mags():
        add_dir(Path("profiling_genomes/drep/dereplicated_genomes"))

    def add_profile():
        add_file(Path("profiling_genomes/final/bases.tsv"))
        add_file(Path("profiling_genomes/final/counts.tsv"))
        add_file(Path("profiling_genomes.tsv"))

    def add_expression():
        add_file(Path("expressing/gene_counts.tsv.xz"))

    def add_bins():
        add_dir(Path("cataloging/final"))

    def add_runner_logs():
        log_files = sorted(base_dir.glob("drakkar_*.yaml"))
        if log_files:
            selected_files.update(log_files)
        else:
            warnings.append("No drakkar_*.yaml run logs found in the local directory.")

    def add_file(rel_path):
        local_path = base_dir / rel_path
        if local_path.is_file():
            selected_files.add(local_path)
        else:
            missing_paths.append(str(local_path))

    def add_dir(rel_path):
        local_path = base_dir / rel_path
        if local_path.is_dir():
            for file_path in local_path.rglob("*"):
                if file_path.is_file():
                    selected_files.add(file_path)
        else:
            missing_paths.append(str(local_path))

    if args.all:
        for file_path in list_files_recursive(base_dir):
            if file_path.is_file():
                selected_files.add(file_path)
        add_runner_logs()
        return selected_files, missing_paths, warnings

    if args.data:
        for file_path in list_files_recursive(base_dir, exclude_snakemake=True):
            if file_path.is_file():
                selected_files.add(file_path)

    if args.results:
        add_annotations()
        add_mags()
        add_profile()
        add_expression()
        add_bins()

    if args.annotations:
        add_annotations()

    if args.mags:
        add_mags()

    if args.profile:
        add_profile()

    if args.expression:
        add_expression()

    if args.bins:
        add_bins()

    add_runner_logs()

    return selected_files, missing_paths, warnings

def build_sftp_batch_commands(files, base_dir, remote_dir):
    remote_root = PurePosixPath(remote_dir)
    mkdirs = set()
    put_cmds = []

    for file_path in sorted(files):
        rel_path = file_path.relative_to(base_dir)
        remote_path = remote_root / PurePosixPath(rel_path.as_posix())
        parent = remote_path.parent
        while parent != remote_root and parent not in mkdirs:
            mkdirs.add(parent)
            parent = parent.parent
        put_cmds.append((file_path, remote_path))

    commands = []
    for directory in sorted(mkdirs, key=lambda p: len(p.parts)):
        if directory == remote_root:
            continue
        commands.append(f'-mkdir "{directory.as_posix()}"')
    for local_path, remote_path in put_cmds:
        commands.append(f'put "{local_path}" "{remote_path.as_posix()}"')

    return "\n".join(commands) + "\n", put_cmds

def run_sftp_transfer(args):
    base_dir = Path(args.local_dir).resolve()
    if args.erda:
        args.host = "io.erda.dk"
        args.user = "antton.alberdi@snm.ku.dk"
    if not args.host or not args.user:
        print(f"{ERROR}ERROR:{RESET} --host and --user are required unless --erda is set.")
        return
    sftp_cmd = ["sftp"]
    if args.port:
        sftp_cmd += ["-P", str(args.port)]
    if args.identity:
        sftp_cmd += ["-i", args.identity]
    sftp_cmd += ["-b", "-", f"{args.user}@{args.host}"]
    check_cmds = f'ls "{args.remote_dir}"\n'
    check_result = subprocess.run(sftp_cmd, input=check_cmds, text=True, capture_output=True)
    if check_result.returncode != 0:
        stderr = check_result.stderr.strip()
        if stderr:
            print(stderr, file=sys.stderr)
        print(f"{ERROR}ERROR:{RESET} Remote directory not found: {args.remote_dir}")
        return
    selected_files, missing_paths, warnings = collect_transfer_files(base_dir, args)

    if warnings:
        for warning in warnings:
            print(f"{INFO}INFO:{RESET} {warning}")

    if missing_paths:
        for missing in missing_paths:
            print(f"{ERROR}MISSING:{RESET} {missing}")

    if not selected_files:
        print(f"{ERROR}ERROR:{RESET} No files selected for transfer.")
        return

    batch_commands, put_cmds = build_sftp_batch_commands(selected_files, base_dir, args.remote_dir)
    if args.verbose:
        for local_path, remote_path in put_cmds:
            print(f"{INFO}PUT:{RESET} {local_path} -> {remote_path.as_posix()}")

    result = subprocess.run(sftp_cmd, input=batch_commands, text=True, capture_output=True)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        mkdir_only_errors = False
        if stderr:
            error_lines = [line for line in stderr.splitlines() if line.strip()]
            mkdir_only_errors = all("Couldn't create directory" in line for line in error_lines)
        if mkdir_only_errors:
            print(f"{INFO}INFO:{RESET} sftp reported existing directories; continuing.")
        else:
            if stderr:
                print(stderr, file=sys.stderr)
            print(f"{ERROR}ERROR:{RESET} sftp transfer failed with exit code {result.returncode}")
