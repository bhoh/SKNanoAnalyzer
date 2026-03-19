from __future__ import annotations

import argparse
import sys
import time
import os
from pathlib import Path
from typing import Sequence
import subprocess
import shlex
import htcondor2 as htcondor
import uproot

from rich.console import Console
from rich.table import Table
from rich.traceback import install

# Support execution via `python onnx_batch/main.py` (no package context)
try:  # pragma: no cover - runtime convenience
    from .config import (
        DEFAULT_CLASS_PATTERN,
        DEFAULT_MODEL_DIR,
        DEFAULT_OUTPUT_SUFFIX,
        DEFAULT_RECO_PATTERN,
        MAX_JETS,
        ModelConfig,
        RuntimeConfig,
    )
    from .onnx_runner import OnnxRunner, OnnxRuntimeUnavailable
    from .writer import OnnxBatchProcessor
except ImportError:  # pragma: no cover - runtime convenience
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.append(str(_Path(__file__).resolve().parent.parent))
    from onnx_batch.config import (
        DEFAULT_CLASS_PATTERN,
        DEFAULT_MODEL_DIR,
        DEFAULT_OUTPUT_SUFFIX,
        DEFAULT_RECO_PATTERN,
        MAX_JETS,
        ModelConfig,
        RuntimeConfig,
    )
    from onnx_batch.onnx_runner import OnnxRunner, OnnxRuntimeUnavailable
    from onnx_batch.writer import OnnxBatchProcessor


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch ONNX inference over SKNanoAnalyzer ROOT trees.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--mode",
        choices=["run", "submit"],
        default="run",
        help="run: process a single file locally. submit: discover inputs and submit Condor jobs directly (no .sub written).",
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        help="Input ROOT file (omit when using --condor-list/--condor-dir)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=f"Output ROOT file (default: <input>{DEFAULT_OUTPUT_SUFFIX})",
    )
    parser.add_argument(
        "--tree",
        action="append",
        help="Only process these trees (repeatable). Defaults to all TTrees.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=RuntimeConfig().chunk_size,
        help="Events to read/write per chunk.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=RuntimeConfig().batch_size,
        help="Events per ONNX batch.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads per ONNX session (sets intra/inter-op threads).",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=DEFAULT_MODEL_DIR,
        help="Directory containing the ONNX models.",
    )
    parser.add_argument(
        "--class-pattern",
        default=DEFAULT_CLASS_PATTERN,
        help="Pattern for CLASSIF models (use {fold} for formatting).",
    )
    parser.add_argument(
        "--reco-pattern",
        default=DEFAULT_RECO_PATTERN,
        help="Pattern for RECO models (use {fold} for formatting).",
    )
    parser.add_argument(
        "--channel",
        choices=["auto", "El", "Mu"],
        default="auto",
        help="Force channel override if tree name is ambiguous.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Debug helper: stop after processing this many events per tree.",
    )
    parser.add_argument(
        "--cpu-only",
        action="store_true",
        help="Disable CUDAExecutionProvider even if available.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    parser.add_argument(
        "--condor-list",
        type=Path,
        help="Text file with one input ROOT path per line for HTCondor submission.",
    )
    parser.add_argument(
        "--condor-dir",
        type=Path,
        help="Directory to scan for ROOT files (recursively) for HTCondor submission.",
    )
    parser.add_argument(
        "--condor-output-dir",
        type=Path,
        help="Output directory for Condor jobs (default: 'SPANET' directory alongside each input).",
    )
    parser.add_argument(
        "--condor-logdir",
        type=Path,
        help="Directory to store Condor log/out/err (default: condor_logs next to submit file).",
    )
    parser.add_argument(
        "--condor-submit",
        action="store_true",
        help="Immediately run condor_submit on the generated submit file.",
    )
    parser.add_argument(
        "--condor-memory",
        type=str,
        default="4GB",
        help="Memory request for Condor jobs.",
    )
    parser.add_argument(
        "--condor-cpus",
        type=int,
        default=4,
        help="CPU request for Condor jobs.",
    )
    return parser


def _get_max_entries(path: Path) -> int | None:
    """Peek ROOT file and return max entries across TTrees; None on failure."""
    try:
        with uproot.open(path, object_cache=None, array_cache=None) as f:
            counts = []
            for k in f.keys(filter_classname="TTree"):
                try:
                    counts.append(f[k].num_entries)
                except Exception:
                    continue
            return max(counts) if counts else None
    except Exception:
        return None


def main(argv: Sequence[str] | None = None) -> int:
    install()  # enable rich tracebacks with line numbers
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    console = Console(force_terminal=True, force_interactive=True)

    if args.threads < 1:
        console.print("[red]--threads must be >= 1[/]")
        return 1
    for var in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OPENBLAS_NUM_THREADS"]:
        os.environ.setdefault(var, str(args.threads))

    # Condor submission mode: generate a submit file (and optionally submit) then exit.
    if args.mode == "submit":
        if not (args.condor_list or args.condor_dir):
            console.print("[red]--mode submit requires --condor-list or --condor-dir.[/]")
            return 1
        inputs = []
        if args.condor_list:
            for line in args.condor_list.read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                inputs.append(Path(line).expanduser())
        if args.condor_dir:
            condor_dir = args.condor_dir.expanduser()
            if not condor_dir.exists():
                console.print(f"[red]Condor input directory not found:[/] {condor_dir}")
                return 1
            inputs.extend(sorted(condor_dir.glob("*.root")))
        # deduplicate while preserving order
        seen = set()
        unique_inputs = []
        for p in inputs:
            if p in seen:
                continue
            seen.add(p)
            unique_inputs.append(p)
        inputs = unique_inputs
        if not inputs:
            console.print(
                "[red]No inputs found for Condor submission "
                "(check --condor-list or --condor-dir).[/]"
            )
            return 1

        repo_root = Path(__file__).resolve().parents[2]
        base_submit = args.condor_list if args.condor_list else args.condor_dir
        logdir = (
            args.condor_logdir
            if args.condor_logdir
            else base_submit.with_name(base_submit.stem + "_logs")
        )
        logdir.mkdir(parents=True, exist_ok=True)

        extra_args = []
        if args.tree:
            for t in args.tree:
                extra_args += ["--tree", t]
        if args.chunk_size != RuntimeConfig().chunk_size:
            extra_args += ["--chunk-size", str(args.chunk_size)]
        if args.batch_size != RuntimeConfig().batch_size:
            extra_args += ["--batch-size", str(args.batch_size)]
        if args.model_dir != DEFAULT_MODEL_DIR:
            extra_args += ["--model-dir", str(args.model_dir)]
        if args.class_pattern != DEFAULT_CLASS_PATTERN:
            extra_args += ["--class-pattern", args.class_pattern]
        if args.reco_pattern != DEFAULT_RECO_PATTERN:
            extra_args += ["--reco-pattern", args.reco_pattern]
        if args.channel != "auto":
            extra_args += ["--channel", args.channel]
        if args.max_events is not None:
            extra_args += ["--max-events", str(args.max_events)]
        if args.cpu_only:
            extra_args += ["--cpu-only"]
        if args.force:
            extra_args += ["--force"]

        env_parts = [
            "PYTHONPATH=python",
            "OMP_NUM_THREADS=$(threads)",
            "MKL_NUM_THREADS=$(threads)",
            "NUMEXPR_NUM_THREADS=$(threads)",
            "OPENBLAS_NUM_THREADS=$(threads)",
            "ORT_INTRA_OP_NUM_THREADS=$(threads)",
            "ORT_INTER_OP_NUM_THREADS=1",
            "PYTHONUNBUFFERED=1",
        ]

        cmd_tokens = [
            str(repo_root / "python/onnx_batch/main.py"),
            "--mode",
            "run",
            "$(input)",
            "--output",
            "$(dst)",
            "--threads",
            "$(threads)",
            *extra_args,
        ]

        submit_ad = htcondor.Submit(
            {
                "universe": "vanilla",
                "getenv": "True",
                "request_cpus": "$(cpus)",
                "request_memory": "$(mem)",
                "initialdir": str(repo_root),
                "executable": sys.executable,
                "environment": ";".join(env_parts),
                "log": f"{logdir}/$(ClusterId).log",
                "output": "$(ClusterId).$(ProcId).out",
                "error": "$(ClusterId).$(ProcId).err",
                "should_transfer_files": "YES",
                "when_to_transfer_output": "ON_EXIT_OR_EVICT",
                "stream_output": "True",
                "stream_error": "True",
                "transfer_input_files": "",
                "transfer_output_remaps": f"$(ClusterId).$(ProcId).out={logdir}/$(ClusterId).$(ProcId).out;$(ClusterId).$(ProcId).err={logdir}/$(ClusterId).$(ProcId).err",
                "arguments": " ".join(cmd_tokens),
            }
        )

        itemdata = []
        for inp in inputs:
            if not inp.exists():
                console.print(f"[yellow]Warning: input not found now (will rely on worker): {inp}[/]")
            if args.condor_output_dir:
                out_path = args.condor_output_dir / f"{inp.stem}{DEFAULT_OUTPUT_SUFFIX}"
            else:
                out_path = args.condor_dir / "SPANET" / f"{inp.stem}{DEFAULT_OUTPUT_SUFFIX}"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            # Auto-scale resources based on entries; fallback to user defaults if counting fails (e.g., backup cycle/offline).
            entries = _get_max_entries(inp)
            if entries is None:
                job_cpus = args.condor_cpus
                job_mem = args.condor_memory
            else:
                min_e, max_e = 1000, 15_340_148
                min_c, max_c = 4, 64
                min_mem, max_mem = 16, 220  # GB
                if entries <= min_e:
                    frac = 0.0
                elif entries >= max_e:
                    frac = 1.0
                else:
                    frac = (entries - min_e) / (max_e - min_e)
                job_cpus = int(round(min_c + frac * (max_c - min_c)))
                job_mem_val = min_mem + frac * (max_mem - min_mem)
                job_mem = f"{int(job_mem_val + 0.999)}GB"
            # To avoid Condor killing jobs for CPU overuse, keep thread count slightly below requested cpus.
            job_threads = max(1, int(job_cpus * 0.8))
            itemdata.append(
                {
                    "input": str(inp),
                    "dst": str(out_path),
                    "cpus": str(job_cpus),
                    "threads": str(job_threads),
                    "mem": str(job_mem),
                }
            )

        schedd = htcondor.Schedd()
        try:
            with schedd.transaction() as txn:
                # queue_with_itemdata expects count per itemdata row; use 1 to avoid N^2 submissions.
                result = submit_ad.queue_with_itemdata(txn, 1, iter(itemdata))
            cluster_id = result.cluster() if hasattr(result, "cluster") else "?"
            console.print(f"[green]Submitted {len(itemdata)} jobs.[/] Cluster: {cluster_id}")
        except Exception as exc:  # pragma: no cover - runtime safety
            console.print(f"[red]Condor submission failed:[/] {exc}")
            return 1
        return 0

    if args.input is None:
        console.print("[red]Input file is required unless using --condor-list/--condor-dir.[/]")
        return 1

    input_path = args.input.expanduser()
    if not input_path.exists():
        console.print(f"[red]Input file not found:[/] {input_path}")
        return 1

    output_path = (
        args.output.expanduser()
        if args.output
        else input_path.with_name(f"{input_path.stem}{DEFAULT_OUTPUT_SUFFIX}")
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and not args.force:
        console.print(
            f"[red]Output file already exists:[/] {output_path}. "
            "Use --force to overwrite."
        )
        return 1

    model_cfg = ModelConfig(
        model_dir=args.model_dir,
        class_pattern=args.class_pattern,
        reco_pattern=args.reco_pattern,
        prefer_cuda=not args.cpu_only,
    )
    runtime_cfg = RuntimeConfig(
        batch_size=args.batch_size,
        chunk_size=args.chunk_size,
        max_events=args.max_events,
        channel_hint=args.channel,
    )

    try:
        runner = OnnxRunner(model_cfg, console, threads=args.threads)
    except OnnxRuntimeUnavailable as exc:
        console.print(f"[red]{exc}[/]")
        return 1
    except Exception as exc:  # pragma: no cover - runtime safety
        console.print(f"[red]Failed to initialize ONNX runtime:[/] {exc}")
        return 1

    processor = OnnxBatchProcessor(runner, runtime_cfg, console, max_jets=MAX_JETS)

    console.rule("[bold green]Starting batched ONNX inference")
    start = time.perf_counter()
    try:
        tree_stats = processor.process_file(input_path, output_path, args.tree)
    except Exception as exc:  # pragma: no cover - runtime safety
        console.print(f"[red]Processing failed:[/] {exc}")
        console.print_exception()
        return 1
    duration = time.perf_counter() - start

    table = Table(title="ONNX batch summary", show_lines=False)
    table.add_column("Tree")
    table.add_column("Events", justify="right")
    table.add_column("Inferred", justify="right")
    table.add_column("Valid", justify="right")
    table.add_column("Walltime [s]", justify="right")

    for stat in tree_stats:
        table.add_row(
            stat.name,
            f"{stat.entries}",
            f"{stat.inferred}",
            f"{stat.valid}",
            f"{stat.duration:.2f}",
        )

    console.print(table)
    console.print(
        f"[bold green]Output:[/] {output_path} | "
        f"Provider order: {', '.join(runner.provider_info.requested)} | "
        f"Models loaded: {len(runner.loaded_models)} | "
        f"Elapsed: {duration:.2f}s"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
