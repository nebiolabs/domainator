"""Split a sequence or profile file into connected components derived from a matrix.

Given one input file and one symmetric matrix, threshold the matrix at the
requested lower bound, find connected components, and write one output file per
cluster. Clusters smaller than the requested minimum size are combined into a
single small_clusters output file.
"""

from collections import Counter
from contextlib import ExitStack
from jsonargparse import ActionConfigFile, ArgumentParser
from os import PathLike
from pathlib import Path
import sys
from typing import Dict, Iterable, Iterator, List, Optional, Tuple, Union

import pyhmmer

from domainator import __version__, RawAndDefaultsFormatter
from domainator import build_ssn
from domainator.Bio import SeqIO
from domainator.data_matrix import DataMatrix
from domainator import select_by_cds
from domainator.utils import DomainatorCDS, get_file_type, parse_seqfiles, pyhmmer_decode, write_genbank


SUPPORTED_INPUT_TYPES = {"genbank", "fasta", "hmm", "cm"}


def _positive_int(value: Union[str, int]) -> int:
    value = int(value)
    if value < 1:
        raise ValueError("--min_cluster_size must be a positive integer")
    return value


def _get_supported_input_type(input_path: Union[str, PathLike]) -> str:
    input_type = get_file_type(input_path)
    if input_type not in SUPPORTED_INPUT_TYPES:
        raise ValueError(
            f"Unsupported input type for {input_path}. Supported file types are: genbank, fasta, hmm, and cm."
        )
    return input_type


def _get_output_extension(input_path: Union[str, PathLike]) -> str:
    extension = Path(input_path).suffix
    if extension == "":
        raise ValueError("Input file must have an extension so output files can preserve the input format.")
    return extension


def _validate_unique_labels(labels: Iterable[str], label_source: str) -> List[str]:
    ordered_labels = list(labels)
    seen = set()
    duplicates = set()
    for label in ordered_labels:
        if label in seen:
            duplicates.add(label)
        seen.add(label)
    if duplicates:
        duplicates_text = ", ".join(sorted(duplicates))
        raise ValueError(f"{label_source} contain duplicate labels: {duplicates_text}")
    return ordered_labels


def _iter_hmm_models(input_path: Union[str, PathLike]) -> Iterator[Tuple[str, pyhmmer.plan7.HMM]]:
    with pyhmmer.plan7.HMMFile(input_path) as hmm_file:
        for model in hmm_file:
            yield pyhmmer_decode(model.name), model


def _iter_cm_blocks(input_path: Union[str, PathLike]) -> Iterator[Tuple[str, str]]:
    record_lines: List[str] = []
    current_name: Optional[str] = None

    with open(input_path, "r", encoding="utf-8") as handle:
        for line in handle:
            record_lines.append(line)
            if line.startswith("NAME"):
                parts = line.split(None, 1)
                if len(parts) != 2:
                    raise ValueError(f"Malformed NAME line in CM file: {line.rstrip()}")
                current_name = parts[1].strip()
            if line.strip() == "//":
                if current_name is None:
                    raise ValueError("Encountered a CM record without a NAME field.")
                yield current_name, "".join(record_lines)
                record_lines = []
                current_name = None

    if record_lines:
        raise ValueError("CM file ended before a record terminator ('//').")


def _read_input_labels(input_path: Union[str, PathLike], input_type: str) -> List[str]:
    if input_type in {"genbank", "fasta"}:
        return _validate_unique_labels(
            (record.id for record in parse_seqfiles([input_path], filetype_override=input_type)),
            "Input records",
        )
    if input_type == "hmm":
        return _validate_unique_labels((name for name, _ in _iter_hmm_models(input_path)), "Input HMM profiles")
    return _validate_unique_labels((name for name, _ in _iter_cm_blocks(input_path)), "Input CM profiles")


def _load_cluster_membership(matrix_path: Union[str, PathLike], lb: float) -> Tuple[DataMatrix, Dict[str, int], Counter]:
    matrix = DataMatrix.from_file(matrix_path, lower_bound=lb)
    if not matrix.symmetric_labels:
        raise ValueError("Input matrix must have symmetric row and column labels.")

    matrix_labels = _validate_unique_labels(matrix.rows, "Matrix rows")
    if len(matrix_labels) == 0:
        raise ValueError("Input matrix does not contain any labels.")

    cluster_labels = build_ssn.cluster_labels_from_graph(matrix, lb)
    label_to_cluster = {label: int(cluster_labels[idx]) for idx, label in enumerate(matrix.rows)}
    return matrix, label_to_cluster, Counter(label_to_cluster.values())


def _validate_label_reconciliation(matrix_labels: Iterable[str], input_labels: Iterable[str]) -> None:
    matrix_label_set = set(matrix_labels)
    input_label_set = set(input_labels)

    matrix_only = sorted(matrix_label_set.difference(input_label_set))
    input_only = sorted(input_label_set.difference(matrix_label_set))
    if matrix_only or input_only:
        details = []
        if matrix_only:
            details.append(f"matrix-only labels: {', '.join(matrix_only[:10])}")
        if input_only:
            details.append(f"input-only labels: {', '.join(input_only[:10])}")
        raise ValueError("Matrix labels and input labels do not match exactly (" + "; ".join(details) + ").")


def _plan_output_paths(
    outdir: Union[str, PathLike],
    extension: str,
    cluster_sizes: Counter,
    min_cluster_size: int,
    write_small_clusters: bool,
) -> Tuple[Dict[int, Path], Optional[Path]]:
    output_dir = Path(outdir)
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"Output path exists and is not a directory: {outdir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_paths: Dict[int, Path] = dict()
    for cluster_id, size in sorted(cluster_sizes.items()):
        if size >= min_cluster_size:
            cluster_paths[cluster_id] = output_dir / f"{cluster_id}{extension}"

    small_clusters_path = None
    if write_small_clusters and any(size < min_cluster_size for size in cluster_sizes.values()):
        small_clusters_path = output_dir / f"small_clusters{extension}"

    for output_path in list(cluster_paths.values()) + ([small_clusters_path] if small_clusters_path is not None else []):
        if output_path.exists():
            raise FileExistsError(f"Refusing to overwrite existing file: {output_path}")

    return cluster_paths, small_clusters_path


def _get_destination_handle(
    cluster_id: int,
    cluster_paths: Dict[int, Path],
    regular_handles: Dict[int, object],
    small_clusters_handle: Optional[object],
) -> Optional[object]:
    if cluster_id in cluster_paths:
        return regular_handles[cluster_id]
    if small_clusters_handle is None:
        return None
    return small_clusters_handle


def _center_record_on_search_hit(record):
    cdss = DomainatorCDS.list_from_contig(record, include_nucleic_acid_annotations=True)
    cdss.sort(key=lambda cds: cds.feature.location.stranded_start)

    for focus_index, cds in enumerate(cdss):
        if cds.domain_search_feature is not None:
            centered_record, _ = select_by_cds.get_cds_neighborhood(
                record,
                cdss,
                focus_index,
                whole_contig=True,
            )
            return centered_record

    raise ValueError(
        f"Record {record.id} does not contain any Domain_Search annotations required by --pad_on_search_hits."
    )


def split_by_cluster(
    input_path: Union[str, PathLike],
    matrix_path: Union[str, PathLike],
    outdir: Union[str, PathLike],
    lb: float = 0,
    min_cluster_size: int = 1,
    write_small_clusters: bool = False,
    pad_on_search_hits: bool = False,
) -> Dict[str, Union[Dict[int, Path], Optional[Path]]]:
    input_type = _get_supported_input_type(input_path)
    extension = _get_output_extension(input_path)

    if pad_on_search_hits and input_type != "genbank":
        raise ValueError("--pad_on_search_hits is only supported for GenBank inputs.")

    matrix, label_to_cluster, cluster_sizes = _load_cluster_membership(matrix_path, lb)
    input_labels = _read_input_labels(input_path, input_type)
    _validate_label_reconciliation(matrix.rows, input_labels)

    cluster_paths, small_clusters_path = _plan_output_paths(
        outdir,
        extension,
        cluster_sizes,
        min_cluster_size,
        write_small_clusters,
    )
    open_mode = "wb" if input_type == "hmm" else "w"

    with ExitStack() as stack:
        regular_handles = {
            cluster_id: stack.enter_context(open(path, open_mode, encoding=None if open_mode == "wb" else "utf-8"))
            for cluster_id, path in cluster_paths.items()
        }
        small_clusters_handle = None
        if small_clusters_path is not None:
            small_clusters_handle = stack.enter_context(
                open(small_clusters_path, open_mode, encoding=None if open_mode == "wb" else "utf-8")
            )

        if input_type in {"genbank", "fasta"}:
            padded_cluster_records = None
            padded_small_cluster_records = None
            if input_type == "genbank" and pad_on_search_hits:
                padded_cluster_records = {cluster_id: list() for cluster_id in cluster_paths}
                padded_small_cluster_records = [] if small_clusters_path is not None else None

            for record in parse_seqfiles([input_path], filetype_override=input_type):
                cluster_id = label_to_cluster[record.id]
                if input_type == "genbank" and pad_on_search_hits:
                    centered_record = _center_record_on_search_hit(record)
                    if cluster_id in cluster_paths:
                        padded_cluster_records[cluster_id].append(centered_record)
                    elif padded_small_cluster_records is not None:
                        padded_small_cluster_records.append(centered_record)
                    continue

                destination_handle = _get_destination_handle(cluster_id, cluster_paths, regular_handles, small_clusters_handle)
                if destination_handle is None:
                    continue
                if input_type == "genbank":
                    write_genbank([record], destination_handle)
                else:
                    SeqIO.write(record, destination_handle, "fasta")

            if input_type == "genbank" and pad_on_search_hits:
                for cluster_id, records in padded_cluster_records.items():
                    select_by_cds.pad_records(records)
                    write_genbank(records, regular_handles[cluster_id])
                if padded_small_cluster_records is not None:
                    select_by_cds.pad_records(padded_small_cluster_records)
                    write_genbank(padded_small_cluster_records, small_clusters_handle)
        elif input_type == "hmm":
            for model_name, model in _iter_hmm_models(input_path):
                cluster_id = label_to_cluster[model_name]
                destination_handle = _get_destination_handle(cluster_id, cluster_paths, regular_handles, small_clusters_handle)
                if destination_handle is None:
                    continue
                model.write(destination_handle)
        else:
            for model_name, cm_text in _iter_cm_blocks(input_path):
                cluster_id = label_to_cluster[model_name]
                destination_handle = _get_destination_handle(cluster_id, cluster_paths, regular_handles, small_clusters_handle)
                if destination_handle is None:
                    continue
                destination_handle.write(cm_text)

    return {"clusters": cluster_paths, "small_clusters": small_clusters_path}


def main(argv):
    parser = ArgumentParser(description=f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Input genbank, fasta, hmm, or cm file to split by cluster.",
    )
    parser.add_argument(
        "--matrix",
        required=True,
        type=str,
        help="Matrix input in dense HDF5, sparse HDF5, or dense text format.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=str,
        help="Directory to write cluster files into. Existing target files cause an error.",
    )
    parser.add_argument(
        "--lb",
        type=float,
        default=0,
        help="Delete edges with weights less than or equal to this value before finding connected components.",
    )
    parser.add_argument(
        "--min_cluster_size",
        type=_positive_int,
        default=1,
        help="Minimum cluster size required to write a dedicated cluster file.",
    )
    parser.add_argument(
        "--write_small_clusters",
        action="store_true",
        default=False,
        help="If set, write clusters smaller than --min_cluster_size into small_clusters.<ext>. Otherwise those records or models are skipped.",
    )
    parser.add_argument(
        "--pad_on_search_hits",
        action="store_true",
        default=False,
        help="For GenBank inputs, center each contig on its first Domain_Search annotation and pad records within each output file so the search hits align. Raises an error if any record lacks Domain_Search.",
    )
    parser.add_argument("--config", action=ActionConfigFile)

    params = parser.parse_args(argv)
    split_by_cluster(
        input_path=params.input,
        matrix_path=params.matrix,
        outdir=params.outdir,
        lb=params.lb,
        min_cluster_size=params.min_cluster_size,
        write_small_clusters=params.write_small_clusters,
        pad_on_search_hits=params.pad_on_search_hits,
    )


def _entrypoint():
    main(sys.argv[1:])


if __name__ == "__main__":
    main(sys.argv[1:])