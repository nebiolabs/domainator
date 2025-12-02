from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

from .app import create_app
from .config import ServerConfig, DEFAULT_DATA_DIR


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Start the Domainator web server",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--port", type=int, default=0, help="Port to bind the server. By default, a random free port is chosen.")
    parser.add_argument("--host", default="127.0.0.1", help="Host interface to bind. Defaults to localhost. To make the server accessible on the network, use 0.0.0.0")
    parser.add_argument("--data-dir", type=Path, default=None, help=f"Data directory for uploads, jobs, and logs. Defaults to {DEFAULT_DATA_DIR}")
    parser.add_argument(
        "--schema-dir",
        type=Path,
        action="append",
        default=None,
        help="Directory containing tool and workflow schema JSON files",
    )
    parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    schema_dirs = [path.expanduser().resolve() for path in args.schema_dir] if args.schema_dir else None
    config = ServerConfig(data_dir=args.data_dir, debug=args.debug, schema_dirs=schema_dirs)
    app = create_app(config)

    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


if __name__ == "__main__":
    sys.exit(main())
