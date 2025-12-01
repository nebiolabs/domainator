from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Mock genome annotation workflow")
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)

    with open(args.config, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    output_path = Path(payload.get("output", "workflow_result.txt"))
    output_path.write_text("workflow complete\n")
    return 0


if __name__ == "__main__":
    main()
