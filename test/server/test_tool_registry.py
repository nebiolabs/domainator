from __future__ import annotations

import json
import pytest
from pathlib import Path

from domainator.server.tool_registry import ToolRegistry


def write_schema(directory: Path, tool_id: str, entry_point: str = "tool.py") -> Path:
    payload = {
        "id": tool_id,
        "runner": "python",
        "entry_point": entry_point,
        "parameters": [],
    }
    path = directory / f"{tool_id}.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    return path


def test_tool_registry_loads_schema(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schemas"
    schema_dir.mkdir()
    write_schema(schema_dir, "example")

    registry = ToolRegistry([schema_dir])
    schema = registry.get("example")
    assert schema is not None
    assert schema["id"] == "example"


def test_tool_registry_duplicate_ids_raise(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schemas"
    schema_dir.mkdir()
    write_schema(schema_dir, "duplicate")
    second = schema_dir / "other.json"
    with second.open("w", encoding="utf-8") as handle:
        json.dump({"id": "duplicate", "entry_point": "tool.py"}, handle)

    with pytest.raises(ValueError):
        ToolRegistry([schema_dir])


def test_tool_registry_ignores_comment_lines(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schemas"
    schema_dir.mkdir()
    schema_path = schema_dir / "commented.json"
    schema_path.write_text(
        "// leading comment\n"
        "{\n"
        "  \"id\": \"commented\",\n"
        "  // mid comment\n"
        "  \"runner\": \"python\",\n"
        "  \"entry_point\": \"tool.py\",\n"
        "  \"parameters\": []\n"
        "}\n",
        encoding="utf-8",
    )

    registry = ToolRegistry([schema_dir])
    schema = registry.get("commented")
    assert schema is not None
    assert schema["runner"] == "python"
