from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

_SCHEMA_EXT = ".json"


@dataclass(slots=True)
class ToolSchema:
    id: str
    payload: dict
    source_path: Path
    category: Optional[str] = None


class ToolRegistry:
    def __init__(self, schema_dirs: Iterable[Path]):
        self._schema_dirs = [Path(directory).expanduser().resolve() for directory in schema_dirs]
        self._schemas: Dict[str, ToolSchema] = {}
        self.refresh()

    @staticmethod
    def _strip_json_comments(raw: str) -> str:
        lines: List[str] = []
        for line in raw.splitlines():
            stripped = line.lstrip()
            if stripped.startswith("//"):
                continue
            lines.append(line)
        return "\n".join(lines)

    def refresh(self) -> None:
        schemas: Dict[str, ToolSchema] = {}
        for directory in self._schema_dirs:
            if not directory.exists():
                continue
            for path in directory.rglob(f"*{_SCHEMA_EXT}"):
                try:
                    with path.open("r", encoding="utf-8") as handle:
                        raw = handle.read()
                    payload = json.loads(self._strip_json_comments(raw))
                except json.JSONDecodeError as exc:
                    raise ValueError(f"Invalid JSON schema: {path}") from exc
                tool_id = payload.get("id")
                if not tool_id:
                    raise ValueError(f"Schema missing 'id': {path}")
                if tool_id in schemas:
                    raise ValueError(f"Duplicate tool id '{tool_id}' from {path}")
                category = payload.get("category")
                schemas[tool_id] = ToolSchema(id=tool_id, payload=payload, source_path=path, category=category)
        self._schemas = schemas

    def list(self) -> List[ToolSchema]:
        return list(self._schemas.values())

    def list_by_category(self, category: str) -> List[ToolSchema]:
        return [schema for schema in self._schemas.values() if schema.category == category]

    def get(self, tool_id: str) -> Optional[dict]:
        schema = self._schemas.get(tool_id)
        return schema.payload if schema else None

    def add_schema_dir(self, directory: Path) -> None:
        path = Path(directory).expanduser().resolve()
        if path not in self._schema_dirs:
            self._schema_dirs.append(path)
            self.refresh()