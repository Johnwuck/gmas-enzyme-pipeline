"""允许 `python -m src.pipeline`。"""

from __future__ import annotations

from .pipeline import main

if __name__ == "__main__":
    raise SystemExit(main())
