from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

__all__ = [
    "configure_log_directory",
    "open_log",
    "emit_log",
    "close_log",
]

_LOGGER = logging.getLogger("nlsl")
_HANDLER: Optional[logging.Handler] = None
_LOG_DIR = Path.home() / ".nlsl"


def configure_log_directory(directory: Path | str) -> None:
    """Configure the base directory used for log files."""
    global _LOG_DIR
    _LOG_DIR = Path(directory)
    _LOG_DIR.mkdir(parents=True, exist_ok=True)


def _remove_handler() -> None:
    global _HANDLER
    if _HANDLER is not None:
        _LOGGER.removeHandler(_HANDLER)
        _HANDLER.close()
        _HANDLER = None


def open_log(filename: str) -> str:
    """Install a file handler for *filename* under the configured directory."""
    global _HANDLER
    _remove_handler()
    path = Path(filename)
    if not path.is_absolute():
        path = _LOG_DIR / path
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(message)s"))
    _LOGGER.addHandler(handler)
    _LOGGER.propagate = False
    _HANDLER = handler
    return str(path)


def emit_log(message: str) -> None:
    """Emit *message* to the configured log handler."""
    _LOGGER.debug(message)


def close_log() -> None:
    """Remove the installed handler, if any."""
    _remove_handler()

