# wmf/__init__.py
from importlib import import_module
from importlib.metadata import version, PackageNotFoundError

__all__ = ["cu", "models", "__version__"]

try:
    __version__ = version("wmf")
except PackageNotFoundError:
    __version__ = "0+local"

# Re-export the high-level function from wmf/wmf.py
try:
    from .wmf import read_map_raster  # noqa: F401
except Exception as _e:
    def read_map_raster(*args, **kwargs):
        raise ImportError("wmf.read_map_raster unavailable") from _e

# Eagerly try to load compiled submodules (built as wmf.cu / wmf.models)
try:
    cu = import_module(".cu", __name__)
except Exception:
    cu = None
try:
    models = import_module(".models", __name__)
except Exception:
    models = None
