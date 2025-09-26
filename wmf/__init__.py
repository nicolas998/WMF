from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

# Version
try:
    __version__ = version("wmf")
except PackageNotFoundError:
    __version__ = "0+local"

# Keep old usage working: from wmf import wmf
from . import wmf as wmf  # noqa: F401

# Try to import compiled submodules so: from wmf import cu, models works
try:
    cu = import_module(".cu", __name__)
except Exception:
    cu = None

try:
    models = import_module(".models", __name__)
except Exception:
    models = None

# Re-export public API from wmf.wmf at package top-level
_public = getattr(wmf, "__all__", None)
if _public is None:
    _public = [n for n in dir(wmf) if not n.startswith("_")]
for _name in _public:
    globals()[_name] = getattr(wmf, _name)
__all__ = ["__version__", "wmf", "cu", "models"] + list(_public)