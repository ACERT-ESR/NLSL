import multiprocessing
import os
import sys
import types

class _QtDummy(types.ModuleType):
    def __getattr__(self, name):
        placeholder = types.SimpleNamespace()
        setattr(self, name, placeholder)
        return placeholder

def qt_signal_stub(*args, **kwargs):
    class _SignalPlaceholder:
        def __call__(self, *call_args, **call_kwargs):
            return None
    return _SignalPlaceholder()

def qt_slot_stub(*args, **kwargs):
    def decorator(fn):
        return fn
    return decorator

def qt_property_stub(*args, **kwargs):
    class _PropertyPlaceholder:
        def __call__(self, *call_args, **call_kwargs):
            return None
    return _PropertyPlaceholder()

class _QtValueConstant:
    def __init__(self, value=0):
        self.value = value
    def __int__(self):
        return int(self.value)

class _QtDummyConstants:
    def __getattr__(self, name):
        return _QtValueConstant()

class _QtConstants:
    def __init__(self):
        self.Key = _QtDummyConstants()
        self.KeyboardModifier = _QtDummyConstants()
        self.CursorShape = _QtDummyConstants()
        self.MouseButton = _QtDummyConstants()
    def __getattr__(self, name):
        return 0

# pyspecdata pulls in Qt bindings that are not available in the test environment.
# We provide lightweight stubs so that the import path succeeds without a GUI.
pyqt_parent = types.ModuleType("PyQt5")
pyqt_parent.QtWidgets = _QtDummy("PyQt5.QtWidgets")
pyqt_parent.QtCore = _QtDummy("PyQt5.QtCore")
pyqt_parent.QtGui = _QtDummy("PyQt5.QtGui")
pyqt_parent.__version__ = "6.5.0"
class _QtWidgetBase:
    def __init__(self, *args, **kwargs):
        pass
pyqt_parent.QtWidgets.QPushButton = _QtWidgetBase
pyqt_parent.QtWidgets.QHBoxLayout = _QtWidgetBase
pyqt_parent.QtWidgets.QGridLayout = _QtWidgetBase
pyqt_parent.QtWidgets.QWidget = _QtWidgetBase
pyqt_parent.QtWidgets.QDialog = _QtWidgetBase
pyqt_parent.QtWidgets.QMainWindow = _QtWidgetBase
pyqt_parent.QtWidgets.QToolBar = _QtWidgetBase
class _FakeLibraryInfo:
    @staticmethod
    def version():
        class _FakeVersion:
            def segments(self):
                return (6, 5, 0)
        return _FakeVersion()
pyqt_parent.QtCore.QLibraryInfo = _FakeLibraryInfo
pyqt_parent.QtCore.Signal = qt_signal_stub
pyqt_parent.QtCore.Slot = qt_slot_stub
pyqt_parent.QtCore.Property = qt_property_stub
pyqt_parent.QtCore.Qt = _QtConstants()
class _QtColor:
    def __init__(self, *args, **kwargs):
        pass
pyqt_parent.QtGui.QColor = _QtColor
sys.modules["PyQt5"] = pyqt_parent
sys.modules["PyQt5.QtWidgets"] = pyqt_parent.QtWidgets
sys.modules["PyQt5.QtCore"] = pyqt_parent.QtCore
sys.modules["PyQt5.QtGui"] = pyqt_parent.QtGui

pyside_parent = types.ModuleType("PySide6")
pyside_parent.QtWidgets = _QtDummy("PySide6.QtWidgets")
pyside_parent.QtCore = _QtDummy("PySide6.QtCore")
pyside_parent.QtGui = _QtDummy("PySide6.QtGui")
pyside_parent.__version__ = "6.5.0"
pyside_parent.QtWidgets.QPushButton = _QtWidgetBase
pyside_parent.QtWidgets.QHBoxLayout = _QtWidgetBase
pyside_parent.QtWidgets.QGridLayout = _QtWidgetBase
pyside_parent.QtWidgets.QWidget = _QtWidgetBase
pyside_parent.QtWidgets.QDialog = _QtWidgetBase
pyside_parent.QtWidgets.QMainWindow = _QtWidgetBase
pyside_parent.QtWidgets.QToolBar = _QtWidgetBase
pyside_parent.QtCore.QLibraryInfo = _FakeLibraryInfo
pyside_parent.QtCore.Signal = qt_signal_stub
pyside_parent.QtCore.Slot = qt_slot_stub
pyside_parent.QtCore.Property = qt_property_stub
pyside_parent.QtCore.Qt = _QtConstants()
pyside_parent.QtGui.QColor = _QtColor
sys.modules["PySide6"] = pyside_parent
sys.modules["PySide6.QtWidgets"] = pyside_parent.QtWidgets
sys.modules["PySide6.QtCore"] = pyside_parent.QtCore
sys.modules["PySide6.QtGui"] = pyside_parent.QtGui
sys.modules["shiboken6"] = types.ModuleType("shiboken6")

from examples import pyspec_differential_evolution as pde

# We force the spawn start method so that every process builds a fresh interpreter.
multiprocessing.set_start_method("spawn", force=True)

def _set_and_return_gxx(queue, new_value):
    ctx = pde._get_ctx()
    ctx["n"]["gxx"] = new_value
    queue.put((os.getpid(), float(ctx["n"]["gxx"])))

def _read_gxx(queue):
    ctx = pde._get_ctx()
    queue.put((os.getpid(), float(ctx["n"]["gxx"])))

def test_nlsl_instances_are_process_isolated():
    # Two workers set distinct values and report back. If the processes shared an
    # NLSL instance we would see interference between the returned values.
    ctx = multiprocessing.get_context("spawn")
    queue = ctx.Queue()

    proc_one = ctx.Process(target=_set_and_return_gxx, args=(queue, 2.5))
    proc_two = ctx.Process(target=_set_and_return_gxx, args=(queue, 2.7))
    proc_one.start()
    proc_two.start()
    proc_one.join()
    proc_two.join()

    results = [queue.get(), queue.get()]
    pids = {item[0] for item in results}
    values = {item[1] for item in results}

    assert len(pids) == 2
    assert values == {2.5, 2.7}

    # A third process should see the initial gxx value, demonstrating that state
    # mutations do not leak between processes.
    proc_read = ctx.Process(target=_read_gxx, args=(queue,))
    proc_read.start()
    proc_read.join()
    read_pid, read_value = queue.get()

    assert read_pid not in pids
    assert read_value == float(pde.initial_params["gxx"])
