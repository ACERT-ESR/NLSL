import os

# Ensure headless Qt operation without requiring callers to set environment
# variables explicitly when running the test suite.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
