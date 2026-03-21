import os
from pathlib import Path


def run_runfile(model, runfile_path, command_callback=None):
    """Execute a legacy runfile, honouring nested ``call`` directives."""

    runfile_path = Path(runfile_path).resolve()
    original_dir = os.getcwd()
    # The classic interpreter expects to run relative to the runfile's
    # directory, so temporarily adjust the working directory.
    os.chdir(runfile_path.parent)

    def execute(path):
        with open(path, "r") as handle:
            for raw_line in handle:
                stripped = raw_line.strip()
                if not stripped:
                    continue
                lowered = stripped.lower()
                if lowered.startswith("echo"):
                    continue
                if lowered.startswith("call "):
                    # ``call`` directives may be nested, so recurse on the
                    # referenced file.
                    target = stripped.split(None, 1)[1]
                    execute((Path(path).parent / target).resolve())
                    continue
                # The Fortran-style runfiles treat leading ``c`` or ``C`` as a
                # comment indicator, which modern ``procline`` ignores.
                if stripped[0].lower() == "c" and (len(stripped) == 1 or stripped[1].isspace()):
                    continue
                if command_callback is not None:
                    command_callback(stripped)
                model.procline(stripped)

    try:
        execute(runfile_path)
    finally:
        os.chdir(original_dir)
    return model
