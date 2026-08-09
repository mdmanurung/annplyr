"""Run ASV while keeping its machine registry in ignored repository state."""

from __future__ import annotations

from pathlib import Path

from asv.machine import MachineCollection

MachineCollection.get_machine_file_path = staticmethod(lambda: str(Path.cwd() / ".asv-machine.json"))


if __name__ == "__main__":
    from asv.main import main

    main()
