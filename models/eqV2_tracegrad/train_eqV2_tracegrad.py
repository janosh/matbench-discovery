"""For backwards compatibility"""

from __future__ import annotations

import sys

sys.path.insert(0, "/home/helxgroup/yinshi/fywang/code/fairchem-tracegrad/src")

if __name__ == "__main__":
    from fairchem.core._cli import main

    main()
