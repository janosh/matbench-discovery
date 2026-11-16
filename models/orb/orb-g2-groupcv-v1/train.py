"""Reproduce the frozen training state from the verified release data archive.

Usage: python train.py --root /path/to/unpacked/source
Optional: --verify-backbone --base-checkpoint /path/to/orb.ckpt
"""

from orb_g2_cv.reproduce import main

if __name__ == "__main__":
    main()
