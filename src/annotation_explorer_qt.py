from __future__ import annotations

import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BACTolomew Annotation Explorer")
    parser.add_argument(
        "--dev",
        action="store_true",
        help="Skip Bakta env/DB setup and use current hardcoded dev Bakta env+DB.",
    )
    args = parser.parse_args()
    from annotation_explorer_polars.app import main

    sys.exit(main(dev_mode=bool(args.dev)))
