from __future__ import annotations

from pipeline.CLI import read_arguments, workflow_cfg
from pipeline.Runner import logger, run_workflow


def main() -> None:
    args = read_arguments()
    CFG = workflow_cfg(args, logger)
    run_workflow(CFG)


if __name__ == "__main__":
    main()
