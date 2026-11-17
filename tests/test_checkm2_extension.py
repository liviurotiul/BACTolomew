import unittest
from pathlib import Path
from unittest import mock

from pipeline import AssemblyQualityWorkflow as aqw


class CheckM2ExtensionTest(unittest.TestCase):
    """CheckM2 alone scans only .fna: genomes are staged under one extension passed with -x."""

    def _cmd(self, extra_args):
        ext = aqw._checkm2_extension(extra_args)
        with mock.patch.object(aqw, "tool_cmd", return_value=["checkm2"]):
            return ext, aqw._checkm2_command(
                input_dir=Path("in"), output_dir=Path("out"), threads=2,
                settings={}, extra_args=extra_args, extension=ext,
            )

    def test_default_and_user_extension(self):
        ext, cmd = self._cmd([])
        self.assertEqual(ext, "fna")
        self.assertEqual(cmd[-2:], ["-x", ".fna"])
        # The user's own -x still decides the staged extension and is passed once.
        ext, cmd = self._cmd(["-x", ".fasta", "--remove_intermediates"])
        self.assertEqual(ext, "fasta")
        self.assertEqual(cmd.count("-x"), 1)
        self.assertEqual(cmd[-3:], ["-x", ".fasta", "--remove_intermediates"])
        self.assertEqual(aqw._checkm2_extension(["--extension=fa"]), "fa")


if __name__ == "__main__":
    unittest.main()
