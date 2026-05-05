from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

from Bio import Phylo


def load_prune_module():
    module_path = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "scripts" / "prune_gtdbtk_tree.py"
    spec = importlib.util.spec_from_file_location("prune_gtdbtk_tree", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class GtdbtkTreePruningTests(unittest.TestCase):
    def test_prune_tree_keeps_only_query_genomes(self) -> None:
        module = load_prune_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batchfile = tmpdir_path / "mag_input.tsv"
            batchfile.write_text(
                "/tmp/MAG_A.fa\tMAG_A\n"
                "/tmp/MAG_B.fa\tMAG_B\n",
                encoding="utf-8",
            )

            input_tree = tmpdir_path / "gtdbtk.backbone.bac120.classify.tree"
            input_tree.write_text(
                "(GB_GCA_000001:0.1,(MAG_A:0.2,GB_GCA_000002:0.3):0.4,MAG_B:0.5);\n",
                encoding="utf-8",
            )

            output_tree = tmpdir_path / "bacteria.tree"
            module.prune_tree(str(input_tree), str(batchfile), str(output_tree))

            pruned_tree = Phylo.read(str(output_tree), "newick")
            terminal_names = {terminal.name for terminal in pruned_tree.get_terminals()}
            self.assertEqual(terminal_names, {"MAG_A", "MAG_B"})

    def test_annotating_workflow_targets_bacteria_tree(self) -> None:
        snakefile = (Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "Snakefile").read_text(encoding="utf-8")
        taxonomy_rules = (
            Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "rules" / "annotating_taxonomy.smk"
        ).read_text(encoding="utf-8")

        self.assertIn('f"{OUTPUT_DIR}/annotating/bacteria.tree"', snakefile)
        self.assertIn("rule gtdbtk_pruned_trees:", taxonomy_rules)
        self.assertIn('archaea_output=f"{OUTPUT_DIR}/annotating/archaea.tree"', taxonomy_rules)


if __name__ == "__main__":
    unittest.main()
