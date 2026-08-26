"""The Newick reader the circular phylogeny is laid out from.

The trees it has to cope with are the ones GTDB-Tk writes and the pruning rule
rewrites through Biopython: quoted internal labels holding the punctuation
Newick otherwise reserves, branch lengths in scientific notation, and — for a
placement made on topology alone — no branch lengths at all.
"""

import unittest

from drakkar.report.newick import NewickError, count_tips, parse


class ParseTests(unittest.TestCase):
    def test_tips_are_returned_in_the_order_the_file_lists_them(self):
        tree = parse("((A:0.1,B:0.2):0.3,C:0.4);")
        self.assertEqual([leaf.name for leaf in tree.leaves()], ["A", "B", "C"])

    def test_branch_lengths_are_read_onto_the_node_they_lead_to(self):
        tree = parse("(A:0.125,B:0.25):0.5;")
        self.assertEqual(tree.length, 0.5)
        self.assertEqual([leaf.length for leaf in tree.leaves()], [0.125, 0.25])

    def test_scientific_notation_is_a_branch_length(self):
        tree = parse("(A:1.5e-06,B:2E3);")
        self.assertEqual([leaf.length for leaf in tree.leaves()], [1.5e-06, 2000.0])

    def test_a_quoted_label_keeps_its_punctuation(self):
        # What Biopython writes back out for a GTDB-Tk internal node: a support
        # value and a lineage, colons and semicolons and all.
        tree = parse("((A,B)'100.0:d__Bacteria; p__Bacillota',C);")
        self.assertEqual(tree.children[0].name, "100.0:d__Bacteria; p__Bacillota")

    def test_a_doubled_quote_inside_a_quoted_label_is_one_quote(self):
        tree = parse("(A,'it''s')  ;")
        self.assertEqual([leaf.name for leaf in tree.leaves()], ["A", "it's"])

    def test_underscores_in_tip_names_are_left_alone(self):
        # The standard reads them as spaces; here they are MAG names, and the
        # taxonomy table is joined on them.
        self.assertEqual(
            [leaf.name for leaf in parse("(MAG_1:0.1,MAG_2:0.2);").leaves()],
            ["MAG_1", "MAG_2"],
        )

    def test_comments_carry_no_phylogeny_and_are_dropped(self):
        tree = parse("(A[&rate=0.1]:0.1,B:0.2)[comment];")
        self.assertEqual([leaf.name for leaf in tree.leaves()], ["A", "B"])
        self.assertEqual(tree.leaves()[0].length, 0.1)

    def test_a_tree_without_branch_lengths_still_parses(self):
        tree = parse("((A,B),(C,D));")
        self.assertEqual(len(tree.leaves()), 4)
        self.assertEqual([leaf.length for leaf in tree.leaves()], [0.0] * 4)

    def test_whitespace_and_a_missing_semicolon_are_tolerated(self):
        tree = parse("  ( A : 0.1 , B : 0.2 )  ")
        self.assertEqual([leaf.name for leaf in tree.leaves()], ["A", "B"])

    def test_walk_visits_parents_before_children(self):
        names = [node.name for node in parse("((A,B)ab,C)root;").walk()]
        self.assertEqual(names, ["root", "ab", "A", "B", "C"])

    def test_a_single_tip_is_a_tree_of_one_leaf(self):
        tree = parse("A:0.1;")
        self.assertEqual([leaf.name for leaf in tree.leaves()], ["A"])

    def test_empty_text_is_not_a_tree(self):
        for value in (None, "", "   "):
            with self.assertRaises(NewickError):
                parse(value)

    def test_an_unclosed_clade_is_not_a_tree(self):
        with self.assertRaises(NewickError):
            parse("((A,B,C);")

    def test_an_unterminated_quote_is_not_a_tree(self):
        with self.assertRaises(NewickError):
            parse("(A,'B);")


class CountTipsTests(unittest.TestCase):
    def test_tips_are_counted_without_the_caller_parsing(self):
        self.assertEqual(count_tips("((A,B),(C,D));"), 4)

    def test_an_unreadable_tree_counts_zero_rather_than_raising(self):
        # The ingest asks for a count to decide whether a file is worth
        # storing, so a broken file has to answer rather than fail.
        self.assertEqual(count_tips("((A,B);"), 0)
        self.assertEqual(count_tips(""), 0)


if __name__ == "__main__":
    unittest.main()
