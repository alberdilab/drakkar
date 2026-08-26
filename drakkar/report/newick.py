"""A small Newick reader for the report's circular phylogeny.

Drakkar already writes pruned GTDB-Tk trees — ``annotating/bacteria.tree`` and,
when the catalogue holds archaea, ``annotating/archaea.tree`` — so the report
only has to read them. Biopython could, but it is a workflow dependency rather
than a reporting one, and the report modules deliberately import nothing beyond
pandas: this file is the whole of what the renderer needs.

The parser is written against what GTDB-Tk actually emits. Internal nodes carry
support values and lineage strings, which Biopython writes back out as
single-quoted labels holding colons, semicolons and doubled quotes, so labels
are unquoted here rather than split on punctuation. Comments in ``[...]`` are
dropped, and branch lengths may be written in scientific notation.
"""

import re

# A branch length as Newick allows it, scientific notation included.
_LENGTH = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

# Everything that ends an unquoted label or a branch length.
_PUNCTUATION = set("(),:;[]")


class NewickError(ValueError):
    """Raised when a string is not a tree this reader can make sense of."""


class Node:
    """One clade: a name, the length of the branch leading to it, children."""

    __slots__ = ("name", "length", "children")

    def __init__(self, name="", length=0.0, children=None):
        self.name = name
        self.length = length
        self.children = children if children is not None else []

    @property
    def is_leaf(self):
        return not self.children

    def leaves(self):
        """Every tip under this node, left to right."""
        found = []
        for node in self.walk():
            if not node.children:
                found.append(node)
        return found

    def walk(self):
        """This node and every node under it, parents before children.

        Iterative, like everything else here: a tree pruned down to a ladder
        is as deep as it is wide, and a recursive walk of one would exhaust the
        interpreter's stack on a catalogue of a few hundred genomes.
        """
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))


def parse(text):
    """Parse one Newick string into its root :class:`Node`.

    Raises :class:`NewickError` on anything unparsable, which the caller treats
    as "this output directory has no usable tree" rather than as a failure.
    """
    if text is None:
        raise NewickError("no tree text")
    source = text.strip()
    if not source:
        raise NewickError("empty tree text")
    reader = _Reader(source)
    root = reader.tree()
    reader.skip_space()
    if reader.peek() == ";":
        reader.advance()
    return root


def count_tips(text):
    """How many tips a Newick string holds, or 0 when it is not readable."""
    try:
        return len(parse(text).leaves())
    except NewickError:
        return 0


class _Reader:
    """A cursor over the Newick string, one clade at a time."""

    def __init__(self, source):
        self.source = source
        self.position = 0

    def peek(self):
        return self.source[self.position] if self.position < len(self.source) else ""

    def advance(self):
        self.position += 1

    def skip_space(self):
        while self.position < len(self.source) and self.source[self.position].isspace():
            self.position += 1

    def skip_comment(self):
        """Drop a ``[...]`` comment, which carries no phylogeny."""
        while self.peek() == "[":
            depth = 0
            while self.position < len(self.source):
                character = self.source[self.position]
                self.advance()
                if character == "[":
                    depth += 1
                elif character == "]":
                    depth -= 1
                    if depth == 0:
                        break
            self.skip_space()

    def tree(self):
        """Read one clade and everything under it, without recursing.

        ``open`` holds the sibling lists of the clades still being read, one
        entry per unclosed ``(``, so a tree nested a thousand deep costs a
        thousand list entries rather than a thousand stack frames.
        """
        open_clades = []
        while True:
            self.skip_space()
            self.skip_comment()
            while self.peek() == "(":
                self.advance()
                open_clades.append([])
                self.skip_space()
                self.skip_comment()
            node = Node(*self.described())
            while True:
                self.skip_space()
                character = self.peek()
                if character == ",":
                    self.advance()
                    if not open_clades:
                        raise NewickError(f"stray comma at character {self.position}")
                    open_clades[-1].append(node)
                    break
                if character == ")":
                    self.advance()
                    if not open_clades:
                        raise NewickError(f"stray ')' at character {self.position}")
                    children = open_clades.pop()
                    children.append(node)
                    self.skip_space()
                    self.skip_comment()
                    node = Node(*self.described(), children=children)
                    continue
                if open_clades:
                    raise NewickError(f"unclosed clade at character {self.position}")
                return node

    def described(self):
        """The name and branch length written after a clade or beside a tip."""
        name = self.label()
        self.skip_space()
        self.skip_comment()
        length = self.length()
        self.skip_comment()
        return name, length

    def label(self):
        """A node name: single-quoted with ``''`` escapes, or bare."""
        if self.peek() == "'":
            self.advance()
            parts = []
            while self.position < len(self.source):
                character = self.source[self.position]
                self.advance()
                if character == "'":
                    if self.peek() == "'":
                        self.advance()
                        parts.append("'")
                        continue
                    return "".join(parts)
                parts.append(character)
            raise NewickError("unterminated quoted label")
        start = self.position
        while self.position < len(self.source):
            character = self.source[self.position]
            if character in _PUNCTUATION or character.isspace():
                break
            self.advance()
        # Underscores are left alone: the standard reads them as spaces, but
        # here the tips are MAG names, and those carry real underscores that
        # the taxonomy table is joined on.
        return self.source[start:self.position].strip()

    def length(self):
        # Trees are wrapped across lines by some writers, so the colon need not
        # sit against the label it belongs to.
        self.skip_space()
        if self.peek() != ":":
            return 0.0
        self.advance()
        self.skip_space()
        match = _LENGTH.match(self.source, self.position)
        if match is None:
            # A colon with no number after it: no length, not a broken tree.
            return 0.0
        self.position = match.end()
        return float(match.group())
