#! /usr/bin/env python3
#
# Copyright 2020 Joshua Watt <JPEWhacker@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import codecs
import collections
import contextlib
import csv
import glob
import itertools
import os
import sys
import unittest

Gene = collections.namedtuple("Gene", ("id", "iso", "foldcount", "p_value"))


def sanitize_id(ident):
    id_fields = ident.split("-")
    if id_fields[0] == "tRNA":
        id_fields = id_fields[1:]
    return "-".join(id_fields[:3])


def load_file(path):
    """
    Parses a CSV file of genes
    """
    genes = dict()
    skipped = 0

    with codecs.open(path, "r", "iso-8859-1") as f:
        reader = csv.reader(f)
        # Skip header
        next(reader)
        for l in reader:
            try:
                fields = l[0 : len(Gene._fields)]
                fields[0] = sanitize_id(fields[0])

                g = Gene(*fields)
                genes[g.id] = g
            except TypeError:
                skipped += 1

    return genes, skipped


class Dataset(object):
    def __init__(self, genes, name):
        self.genes = genes
        self.name = name
        self.ids = frozenset(genes.keys())

    def __repr__(self):
        return "Dataset %s" % self.name

    def __str__(self):
        return self.name


class DataSlice(object):
    def __init__(self, keep, exclude, collection):
        self.keep = keep
        self.exclude = exclude
        self.name = "(%s) - (%s)" % (
            " & ".join(d.name for d in keep),
            " | ".join(d.name for d in exclude),
        )
        self.ids = self._calc_ids(keep, exclude)
        self.collection = collection

    def __str__(self):
        return self.name

    def _calc_ids(self, keep_datasets, exclude_datasets):
        keep_ids = [d.ids for d in keep_datasets]
        exclude_ids = [d.ids for d in exclude_datasets]

        return sorted(
            list(keep_ids[0].intersection(*keep_ids[1:]).difference(*exclude_ids))
        )

    def fields(self):
        fields = list(Gene._fields[0:2])
        for f in Gene._fields[2:]:
            for d in sorted(self.keep, key=lambda d: d.name):
                fields.append("%s %s" % (d.name, f))
        return tuple(fields)

    def genes(self):
        for i in self.ids:
            name = self.collection.gene_names.get(i, "")
            values = [i, name]
            for f in Gene._fields[2:]:
                for d in sorted(self.keep, key=lambda d: d.name):
                    values.append(getattr(d.genes[i], f))
            yield tuple(values)


class DatasetCollection(object):
    def __init__(self):
        self.datasets = set()
        self.gene_names = {}

    def add(self, dataset):
        self.datasets.add(dataset)

    def slices(self):
        for l in range(len(self.datasets)):
            combos = itertools.combinations(self.datasets, l + 1)
            for c in combos:
                keep = set(c)
                exclude = self.datasets - keep
                yield DataSlice(keep, exclude, self)

    def __iter__(self):
        return iter(self.datasets)

    def __len__(self):
        return len(self.datasets)


class TestCalc(unittest.TestCase):
    def test_calc(self):
        A, _ = load_file(os.path.join("test", "2", "A.csv"))
        B, _ = load_file(os.path.join("test", "2", "B.csv"))

        collection = DatasetCollection()
        collection.add(Dataset(A, "A"))
        collection.add(Dataset(B, "B"))

        values = set(tuple(s.genes()) for s in collection.slices())

        self.assertEqual(
            values,
            {
                (
                    (
                        "Alpha-ABC-1",
                        "",
                        "B-Alpha-fc",
                        "B-Alpha-p",
                    ),
                ),
                (
                    (
                        "Delta-ABC-1",
                        "",
                        "A-Delta-fc",
                        "A-Delta-p",
                    ),
                ),
                (
                    (
                        "Beta-DEF-2",
                        "",
                        "A-Beta-fc",
                        "B-Beta-fc",
                        "A-Beta-p",
                        "B-Beta-p",
                    ),
                ),
            },
        )


def write_slice(out, s, gene_names):
    out.write("%s,Count,%d\n" % (s.name, len(s.ids)))

    out.write("%s\n" % ",".join(s.fields()))
    for g in s.genes():
        out.write("%s\n" % ",".join(g))
    out.write("\n")


def main():
    @contextlib.contextmanager
    def get_output(name):
        if name == "-":
            yield sys.stdout
        else:
            with open(name, "w") as f:
                yield f

    parser = argparse.ArgumentParser(description="Gene comparison")
    parser.add_argument("output", help="Output CSV file ('-' for stdout)")
    parser.add_argument(
        "--files", nargs="*", default=[], help="Add files based on a glob"
    )

    args = parser.parse_args()

    collection = DatasetCollection()
    files = list()

    for g in args.files:
        files.extend((p, os.path.basename(p)) for p in glob.iglob(g))

    for (path, name) in files:
        genes, skipped = load_file(path)
        d = Dataset(genes, name)
        collection.add(d)
        print("Loaded %d genes from %s (skipped %d)" % (len(d.genes), path, skipped))

    print("----------------------")
    print(
        "Loaded %d total genes from %d files"
        % (sum(len(d.genes) for d in collection), len(collection))
    )

    with get_output(args.output) as f:
        for s in collection.slices():
            write_slice(f, s, collection.gene_names)

    return 0


if __name__ == "__main__":
    sys.exit(main())
