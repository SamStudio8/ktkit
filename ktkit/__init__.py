import sys
import os

from .ncbiTaxonomyTree import NcbiTaxonomyTree

def _load_dump(path):
    tree = NcbiTaxonomyTree(nodes_filename=os.path.join(path, "nodes.dmp"), names_filename=os.path.join(path, "names.dmp"))

def count(args):
    tree = _load_dump(args.dump)

    counts = {}
    total_bp = 0
    with open(args.input) as kh:
        for line in kh:
            fields = line.strip().split("\t")
            try:
                hit_tax = int(fields[2])
            except ValueError:
                # Get the taxid from the --use-names output
                hit_tax = int(fields[2].split("taxid ")[1][:-1])
            hit_len = int(fields[3])

            if hit_tax not in counts:
                counts[hit_tax] = 0
            counts[hit_tax] += hit_len
            total_bp += hit_len

    for tax_id in counts:
        print("\t".join([str(x) for x in [
            tax_id,
            counts[tax_id],
            counts[tax_id]/total_bp*100.0
        ]]))

def check(args):
    if os.path.exists(os.path.join(args.dump, "nodes.dmp")) and os.path.exists(os.path.join(args.dump, "names.dmp")):
        return True
    return False

def cli():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mode")
    parser.add_argument("input")
    parser.add_argument("--dump", help="Path to NCBI taxonomy dump [default: ~/.ktkit/]", default="~/.ktkit")

    args = parser.parse_args()
    args.dump = os.path.expanduser(args.dump)

    if not check(args):
        sys.stderr.write("NCBI dump not found in %s\n" % args.dump)
        sys.stderr.write("mkdir -p %s; cd %s; wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz; tar xvf taxdump.tar.gz\n" % (args.dump, args.dump))
        return

    modes = {
        "count": count,
    }
    if args.mode in modes:
        modes[args.mode](args)
    else:
        print("Command not found. Meow.")

