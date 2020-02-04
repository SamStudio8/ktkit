import sys
import os

from .ncbiTaxonomyTree import NcbiTaxonomyTree

def _load_dump(path):
    return NcbiTaxonomyTree(nodes_filename=os.path.join(path, "nodes.dmp"), names_filename=os.path.join(path, "names.dmp"))

def _get_tid_for_rank(tree, rank, tid):
    if tid == 0:
        return 0

    try:
        lineage = tree.getAscendantsWithRanksAndNames([tid])[tid]
    except KeyError:
        sys.stderr.write("[WARN] Taxon ID %d not found in NCBI Dump\n" % tid)
        return tid

    for node in lineage:
        if node.rank == rank:
            return node.taxid
    sys.stderr.write("Rank %s not found for taxon %d\n" % (rank, tid))
    return 1


def count(args):
    tree = _load_dump(args.dump)
    cache_map = {}

    mask = [_get_tid_for_rank(tree, args.rank, x) for x in args.mask]
    mask.extend(args.mask)
    mask.extend([0, 1])
    sys.stderr.write("[NOTE] Masking: %s\n" % str(mask))

    counts = {}
    total_bp = 0
    total_unmasked_bp = 0
    with open(args.input) as kh:
        for line in kh:
            fields = line.strip().split("\t")
            try:
                hit_tax = int(fields[2])
            except ValueError:
                # Get the taxid from the --use-names output
                hit_tax = int(fields[2].split("taxid ")[1][:-1])

            if hit_tax not in cache_map:
                cache_map[hit_tax] = _get_tid_for_rank(tree, args.rank, hit_tax)
            hit_tax = cache_map[hit_tax]

            hit_len = int(fields[3])

            if hit_tax not in counts:
                counts[hit_tax] = 0
            counts[hit_tax] += hit_len
            total_bp += hit_len

            if hit_tax not in mask:
                total_unmasked_bp += hit_len

    for tax_id in counts:

        if tax_id == 0:
            name = "unclassified"
        else:
            try:
                name = tree.getName([tax_id])[tax_id].replace(" ", "_")
            except KeyError:
                name = tax_id

        unmasked_total = 0
        if tax_id not in mask:
            unmasked_total = counts[tax_id]/total_unmasked_bp*100.0
        print("\t".join([str(x) for x in [
            tax_id,
            name,
            counts[tax_id],
            counts[tax_id]/total_bp*100.0,
            unmasked_total,
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
    parser.add_argument("--rank",help="Rank to output [default: species]", default="species")
    parser.add_argument("--dump", help="Path to NCBI taxonomy dump [default: ~/.ktkit/]", default="~/.ktkit")
    parser.add_argument("--mask", help="NCBI Taxon IDs to suppress [default: 9606]", default="9606")

    args = parser.parse_args()
    args.dump = os.path.expanduser(args.dump)
    args.mask = [int(x) for x in args.mask.split(",")]

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

