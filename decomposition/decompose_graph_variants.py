#!/usr/bin/env python3
import re
import argparse
#from bdsg.bdsg import ODGI
from bdsg.bdsg import HashGraph
from cyvcf2 import VCF, Writer

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", required=True)
#parser.add_argument("graph", help="graph in ODGI format")
parser.add_argument("graph", help="graph in HashGraph format")
parser.add_argument("vcffile", help="input VCF file")
args = parser.parse_args()

def get_steps(path):
    return re.findall("[><]\d+", path)

def reverse_path(path):
    d = str.maketrans("><", "<>")
    rev_path = []
    for step in get_steps(path)[::-1]:
        rev_path.append(step.translate(d))
    return "".join(rev_path)

def decompose_traversal(query, target):
    query_steps = get_steps(query)
    target_steps = get_steps(target)
    assert(target_steps[0] == query_steps[0])
    assert(target_steps[-1] == query_steps[-1])
    source = (0, 0)
    for i, target_step in enumerate(target_steps):
        if i > 0 and target_step in query_steps[source[1]+1:]:
            j = source[1] + 1 + query_steps[source[1]+1:].index(target_step)
            target_interval = i - source[0] - 1
            query_interval = j - source[1] - 1
            if target_interval > 0 or query_interval > 0:
                sink = (i, j)
                target_path = "".join(target_steps[source[0]:sink[0]+1])
                query_decomposed_steps = query_steps[source[1]:sink[1]+1]
                query_path = "".join(query_decomposed_steps)
                rev_target_steps = get_steps(reverse_path(target_path))
                source = (i, j)
                has_inv = False
                if set(rev_target_steps) & set(query_decomposed_steps):
                    has_inv = True
                yield query_path, target_path, has_inv
            else:
                source = (i, j)

def get_allele_seq(path, graph):
    seq = ""
    # Exclude the first and last steps because they are the same as REF
    for step in get_steps(path)[1:-1]:
        strand = step[0]
        node = int(step[1:])
        if strand == ">":
            seq += graph.get_sequence(graph.get_handle(node, False))
        else:
            seq += graph.get_sequence(graph.get_handle(node, True))
    return seq

"""
def get_allele_seq(query, target, graph):
    target_steps = get_steps(target)
    query_steps = get_steps(query)
    assert(target_steps[0] == query_steps[0])
    assert(target_steps[-1] == query_steps[-1])
    target_seq = ""
    for target_step in target_steps[1:-1]:
        target_node = int(target_step[1:])
        if target_step[0] == ">":
            target_seq += graph.get_sequence(graph.get_handle(target_node, False))
        else:
            target_seq += graph.get_sequence(graph.get_handle(target_node, True))
    query_seq = ""
    for query_step in query_steps[1:-1]:
        query_node = int(query_step[1:])
        if query_step[0] == ">":
            query_seq += graph.get_sequence(graph.get_handle(query_node, False))
        else:
            query_seq += graph.get_sequence(graph.get_handle(query_node, True))
    return query_seq, target_seq
"""

#graph = ODGI()
graph = HashGraph()
graph.deserialize(args.graph)

vcf = VCF(args.vcffile)
#w = Writer(args.output, vcf)
with open(args.output, "w") as outfile:
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        id = variant.ID
        gt = set(variant.genotypes[0][:2])
        if (1 in gt) or (2 in gt):
            ats = variant.INFO.get("AT").split(",")
            ref_at = ats[0]
            alt_ats = []
            for i in sorted(gt):
                if i != 0:
                    alt_ats.append(ats[i])
            for alt_at in alt_ats:
                for i, (alt_path, ref_path, has_inv) in enumerate(decompose_traversal(alt_at, ref_at)):
                    if has_inv:
                        outfile.write(f"{ref_path}\t{alt_path}\t1\n")
                    else:
                        outfile.write(f"{ref_path}\t{alt_path}\t0\n")
                    #alt_seq, ref_seq = get_allele_seq(alt_path, ref_path, graph)
                    ref_seq = get_allele_seq(ref_path, graph)
                    alt_seq = get_allele_seq(alt_path, graph)
                    if alt_seq == ref_seq:
                        print("There is a redundant call!")
                    else:
                        if len(alt_seq) == len(ref_seq):
                            print("This is a record without preceding base.")
                        else:
                            print("This is a record with preceding base.")
                            
                    if has_inv:
                        print(f"{id}-{i+1}: {ref_path}, {alt_path} (INV)")
                    else:
                        print(f"{id}-{i+1}: {ref_path}, {alt_path}")

    vcf.close()
    #w.close()
