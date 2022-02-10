#!/usr/bin/env python3
import re
import argparse
from cyvcf2 import VCF, Writer

parser = argparse.ArgumentParser()
parser.add_argument("vcffile")
args = parser.parse_args()

def get_source_snarl(snarls, id):
    ps_id = snarls[id]["parent"]
    while ps_id in snarls:
        id = ps_id
        ps_id = snarls[ps_id]["parent"]
    return id

# First round to get the source snarls for all nested snarls
snarls = {}
not_top_snarl_ids = set()
vcf = VCF(args.vcffile)
for v in vcf:
    id = v.ID
    level = v.INFO.get("LV")
    parent = "*"
    if level > 0:
        parent = v.INFO.get("PS")
        not_top_snarl_ids.add(id)
    if id not in snarls:
        snarls[id] = {"level": level, "parent": parent}
vcf.close()

map2source = {}
for id in sorted(not_top_snarl_ids):
    source_id = get_source_snarl(snarls, id)
    if source_id != id:
        map2source[id] = source_id

# Second round to get the source ALT traversals
vcf = VCF(args.vcffile)
source_alts = {}
for v in vcf:
    id = v.ID
    source_ids = set(map2source.values())
    if id in source_ids:
        # Drop the REF AT, which is always the first one in AT field
        alts = v.INFO.get("AT").split(",")[1:]
        if id not in source_alts:
            source_alts[id] = set(alts)
        else:
            source_alts[id].update(alts)
vcf.close()

# Final round to check if child traversal is nested in parent's 
vcf = VCF(args.vcffile)
prefix = re.search("(\S+)\.vcf\.gz", args.vcffile)[1]
w = Writer(f"{prefix}.consistent.vcf.gz", vcf)
e = Writer(f"{prefix}.inconsistent.vcf.gz", vcf)
for v in vcf:
    id = v.ID
    if id in map2source:
        source_id = map2source[id]
        nested_count = 0
        alts = v.INFO.get("AT").split(",")[1:]
        for alt in alts:
            for source_alt in source_alts[source_id]:
                if re.search(alt + "\D+", source_alt + "$"):
                    nested_count += 1
                    break
        if nested_count == len(alts):
            w.write_record(v)
        else:
            e.write_record(v)
    else:
        w.write_record(v)
vcf.close()
w.close()
e.close()
