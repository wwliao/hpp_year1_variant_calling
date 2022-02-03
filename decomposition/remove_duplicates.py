#!/usr/bin/env python3
import re
import gzip
import argparse
from os.path import basename

parser = argparse.ArgumentParser()
parser.add_argument("vcffile", help="input sorted VCF file")
args = parser.parse_args()

def get_alt_count(genotype):
    count = 0
    for allele in genotype.split("|"):
        if allele == ".":
            count -= 1
        else:
            count += int(allele)
    return count

prefix = re.search("(\S+).vcf.gz", basename(args.vcffile))[1]
with gzip.open(args.vcffile, "rt", encoding="utf-8") as infile:
    with gzip.open(f"{prefix}.rmdup.vcf.gz", "wt") as outfile:
        seen_id = ""
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                cols = line.strip().split("\t")
                id = (cols[0], cols[1], cols[3], cols[4])
                gt = cols[-1]
                if seen_id == "":
                    seen_id = id
                    record_cols = None
                    if "RS=" in cols[7]:
                        record_cols = cols[:-1]
                    genotypes = [gt]
                elif seen_id == id:
                    if "RS=" in cols[7]:
                        record_cols = cols[:-1]
                    genotypes.append(gt)
                else:

                    # There are some nested traversals don't exist in root traversals
                    # Skip those record for now
                    if record_cols is not None:
                        outfile.write("\t".join(record_cols))
                        genotype = sorted(genotypes, key=get_alt_count, reverse=True)[0]
                        outfile.write(f"\t{genotype}\n")
                    seen_id = id
                    record_cols = None
                    if "RS=" in cols[7]:
                        record_cols = cols[:-1]
                    genotypes = [gt]

        # There are some nested traversals don't exist in root traversals
        # Skip those record for now
        if record_cols is not None:
            outfile.write("\t".join(record_cols))
            genotype = sorted(genotypes, key=get_alt_count, reverse=True)[0]
            outfile.write(f"\t{genotype}\n")
