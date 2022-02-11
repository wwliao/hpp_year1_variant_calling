#!/usr/bin/env python3
import re
import gzip
import argparse
from os.path import basename
from collections import defaultdict

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

prefix = re.search("(\S+)\.vcf(?:\.gz)?", basename(args.vcffile))[1]
with gzip.open(args.vcffile, "rt", encoding="utf-8") as infile:
    with open(f"{prefix}.rmdup.vcf", "w") as outfile:
        seen_id = ""
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                cols = line.strip().split("\t")
                id = (cols[0], cols[1], cols[3], cols[4])
                info_fields = []
                for field in cols[7].split(";"):
                    if field.startswith("LV="):
                        lv = int(field.split("=")[1])
                    elif field.startswith("SS="):
                        ss = field.split("=")[1]
                    else:
                        info_fields.append(field)
                cols[7] = ";".join(info_fields)
                gt = cols[-1]
                if seen_id == "":
                    seen_id = id
                    record_cols = cols[:-2]
                    genotypes = defaultdict(list)
                    genotypes[(lv, ss)].append(gt)
                elif seen_id == id:
                    genotypes[(lv, ss)].append(gt)
                else:
                    outfile.write("\t".join(record_cols))
                    t = sorted(genotypes, reverse=True)[0]
                    lowest_lv, lowest_ss = t
                    outfile.write(f";LV={lowest_lv};SS={lowest_ss}\tGT")
                    if "0|1" in genotypes[t] and "1|0" in genotypes[t]:
                        genotypes[t].append("1|1")
                    genotype = sorted(genotypes[t], key=get_alt_count, reverse=True)[0]
                    outfile.write(f"\t{genotype}\n")
                    seen_id = id
                    record_cols = cols[:-2]
                    genotypes = defaultdict(list)
                    genotypes[(lv, ss)].append(gt)

        outfile.write("\t".join(record_cols))
        t = sorted(genotypes, reverse=True)[0]
        lowest_lv, lowest_ss = t
        outfile.write(f";LV={lowest_lv};SS={lowest_ss}\tGT")
        if "0|1" in genotypes[t] and "1|0" in genotypes[t]:
            genotypes[t].append("1|1")
        genotype = sorted(genotypes[t], key=get_alt_count, reverse=True)[0]
        outfile.write(f"\t{genotype}\n")
