# Graph variant decomposition (work in progress)

1. Extract a given sample from the VG-deconstructed VCF

    ```sh
    PREFIX=$(basename $VCF .vcf.gz)
    bcftools view -a -I -s $SAMPLE -Ou $VCF \
        | bcftools view -e 'GT="ref" | GT="0|." | GT=".|0" | GT=".|." | GT="." | GT="0/." | GT="./0" | GT="./."' \
                        -Oz -o $SAMPLE.$PREFIX.vcf.gz \
	    && bcftools index -t $SAMPLE.$PREFIX.vcf.gz
    ```

2. Drop any sites whose traversals aren’t nested in their parents’

	```sh
	drop_inconsistent_sites.py $SAMPLE.$PREFIX.vcf.gz
	```

3. Split multiallelic sites into biallelic records

    We need to compare traversals between reference and alternate to decompose graph variants.
	Therefore, it is easier to work with if the VCF file is biallelic instead of multiallelic.

	```sh
    bcftools norm -m -any \
	              -Oz -o $SAMPLE.$PREFIX.consistent.biallelic.vcf.gz \
				  $SAMPLE.$PREFIX.consistent.vcf.gz
	```

4. Convert the format of HPRC pangenome graphs from GFA to HashGraph/ODGI/PackedGraph

    There are two HPRC pangenome graphs:
    
    - Minigraph-Cactus: [hprc-v1.0-mc-grch38.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gfa.gz)
    - PGGB: [hprc-v1.0-pggb.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/hprc-v1.0-pggb.gfa.gz)

    Here, we chose to convert GFA to PackedGraph because of its low-memory footprint:

    ```sh
    vg convert -g -p -t 36 <input GFA> > <output PackedGraph>
    ```

5. Decompose graph variants through comparing traversals between reference and alternate

    `decompose_graph_variants.py` requires the following packages to be installed:

    - [`libbdsg`](https://github.com/vgteam/libbdsg): to access the HPRC pangenome graph in HashGraph/ODGI/PackedGraph format
    - [`cyvcf2`](https://github.com/brentp/cyvcf2): to parse and write the graph variants in VCF format
    

    ```sh
    python3 decompose_graph_variants.py -o $SAMPLE.$GRAPH.decomposed.unsorted.vcf.gz \
                                           $GRAPH.pg \
                                           $SAMPLE.$GRAPH.vcf.gz \
        && bcftools sort -m 10G \
                         -T $SAMPLE_sort_tmp/ \
                         -Oz -o $SAMPLE.$GRAPH.decomposed.vcf.gz \
                         $SAMPLE.$GRAPH.decomposed.unsorted.vcf.gz \
        && bcftools index -t $SAMPLE.$GRAPH.decomposed.vcf.gz \
        && rm $SAMPLE.$GRAPH.decomposed.unsorted.vcf.gz
    ```
6. Drop homozygous reference records

	Due to the redundant structures in the pangenome graph,
	there are records with GT != 0|0 or 0/0 but their traversals spell the same REF alleles (see below for 3 examples).

	```
	chr1 1249256 . ACTC  ACTC  60 . AT=>1780096>1780121>1780122>1780123,>1780096>1780097>1780098>1780123 GT 1|0
	chr1 7594118 . AGAAA AGAAA 60 . AT=>3716446>3716444>3716441,>3716446>3716443>3716442>3716441         GT 0|1
	chr1 8228919 . AAA   AAA   60 . AT=>3699499>3699498>3699488,>3699499<3699490<3699489>3699488         GT 1|1
	```

	Therefore, we need to drop these hidden homozygous reference records.

	```sh
	bcftools -e 'REF=ALT' \
	         -Oz -o $SAMPLE.$GRAPH.decomposed.concise.vcf.gz \
			 $SAMPLE.$GRAPH.decomposed.vcf.gz \
        && bcftools index -t $SAMPLE.$GRAPH.decomposed.concise.vcf.gz
	```

7. Update GTs and remove duplicates

    Records decomposed from different level of snarls might represent exactly the same variant but with different genotypes.

    ```sh
    python3 remove_duplicates.py $SAMPLE.$GRAPH.decomposed.concise.vcf.gz
    ```

# To-do
- Figure out why there are still redundant structures in Minigraph-Cactus

# Known issues
- VG-deconstructed AT fields not necessarily consistent with different levels (see [issue #3541](https://github.com/vgteam/vg/issues/3541))
