# Graph variant decomposition (work in progress)

1. Extract a given sample from the VG-deconstructed VCF

    ```sh
    PREFIX=$(basename $VCF .vcf.gz)
    bcftools view -a -I -s $SAMPLE -Ou $VCF \
        | bcftools view -e 'GT="ref" | GT="0|." | GT=".|0" | GT=".|." | GT="."' -Ou \
        | bcftools norm -m - -Oz -o $SAMPLE.$PREFIX.vcf.gz \
	    && bcftools index -t $SAMPLE.$PREFIX.vcf.gz
    ```

2. Convert the format of HPRC pangenome graphs from GFA to HashGraph

    There are two HPRC pangenome graphs:
    
    - Minigraph-Cactus: [hprc-v1.0-mc-grch38.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gfa.gz)
    - PGGB: [hprc-v1.0-pggb.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/hprc-v1.0-pggb.gfa.gz)

    Using VG to convert GFA to HashGraph

    ```sh
    vg convert -g -a -t 36 <input GFA> > <output HashGraph>
    ```

3. Decompose graph variants

    `decompose_graph_variants.py` requires the following packages to be installed:

    - [`libbdsg`](https://github.com/vgteam/libbdsg): to access the HPRC pangenome graph in HashGraph format
    - [`cyvcf2`](https://github.com/brentp/cyvcf2): to parse and write the graph variants in VCF format
    

    ```sh
    python3 decompose_graph_variants.py -o $SAMPLE.$GRAPH.decomposed.vcf.gz \
                                           $GRAPH.hg \
                                           $SAMPLE.$GRAPH.vcf.gz \
        && bcftools sort -m 10G \
                         -T $SAMPLE_sort_tmp/ \
                         -Oz -o $SAMPLE.$GRAPH.decomposed.sorted.vcf.gz \
                         $SAMPLE.$GRAPH.decomposed.vcf.gz \
        && bcftools index -t $SAMPLE.$GRAPH.decomposed.sorted.vcf.gz
    ```

4. Remove duplicates

    Records decomposed from different level of snarls might represent exactly the same variant but different genotypes.
    We only keep those derived from the root snarls (i.e. records with INFO/RS) with a genotype that best fit the corresponding record.

    ```sh
    python3 remove_duplicates.py $SAMPLE.$GRAPH.decomposed.sorted.vcf.gz
    ```

# TODO

- Switch from HashGraph to PackedGraph to reduce the memory usage (PGGB in HashGraph: >500GB memory)
- Remove duplicated records and update genotypes
- Figure out why there are still redundant structures in Minigraph-Cactus
- Some lower-level allele traversals don't exist in higher-level allele traversals. For example, check the top-level snarl `>2539380>2539432` and its nested snarls.
