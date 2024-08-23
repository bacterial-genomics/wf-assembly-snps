# wf-assembly-snps: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.3 - August 23, 2024

### `Added`

### `Fixed`

- [#6](https://github.com/bacterial-genomics/wf-assembly-snps/pull/6) [commit](https://github.com/bacterial-genomics/wf-assembly-snps/commit/d741e2bdfd7158ff6f53009f47ee43572abf2d6d) Decreased input filesize requirements for test and test_full configs by @chrisgulvik
- [#6](https://github.com/bacterial-genomics/wf-assembly-snps/pull/6) [commit](https://github.com/bacterial-genomics/wf-assembly-snps/commit/8b280aebe36706436dc602ef3f09be6f8a02a0ff) Revised ".ref" replacement in Parsnp.tree and Parsnp.SNPs.fa outfiles, fixes #1 issue. This impacts users with ".ref" in >1 input supplied file to now enable recombination (e.g., Gubbins) to succeed. by @chrisgulvik

### `Updated`

### `Deprecated`

- [#6](https://github.com/bacterial-genomics/wf-assembly-snps/pull/6) [commit](https://github.com/bacterial-genomics/wf-assembly-snps/commit/44e0f25bd71bf2c1191bbcb133fc917ccf039281) No longer create VCF outfile (easily made from `harvesttools -i parsnp.ggr -V out.vcf`) but we rarely need this and slows down larger panel analyses. by @chrisgulvik

## v1.0.2 - August 21, 2024

### `Added`

### `Fixed`

- [#3](https://github.com/bacterial-genomics/wf-assembly-snps/pull/3) Refactored profiles to align with Seqera Tower requirements and added a Tower yaml for reporting by @slsevilla

### `Updated`

### `Deprecated`

## v1.0.1 - August 21, 2024

### `Added`

### `Fixed`

- [#2](https://github.com/bacterial-genomics/wf-assembly-snps/pull/2) Lint fixes (@gregorysprenger)

### `Updated`

- [#2](https://github.com/bacterial-genomics/wf-assembly-snps/pull/2) Updated docs (@gregorysprenger)

### `Deprecated`

## v1.0.0 - April 1, 2024

Initial release of wf-assembly-snps.
