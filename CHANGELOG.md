# Changelog

## [0.3.0](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.2.5...v0.3.0) (2025-11-19)


### Features

* downsample bam file to bamsnap ([1c89f9c](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/1c89f9c566c43518faa0bd4e7b69c68c1234dc7b))
* downsampling deduplicated bam files that are input to bamsnap ([f06b1fc](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/f06b1fca496887c5f62770d737588c94daebfe18))
* downsampling merged deduplicated bam files that are input to bamsnap ([24ee530](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/24ee5307ce81598f756765bfaa76a5c6b7b636b1))


### Bug Fixes

* cleaned unwanted rebased sections ([5c04ced](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/5c04ced1a09b89ec139f75c639cc1796fb9494fe))
* container for bamsnap ([7245c3f](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/7245c3f98f017464091b452dbfcc3e0f41cbcc5c))
* input and output file extensions in bamsnap rule ([67410b9](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/67410b9e6a2ac31a01837ec597d268c89c0cf5e8))
* linting and poppy uu version ([e1ce1ba](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e1ce1baa06dd866d54d46d1b18aeedfe390089b6))
* linting and poppy uu version ([c7ab302](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/c7ab302e2f2911dc26ac5eaf51abf69484bee469))
* linting and poppy uu version ([b345193](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/b3451934d8136f76b16f6f894a4f40868aa8e458))
* missing white space ([51901be](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/51901be447be9425641278aaef5cce22d2645c0b))
* poppies versions in snakemake-dry-run.yaml ([de47d4a](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/de47d4aa8a403d6db41a9496a21012784630e00f))
* Reduce max_reads in bamsnap_downsample_bam config ([10f4f1f](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/10f4f1f39a6b136ab2cef23249ec9af7d81e826d))
* Reduce max_reads in bamsnap_downsample_bam config ([5203c18](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/5203c1834802af1e72390307515093a626702f3e))
* remove echo for poppy_uppsala version in version.smk ([2f46f27](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/2f46f277f014cb551c1d32e31b127eedd8e9b967))
* smkfmt ([585bc1e](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/585bc1ed83d3d0974d0d9d448c26ad4c605dca21))
* smkfmt ([1e4425a](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/1e4425a1d92d25816091ac36f9280f5a2ede0329))
* use downsampled files in bamsnap ([6cf9cdf](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/6cf9cdfaf2e9bdb524eabf8d4327db79eefbaa2f))

## [0.2.5](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.2.4...v0.2.5) (2025-10-09)


### Bug Fixes

* padded region in bed file ([786e731](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/786e73137f65e5094d2f257953e427804dd8b1d4))

## [0.2.4](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.2.3...v0.2.4) (2025-09-10)


### Bug Fixes

* adjust ruleorder ([63285b2](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/63285b232185d7f5d83ea3227ef61674c88ba006))
* localrule and index file ([22c1f8d](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/22c1f8d3fa0e64b9a6d4abb6dceebbd63c44ca28))
* ruleorder ([b49fdbd](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/b49fdbdb22970af7737afa12b48cf1ee528822a2))

## [0.2.1](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.2.0...v0.2.1) (2025-06-03)


### Bug Fixes

* hotspot genes in excel report ([e57ce10](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e57ce107fe4a0c8e816a8eb871a0a962046e5691))
* removed pureCN ([3bc3a16](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/3bc3a166265c5f9ba96462c7bfb22bcdc7f4903c))

## [0.2.0](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.1.2...v0.2.0) (2025-04-25)


### Features

* container for venv with bamsnap ([cb824d5](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/cb824d5ae6fb1e657b1790b554760cfaaf1c70be))


### Bug Fixes

* linting container directive ([351504c](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/351504c0761bc746ac236fb0eeb6bf29f82d38c4))
* linting snakefmt ([3f33f43](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/3f33f43052a50452a68325ea133e288e52dafe41))
* linting snakefmt ([eb913af](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/eb913afb7e49ff17e55631953cc490106e175b18))
* linting snakefmt ([beca1f4](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/beca1f43909b03096588724893b58f79072cbbbd))
* merge AF after normalization if 2 variants at same pos ([e328e9e](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e328e9e51b70d95bd29c3c146b3c8ed65fcda5e2))

## [0.1.2](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.1.1...v0.1.2) (2025-02-19)


### Bug Fixes

* format of the sample names ([096ddab](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/096ddab4cd0e52c65db208ffe34ccd16fa1c0aa7))
* format of the sample names ([99f16d1](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/99f16d1d3d2cc2adef14831486b485b5c264ad62))
* make multiqc output comply with stackstorm setup ([e096040](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e0960401dff1fb1367913f4966f924d7474ef684))

## 0.1.0 (2025-01-31)


### Features

* add artifacts ([41be3f8](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/41be3f8ce8046285900fd6a7476240fb7828e9f2))
* add xlsx snv, indels, known sheets ([c560975](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/c560975955ca3f3a620647bb7c6787dd61b869ab))
* merge 'develop' ([2d0e37e](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/2d0e37e5a48d5d0a7bb0a93ee89703f647788ba5))
* remove shortlist and replace with panels using bcftools include ([c6a010d](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/c6a010d4ac056ce9493667e6179ad1ef356f9e63))
* start xlsxfile ([631a733](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/631a733100e451d9b7372e4146a06c603418cb3f))


### Bug Fixes

* add missig panels[panel] ([5cdcc24](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/5cdcc247826cbaaf2e79432e134d7db13b729036))
* add missing lib, and change varible chr to chrom ([e0d4cc9](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e0d4cc926c5719eac0dd30437cab74a42be6a493))
* change from intervals to interval_list ([ee9d017](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/ee9d017d19320b885e9f46f9ea29f17a0ef1e7c5))
* correct config file  for snakemake-dry-run.yaml ([453a83c](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/453a83c01835458c351359282d0bbd8f0cfc02e9))
* lock numpy to v1.26.4 ([ef04cc8](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/ef04cc8f8e3d4da026f0d7e5b78a3143818a7973))
* make sure copy results can handle folders ([449fd98](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/449fd983087ba3da8a960aae63e82f59f3b57a20))
* update common to match poppys ([dcb6cc6](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/dcb6cc654e42c4fb6286f69b66634685ec9b0554))
* update linting test ([b8f8d55](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/b8f8d555db9f968b12f1f65fa70c25c340955760))
* update requirements.test.txt ([aac38c9](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/aac38c961dc9a3e0ae858153b8853186271534df))
* update requirements.txt ([192099a](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/192099a598e3a28b4c3b668e0204ddafdca4954e))
* update smart_open version to lower than 7.0.0 ([8ece5c9](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/8ece5c9601fa31316059cd73bf17291d5a3c26c5))
* update snakefmt.yaml ([959033b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/959033be8a356f3d8629b16d522a7c2da51dfbf5))
* update so not all col set to str when loading samples.tsv ([36207e3](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/36207e3a45ad9fe671e5c5fb683032f5d35d573a))


### Documentation

* add default addded hydra docs ([0976c04](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/0976c04ab3667d1d28beb393943d9850eb4f747e))
* schema updates ([3580b41](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/3580b419d633a8353aa23b2899c2d164f7c6b437))
