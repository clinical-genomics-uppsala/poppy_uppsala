# Changelog

## [0.4.0](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.3.1...v0.4.0) (2025-12-15)


### Features

* downsample bam file to bamsnap ([cbbd30b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/cbbd30b8c34aeb569b1466b012116163cb01d7c6))
* downsampling deduplicated bam files that are input to bamsnap ([e1c3291](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e1c3291c27cfbcab6a5d9a2238ea52ae713fdff0))
* make the pipeline executable offline on miarka (BREAKING CHANGE) ([bf2f46b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/bf2f46b78f91f1247595cbd40b020de5dd768076))
* package pipeline for usage on miarka ([36bac1d](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/36bac1db30364cdbce3462bf35747574ec5dff39))


### Bug Fixes

* --containall singularity argument with proper binding ([57ef6ed](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/57ef6edc7d18bd8a4e2b4a1aea1f2e5c2d40c91d))
* clean comments ([a676733](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/a67673300c86dfef6656e25c888e3fe235b274b3))
* cleaned unwanted rebased sections ([c888dea](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/c888deacbaea471e0608a86dfd2730d11287e2df))
* config for bamsnap downsampling ([955544b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/955544bb828f131d107b42508527d447f2fd9a85))
* container for bamsnap ([a6d8fec](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/a6d8fec3a2554ce19bfee4f2ed8b6a356938af73))
* dockerhub paths to containers in config ([6807815](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/680781541ce14f56209de1a1ca43cfabbe93602d))
* incorporate changes for bamsnap downsampling ([6bf7c18](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/6bf7c1802bae0686ca3da05db20eea2451a1c8b2))
* indentation in blankline linting ([b2ef976](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/b2ef9769af9feccfeddc092e4e8061ab499efec7))
* indentation linting ([0a03b4b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/0a03b4b9f1935cdde32c6fff82378b46c077a837))
* keep intronic variants in UTBF ([b2792af](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/b2792af91b566bbe130695f38c50a7343d0294f0))
* linting and poppy uu version ([5f8ad54](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/5f8ad545f086244816170448123f6e154d4b30b4))
* linting and poppy uu version ([2a78b87](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/2a78b877dfef3fa362187be892696eca5aed9587))
* linting smkfmt ([e71c59b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/e71c59b24dc895a6a7597871d9a5ccccadb6096f))
* linting smkfmt ([c7a7c06](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/c7a7c063978fe7175974eb503b39a8765c3328d7))
* merge conflict ([dcad331](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/dcad331817b32e46e605198be43a367ea8959ebe))
* missing white space ([a597258](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/a597258a1c5d73451b5c951da5980e1833ce6855))
* Reduce max_reads in bamsnap_downsample_bam config ([9cd632d](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/9cd632d5b7b924a8bd5565cf20f2786bdef58a59))
* Reduce max_reads in bamsnap_downsample_bam config ([446fe5e](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/446fe5e576636136cc3a9f9a77e6cc9a2506cca8))
* remove echo for poppy_uppsala version in version.smk ([8adb83b](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/8adb83be75239d059fe6670fa22fd601f2ff54bd))
* requested changes ([ca5b594](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/ca5b594af803825752e473a417cbf8183f75cfc0))
* smkfmt ([ab3fe81](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/ab3fe8169d85a01ff6dfe34f869e243e98b5941b))
* time resources for pindel2vcf with containall ([f31cb54](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/f31cb544340098ebb96595e26b538d6a7b584bb5))
* try/except path to poppy gms ([4e8f972](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/4e8f97238065493ddd008ba9915ce05e03ec1646))
* update super-linter version in workflow ([3498ba6](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/3498ba64e524808e5437a3e5750da29bff56302c))
* updated poppy GMS version to use ([abf0bcf](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/abf0bcf2a535c9875c1e32cf8a238cb39e7929da))
* use downsampled files in bamsnap ([6f8773a](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/6f8773ac363103cda2737638ca833ff511e32c6b))
* versions of hydra-modules required in config ([16ec082](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/16ec0822dda7da681bbc0f32cfe4855f069f3e57))
* versions of hydra-modules required in config ([a198803](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/a198803e8ffd55e5674e2b07cb8ad0343a9d1a77))

## [0.3.1](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.3.0...v0.3.1) (2025-11-26)


### Bug Fixes

* force all similar jobs to using the same copy of the container ([fa7b546](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/fa7b546c99b1fa8b1de0e40d38dd0979561d52a1))

## [0.3.0](https://github.com/clinical-genomics-uppsala/poppy_uppsala/compare/v0.2.5...v0.3.0) (2025-11-19)


### Features

* downsampling merged deduplicated bam files that are input to bamsnap ([24ee530](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/24ee5307ce81598f756765bfaa76a5c6b7b636b1))


### Bug Fixes

* container for bamsnap ([7245c3f](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/7245c3f98f017464091b452dbfcc3e0f41cbcc5c))
* input and output file extensions in bamsnap rule ([67410b9](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/67410b9e6a2ac31a01837ec597d268c89c0cf5e8))
* poppies versions in snakemake-dry-run.yaml ([de47d4a](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/de47d4aa8a403d6db41a9496a21012784630e00f))
* Reduce max_reads in bamsnap_downsample_bam config ([10f4f1f](https://github.com/clinical-genomics-uppsala/poppy_uppsala/commit/10f4f1f39a6b136ab2cef23249ec9af7d81e826d))
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
