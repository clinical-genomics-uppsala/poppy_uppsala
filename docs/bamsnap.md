# Bamsnap

[Bamsnap](https://github.com/parklab/bamsnap) is a tool available as a Python package that generates screenshots in IGV for specific positions.
In poppy_uppsala, some adjustments were made to the original code and the custom bamsnap is [containerized](https://hub.docker.com/layers/hydragenetics/bamsnap/0.2.19/images/sha256-934518c699d724a3e949e0855c8ca40a48013e49c93b7d15c6fcf99cd2e41b5a) with Docker.

The settings for bamsnap are such that SNVs which are called in the panel and that have a VAF > 5%
are automatically treated by bamsnap.

Example of image:

![chr3_78614806](https://raw.githubusercontent.com/clinical-genomics-uppsala/poppy_uppsala/patch-readthedocs/images/chr3_78614806.png)