
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Description of all data included in this folder

**Fish stable isotope data from Moorea**
[sia_fish.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_fish.csv):

- `id`: unique id
- `year`: year of sampling
- `type`: broad sample category (here only “Fishes”)
- `species`: species name
- `tl`: total length (cm)
- `sl`: standard length (cm)
- `weight`: body mass (g)
- `d15N`: δ<sup>15</sup>N (‰)
- `d13C`: δ<sup>13</sup>C (‰)

**Sources stable isotope data from Moorea**
[sia_sources.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_sources.csv)
(includes both basal sources and invertebrates):

- `id`: unique id
- `type`: broad sample category
- `names`: taxonomic group (phylum or class) or POM (particulate organic
  matter)
- `d15N`: δ<sup>15</sup>N (‰)
- `d13C`: δ<sup>13</sup>C (‰)

**Fish trophic guilds**
[trophic_guilds.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/trophic_guilds.csv):
trophic guilds for all species in
[sia_fish.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_fish.csv):

- `family`: family name
- `species`: species name
- `trophic_guild`: fish trophic guild

To find out how the trophic guilds were assigned, see the appendix to
the paper.
