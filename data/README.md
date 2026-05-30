
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Description of all data included in this folder

**Fish stable isotope data from Moorea**
[sia_fish.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_fish.csv):

- `sample_id`: unique id
- `year`: year of sampling
- `sample_type`: broad sample category (here only “Fishes”)
- `species`: species name
- `tl`: total length (cm)
- `sl`: standard length (cm)
- `weight`: body mass (g)
- `d15N`: δ<sup>15</sup>N (‰)
- `d13C`: δ<sup>13</sup>C (‰)
- `N_percent`: nitrogen concentration in the sample (%)
- `C_percent`: carbon concentration in the sample (%)
- `CN_ratio`: C:N ratio

**Basal sources stable isotope data from Moorea**
[sia_sources.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_sources.csv):

- `sample_id`: unique id
- `sample_type`: broad sample category
- `source`: source name (POM = particulate organic matter)
- `species`: species name
- `d13C`: δ<sup>13</sup>C (‰)
- `d15N`: δ<sup>15</sup>N (‰)
- `C_percent`: carbon concentration in the sample (%)
- `N_percent`: nitrogen concentration in the sample (%)
- `CN_ratio`: C:N ratio

**Baselines stable isotope data from Moorea**
[sia_baselines.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_baselines.csv):

- `sample_id`: unique id
- `sample_type`: broad sample category (here only “Invertebrates”)
- `baseline`: baseline name
- `species`: species name
- `d13C`: δ<sup>13</sup>C (‰)
- `d15N`: δ<sup>15</sup>N (‰)
- `C_percent`: carbon concentration in the sample (%)
- `N_percent`: nitrogen concentration in the sample (%)
- `CN_ratio`: C:N ratio

**Fish trophic guilds**
[trophic_guilds.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/trophic_guilds.csv):
trophic guilds for all species in
[sia_fish.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/sia_fish.csv):

- `family`: family name
- `species`: species name
- `trophic_guild`: fish trophic guild

To find out how the trophic guilds were assigned, see the appendix to
the paper (Appendix S5).

**Priors**
[priors.csv](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data/priors.csv):
Informative Dirichlet prior for each trophic guild, including
percentages for each basal source:

- `trophic_guild`: fish trophic guild
- `algae`: percentage of green, brown, and red algae combined
- `cyanobacteria`: percentage of cyanobacteria
- `pom`: percentage of oceanic particulate organic matter

The percentage of algae was equally divide beteen green-brown algae and
red algae, and all percentages were rescaled so that their sum was equal
to the number of sources (i.e., n = 4).
