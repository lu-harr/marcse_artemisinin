
## Prevalence of molecular markers of artemisinin and partner drug resistance in Africa

### Data usage

Data files in `data/clean` are provided only to allow replication of the
results reported in this repository. This data is **not** licensed for
general reuse. Users wishing to use this data for other purposes must
obtain it from: <https://surveyor.iddo.org> and [this
repo](https://github.com/Stephanie-van-Wyk/MARC_SEA_dashboard). The data
files supporting this manuscript will not be updated following
manuscript submission. The IDDO/WWARN Molecular Surveyor is a living
systematic review and is updated monthly on the Surveyor dashboard
platform.

### What the code does

There are ten models in here:

- aggregate *PfKelch13* mutation prevalence (aggregating WHO-designated
  *validated* and *candidate* markers for artemisinin partial
  resistance), together with five models of the prevalences of
  individual *PfKelch13* mutations: C469Y, A675V, R561H, P441L, and
  R622I
- *Pfcrt-K76T* prevalence
- *Pfmdr1-N86Y* prevalence
- *Pfmdr1-Y184F* prevalence
- *Pfmdr1-D1246Y* prevalence

All models take as input published data from the **IDDO Molecular
Surveyors** ([Kelch 13](https://surveyor.iddo.org/map/k13); [partner
drug
markers](https://www.iddo.org/wwarn/tracking-resistance/act-partner-drug-molecular-surveyor))
which collate published records of molecular surveillance in *Plasmodium
falciparum* malaria. The Kelch 13 model is supplemented with unpublished
data from the [WHO Malaria Threats
Map](https://apps.who.int/malaria/maps/threats/) and other sources,
retrieved from the [MARC SE-Africa antimalarial resistance
dashboard](https://www.marcse-africa.org/antimalarial-resistance-dashboard)
of Kelch 13 surveillance (see [this
repo](https://github.com/Stephanie-van-Wyk/MARC_SEA_dashboard) and [this
preprint](https://doi.org/10.1101/2025.01.07.25320158)).

All the modelling is done in R’s `greta` and `greta.gp`. I’ve added my
own kernels to `greta.gp` in this
[fork](https://github.com/lu-harr/greta.gp.st).

#### Model targets

The model of Kelch 13 mutation prevalence estimates prevalence of any
[validated and candidate
markers](https://www.who.int/tools/compendium-of-molecular-markers-for-antimalarial-drug-resistance)
of artemisinin partial resistance in the propeller region of the Kelch
13 gene of *P. falciparum*.

The partner drug marker models estimate the prevalence of mutant (e.g.,
*Pfcrt-76T*) genotypes; mutant and wildtype (e.g., *Pfcrt-K76*) are
inversely associated with reduced susceptibility to amodiaquine and
lumefantrine.

### What’s in here

There are a number of scripts in `code/`:

- `setup.R` sets up R packages, covariate data and code for prevalence
  data formatting
- `fit_betabinomial.R` runs betabinomial model fitting (written to be
  called from the command line or bash/slurm scripts in `slurm/`, with
  arguments for marker (e.g. “crt76”) and seed (an integer))
- `predict_slurm.R` makes predictions from using the outputs of
  inference (written to be called from the command line or bash/slurm
  scripts with arguments for marker, seed, model (e.g., “bb_gne”, a
  beta-binomial model with gneiting kernel), etc.)
- `vis/visualise*.R` generate figures

Secondary to these there are:

- `build_design_matrix.R`, `wrap_fit.R`, and `predict_to_raster.R` wrap
  up some functions for different bits of the workflow.
- `surveillance_effort_slurm.R` wraps up the calculation of kernel
  density estimates for smoothed annual/aggregated surveillance effort
  rasters for each model (intended to be called from the command line or
  bash/slurm scripts with arguments for marker, model, etc.).
- `stable_transmission_mask.R` cooks up a mask for prediction using
  MAP’s estimate of *P. falciparum* parasite rate for 2022.

All code to clean raw Surveyor/MARCSE dashboard data are in
`data/clean/`.

Slurm scripts for fitting, prediction, and validation can be found in
`slurm/`.

All manuscript figures can be found in `figures/`.

### See also

- Flegg et al., 2022
- Flegg et al., 2024
- Foo et al., 2024

### Supporting data

- **Proportion of Children 2 to 10 years of age showing detectable
  Plasmodium falciparum parasite**, Malaria Atlas Project, version
  202508
  - *This layer map provides estimates of the proportion of children 2
    to 10 years of age (PfPR2-10) showing detectable Plasmodium
    falciparum parasite. To produce these estimates, we apply
    geostatistical models to response datasets consisting of PR points
    and routine surveillance reports, and a rich set of geospatial
    covariates that characterise habitat for Anopheles mosquitos that
    spread the disease. The maps cover all malaria-endemic countries and
    are produced annually from year 2000 at a 5x5 km resolution.*
  - Annual gridded estimates, 2000–2024.
  - Accessed through the `malariaAtlas` R package. See [this
    paper](https://doi.org/10.1016/S0140-6736(25)00038-8).
- **Number of newly diagnosed Plasmodium falciparum cases**, Malaria
  Atlas Project, version 202508
  - *This layer map provides estimates of newly diagnosed Plasmodium
    falciparum cases, on a given year. To produce these estimates, we
    apply geostatistical models to response datasets consisting of PR
    points and routine surveillance reports, and a rich set of
    geospatial covariates that characterise habitat for Anopheles
    mosquitos that spread the disease. The maps cover all
    malaria-endemic countries and are produced annually from year 2000
    at a 5x5 km resolution.*
  - Annual gridded estimates, 2000–2024.
  - Accessed through the `malariaAtlas` R package. See [this
    paper](https://doi.org/10.1016/S0140-6736(25)00038-8).

### TODO

- Add proper refs to MAP pkg/data
- Clean up validation code

### R session info
