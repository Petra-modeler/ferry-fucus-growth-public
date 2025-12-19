# Growth responses of Fucus vesiculosus to ferry waves

This repository contains the analysis code and documentation for a study examining tip growth dynamics of *Fucus vesiculosus* in relation to treatment (impact/control), genotype, genotype origin (impact/control).

The project is part of an ongoing research effort and associated manuscript preparation.

---

## Project overview

The aim of this study is to investigate:
- Whether tip growth (June → July) differs between **Impact** and **Control** areas
- Whether **genotype origin** influences growth patterns
- How much of the variation in growth is explained by **initial tip length**, **site**, **plot**, and **individual genotype**

Analyses are conducted using linear mixed-effects models implemented in `glmmTMB`.

---


Raw data files are **not included** in this repository (see Data availability below).

---

## Methods (brief)

- Response variables:
  - Tip weight, length and number of apical tips in July
- Fixed effects:
  - Tip weight, length and number of apical tips in June
  - Treatment (Impact vs Control)
  - Genotype origin
  - Interaction of the previous
- Random effects:
  - Site
  - Tile
  - Genotype

Models were fitted assuming a Gaussian error distribution.

---

## Data availability

Raw data are **not publicly available at this stage** of the project, as the dataset is part of an ongoing study and will be released upon publication of the associated manuscript.

---

## Software

Analyses were conducted in **R** using, among others, the following packages:
- `glmmTMB`
- `emmeans`
- `DHARMa`
- `ggeffects`
- `tidyverse`
- `ggplot2`

---

## Status

This repository reflects an **ongoing research project**.  
Content, analyses, and interpretations may change prior to publication.

---

## Contact

For questions regarding the analyses or collaboration inquiries, please contact the repository owner.
