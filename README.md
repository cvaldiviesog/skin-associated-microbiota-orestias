# Skin microbiota of Andean *Orestias* from extreme saline springs

This repository contains the R scripts associated with the study:

**Valdivieso, C., Chavez, F. P., & Allende, M. L.**  
*Skin microbiota of Andean Orestias from extreme saline springs differs in community organization and predicted functional potential.*

## Overview

This study characterizes the skin-associated microbiota of two endemic Andean fish species, *Orestias ascotanensis* and *Orestias gloriae*, inhabiting distinct high-altitude saline springs in the Chilean Altiplano.

Skin-associated microbial communities were characterized using V3–V4 16S rRNA gene amplicon sequencing. Amplicon reads were processed using the QIIME 2 amplicon-2024.10 platform, including DADA2 denoising and taxonomic assignment against the SILVA 138 SSU database.

The downstream analyses integrate three complementary levels of microbial organization:

- **Taxonomic composition and diversity**
- **Microbial co-occurrence network structure**
- **Predicted functional potential**

Functional profiles were inferred with PICRUSt2 and analyzed with ggpicrust2. Co-occurrence networks were inferred with SPIEC-EASI and further analyzed in Cytoscape. Finally, taxonomic, network, and functional matrices were compared using Mantel and Procrustes analyses.

## Repository contents

```text
skin-associated-microbiota-orestias/
├── README.md
├── pipeline.R
└── figures.R
