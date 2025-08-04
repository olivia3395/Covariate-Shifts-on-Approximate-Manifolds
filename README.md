
# Phase Transition in Nonparametric Minimax Rates for Covariate Shifts on Approximate Manifolds


This repository contains simulation code for our study of **covariate shift under approximate manifold structure**, where we explore the **minimax phase transitions** in nonparametric regression across varying geometric and distributional configurations of source and target domains.


## Repository Contents

| File                    | Description                                                                                                                                                               |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `functions.R`           | Defines core oracle-style functions assuming known \$\beta\$ and \$d\$: local polynomial estimators, data generation, MSE computation.                                    |
| `SAME_DIM.R` | Simulates the **same-dimension baseline** where both source and target are in \$\mathbb{R}^D\$.                                                                           |
| `MANIFOLD.R` | Simulates the **manifold setting**, where the **target** lies on a low-dimensional manifold embedded in \$\mathbb{R}^D\$, while the **source** lies in the ambient space. |
| `ADAPTIVE_PARALLEL.R` | Runs **adaptive Lepskii-type procedures** to estimate smoothness \$\beta\$ and intrinsic dimension \$d\$ from target samples.                                             |
| `figures/`              | Directory for storing all generated simulation plots.                                                                                                                     |



## Experimental Settings

We consider three simulation regimes to investigate the impact of geometry and adaptivity under covariate shift:



### 1. Same-Dimension Setting

* **Distribution**: Both source and target samples are drawn independently from the same ambient space \$\mathbb{R}^D\$.
* **Structure**: No manifold or low-dimensional structure is assumed or exploited.



### 2. Manifold Setting (Oracle)

* **Distribution**: The **target distribution** lies on a smooth low-dimensional manifold \$\mathcal{M} \subset \mathbb{R}^D\$, generated via a mapping \$\phi(z)\$ where \$z \in \mathbb{R}^d\$.
* **Source samples**: Drawn from the full ambient space \$\mathbb{R}^D\$, potentially off-manifold.
* **Assumption**: The intrinsic dimension \$d\$ and the smoothness level \$\alpha\$ of the regression function are **known a priori**.




### 3. Adaptive Manifold Setting (Unknown \$d\$ and \$\beta\$)

* **Distribution**: Same as in the manifold setting, where the target lies on a smooth manifold and the source in the ambient space.
* **Challenge**: The intrinsic dimension \$d\$ and smoothness \$\beta\$ are **unknown and must be estimated** from the data.
* **Method**: We employ a **fully adaptive procedure**, using:

  * **Intrinsic dimension estimation** (e.g., via nearest neighbor distances),
  * **Lepskii's method** for automatic bandwidth and smoothness selection.



