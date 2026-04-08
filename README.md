<div align="center">

<img src="https://img.shields.io/badge/arXiv-2507.00889-b31b1b?style=for-the-badge&logo=arxiv&logoColor=white" alt="arXiv"/>
<img src="https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R"/>
<img src="https://img.shields.io/badge/Status-Active-brightgreen?style=for-the-badge" alt="Status"/>
<img src="https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge" alt="License"/>

<br/><br/>

# Phase Transition in Nonparametric Minimax Rates for Covariate Shifts on Approximate Manifolds

<br/>

> **Wang, Deb & Mukherjee (2025)**  
> *Boston University · University of Chicago*

<br/>

[📄 Read the Paper](https://arxiv.org/abs/2507.00889) 

<br/>



</div>

## 🔍 Overview

This repository contains the simulation code accompanying our paper on **nonparametric regression under covariate shift** with structured, low-dimensional data. We study a fundamental setting where:

-  A **small labeled target dataset** lies near a low-dimensional manifold
-  A **large labeled source dataset** covers the ambient high-dimensional space
-  Standard density-ratio methods **fail** because the density ratio does not exist

We establish new **minimax rates** and identify a sharp **phase transition** governed by the distance $\rho_n$ between the target domain and the underlying manifold.

<br/>



## 💡 The Phase Transition at a Glance

<div align="center">

The minimax rate of estimation transitions at the critical threshold $\kappa_n^* = \bigl(n_P^{\frac{2\beta+d}{2\beta+D}} + n_Q\bigr)^{-\frac{1}{2\beta+d}}$:

</div>

<br/>

<div align="center">

| Regime | Condition | Minimax Rate | Benefit over target-only |
|:---:|:---:|:---:|:---:|
| **Near-manifold** | $\rho_n \ll \kappa_n^*$ | $\bigl(n_P^{\frac{2\beta+d}{2\beta+D}} + n_Q\bigr)^{-\frac{2\beta}{2\beta+d}}$ | ✅ Full transfer |
| **Far-from-manifold** | $\rho_n \gg \kappa_n^*$, $d < D$ | $\bigl(n_P + n_Q\,\rho_n^{d-D}\bigr)^{-\frac{2\beta}{2\beta+D}}$ | ✅ Partial transfer |
| **Same dimension** | $d = D$ | $(n_P + n_Q)^{-\frac{2\beta}{2\beta+D}}$ | ✅ Pooled samples |

</div>

<br/>

> 🔑 **Key insight**: Even when the target is *not exactly* on a manifold, the estimator recovers the same optimal rate as the exact manifold case — as long as $\rho_n$ is below the threshold $\kappa_n^*$.

<br/>



## 📊 Key Results 

### Experiment 1 — Same Dimension Setting ($d = D = 5$)

> Source and target both live in $\mathbb{R}^5$. The proposed estimator (blue) significantly outperforms the target-only estimator (red) at interior test points, with no degradation at exterior test points (no negative transfer).

```
figures/same_dim_mse.pdf        ← MSE curves vs. nQ
figures/same_dim_boxplot.pdf    ← Boxplots at nQ = 1000
```

<div align="center">

| Interior point $x^*_{\text{Int}}$ | Exterior point $x^*_{\text{to}}$ |
|:---:|:---:|
| ✅ Proposed < Target-only | ✅ No negative transfer |
| Gap closes as $n_Q \uparrow \infty$ | Both estimators match |

</div>

<br/>

### Experiment 2 — Manifold Setting with Oracle $(d, \beta)$ ($d=2$, $D=5$)

> Target covariates lie on a 2D manifold embedded in $\mathbb{R}^5$ via a nonlinear map $\phi(z_1, z_2)$. With known smoothness $\beta = 2.5$ and intrinsic dimension $d = 2$, the proposed estimator achieves the optimal manifold rate.

```
figures/manifold_mse.pdf        ← MSE curves vs. nQ
figures/manifold_boxplot.pdf    ← Boxplots at nQ = 1000
```

**Embedding map used:**
$$\phi(z_1, z_2) = \left(\frac{z_1+1}{2},\ \frac{z_2+1}{2},\ z_1^2,\ z_2^2,\ \frac{(z_1+1)(z_2+1)}{4}\right)$$

<br/>

### Experiment 3 — Adaptive Estimation (Unknown $d$ and $\beta$)

> Same manifold setting, but $d$ and $\beta$ are estimated from data using k-NN dimension estimation and Lepskii's method. The adaptive estimator closely tracks the oracle performance.

```
figures/adaptive_mse.pdf        ← MSE comparison: Adaptive vs. Oracle vs. Target-only
```

<div align="center">

| Estimator | $d$, $\beta$ | Uses source data |
|:---:|:---:|:---:|
| **Adaptive (Pooled LPR)** | Estimated | ✅ Yes |
| Oracle (Pooled LPR) | Known | ✅ Yes |
| Oracle (Target only) | Known | ❌ No |

</div>

<br/>



## 📂 Repository Structure

```
.
├── functions.R              # Core estimation utilities (LPR, data generation, MSE)
├── SAME_DIM.R               # Experiment 1: same-dimension baseline
├── MANIFOLD.R               # Experiment 2: manifold setting with oracle (d, β)
├── ADAPTIVE_PARALLEL.R      # Experiment 3: adaptive Lepskii procedure (parallelized)
├── Loglog.R                 # Log-log convergence rate plots
└── figures/                 # Output directory for all generated plots
```

### File Descriptions

| File | Description |
|---|---|
| `functions.R` | Defines the local polynomial regression (LPR) estimator, data-generating processes, bandwidth selection, and MSE computation. Assumes known $\beta$ and $d$. |
| `SAME_DIM.R` | Simulates the same-dimension baseline ($d = D$) with source and target both uniform on hypercubes in $\mathbb{R}^D$. Tests interior and exterior evaluation points. |
| `MANIFOLD.R` | Simulates the manifold setting ($d < D$) with oracle knowledge of $d$ and $\beta$. Target generated via a nonlinear embedding $\phi: \mathbb{R}^d \to \mathbb{R}^D$. |
| `ADAPTIVE_PARALLEL.R` | Fully adaptive procedure: estimates $d$ via $k$-NN distances and selects $\beta$ via Lepskii's method. Parallelized over replications for efficiency. |
| `Loglog.R` | Produces log-log plots of MSE vs. effective sample size, overlaying empirical and theoretical convergence slopes. |

<br/>


## 🚀 Quick Start

### Prerequisites

```r
# Install required R packages
install.packages(c(
  "parallel",    # Parallel computation (ADAPTIVE_PARALLEL.R)
  "ggplot2",     # Plotting
  "dplyr",       # Data manipulation
  "tidyr"        # Data tidying
))
```

### Running the Experiments

**Experiment 1 — Same Dimension:**
```r
source("functions.R")
source("SAME_DIM.R")
# Outputs: figures/same_dim_mse.pdf, figures/same_dim_boxplot.pdf
```

**Experiment 2 — Manifold (Oracle):**
```r
source("functions.R")
source("MANIFOLD.R")
# Outputs: figures/manifold_mse.pdf, figures/manifold_boxplot.pdf
```

**Experiment 3 — Adaptive (parallelized):**
```r
source("functions.R")
source("ADAPTIVE_PARALLEL.R")
# Outputs: figures/adaptive_mse.pdf
# Note: set num_cores in the script to match your machine
```

**Log-log convergence plots:**
```r
source("Loglog.R")
# Outputs: figures/loglog_*.pdf
```

<br/>



## 🔧 Key Parameters

| Parameter | Symbol | Default | Description |
|---|---|---|---|
| Ambient dimension | $D$ | `5` | Dimension of the full covariate space |
| Intrinsic dimension | $d$ | `2` | Dimension of the target manifold |
| Hölder smoothness | $\beta$ | `2.5` | Smoothness of the regression function $f^*$ |
| Source sample size | $n_P$ | `{100, 1000, 5000, 10000}` | Number of source observations |
| Target sample size | $n_Q$ | `{100, ..., 100000}` | Number of target observations |
| Monte Carlo reps | — | `100` | Number of simulation repetitions |

<br/>



## Methodology

### Local Polynomial Regression Estimator

The estimator solves a weighted least-squares problem at the test point $x^*$:


### Optimal Bandwidth

### Adaptive Procedure (Algorithm 1)

1. **Estimate $d$** using $k$-nearest neighbor distances:

2. **Estimate $\beta$** using Lepskii's method over a discrete grid $\mathcal{B}$:

<br/>



## 📎 Citation

If you use this code in your research, please cite:

```bibtex
@article{wang2025phase,
  title     = {Phase Transition in Nonparametric Minimax Rates for Covariate Shifts on Approximate Manifolds},
  author    = {Wang, Yuyao and Deb, Nabarun and Mukherjee, Debarghya},
  journal   = {arXiv preprint arXiv:2507.00889},
  year      = {2025}
}
```

<br/>



## 📬 Contact

| Author | Affiliation | Email |
|---|---|---|
| Yuyao Wang | Boston University | yuyaow@bu.edu |
| Nabarun Deb | University of Chicago (Booth) | nabarun.deb@chicagobooth.edu |
| Debarghya Mukherjee | Boston University | mdeb@bu.edu |

<br/>



<div align="center">



</div>
