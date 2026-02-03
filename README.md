# universal-waveguide-mass-energy

MATLAB scripts for reproducing the figures in the manuscript  
**“A universal waveguide mass–energy relation for lossy one-dimensional waves in nature”**, including **4 main-text Figures** and **25 Extended Data Figures**.

## Overview
This repository accompanies the manuscript and provides reproducible MATLAB workflows for:
- Waveguide mass–energy geometry and the Cai–Smith chart representation.
- Resonance/feedback control on the complex loop plane and its power-transfer bounds.
- Waveguide-theoretic reconstruction and visualization of optical CPA / CPA-EP behavior.
- Waveguide-invariant mapping of electrochemical polarization kinetics (Au(111), Pt(111)).

## Requirements
- MATLAB (R2020a or newer recommended).
- Toolboxes: Optimization Toolbox may be required for some fitting / reconstruction scripts.
- OS: Windows / macOS / Linux should all work.

## Quick start
1. Clone/download this repository.
2. In MATLAB, set the repo root as your working directory and add paths:
   ```matlab
   addpath(genpath(pwd));
````

3. Run figure scripts (recommended naming convention):

   * Main figures: `Fig01_*.m`, `Fig02_*.m`, `Fig03_*.m`, `Fig04_*.m`
   * Extended Data: `EDFig01_*.m` ... `EDFig25_*.m`

> If your local filenames differ, simply map your scripts to the figure list below.

## Reproducing the main-text figures (Figures 1–4)

* **Fig. 1** — Generalized waveguide mass–energy relation and Cai–Smith representation with intrinsic asymmetry.
* **Fig. 2** — Feedback control of the absorbed power partition on the (κ, Θ) plane.
* **Fig. 3** — Unified resonance geometry on the Cai–Smith disc (loop-plane representation, pole proximity, ring-down, delay proxy).
* **Fig. 4** — Waveguide-invariant mapping of electrochemical polarization kinetics on Au(111) across pH (with analogous Pt(111) mapping in Extended Data).

## Reproducing the Extended Data figures (Extended Data Figs. 1–25)

### Extended Data Figs. 1–11: waveguide dispersion, bounds, and resonance archetypes

* **ED Fig. 1** — Dispersion classification based on the dimensionless RLGC electrical length.
* **ED Fig. 2** — Asymmetry-controlled snapshots and extreme envelopes in a lossy waveguide.
* **ED Fig. 3** — Universal error scaling of the low-mass approximation.
* **ED Fig. 4** — Waveguide relativity: the mass–energy relation, CPT symmetry, and impedance.
* **ED Fig. 5** — Power absorption bound in the Cai–Smith representation.
* **ED Fig. 6** — Resonance-aware useful load power landscape and its analytic maximum under fixed attenuation.
* **ED Fig. 7** — Useful load power versus one-port absorptivity on the Cai–Smith chart.
* **ED Fig. 8** — Four fundamental laws for power absorption and emissions in finite, linear, passive waveguides.
* **ED Fig. 9** — Cai–Smith map of the ring-down time and quality factor.
* **ED Fig. 10** — Resonator archetypes and the waveguide laws in the Cai–Smith representation.
* **ED Fig. 11** — Waveguide-theoretic reconstruction of optical CPA-EPs under fixed coherent driving versus SVD eigenchannel probing.

### Extended Data Figs. 12–13: optical reconstruction details (protocol + parameter inference)

* **ED Fig. 12** — Frequency-by-frequency reconstruction diagnostics for the effective propagation factor and boundary reflections (e.g., z(ω), K_eff(ω), Γ_S,eff(ω), Γ_L,eff(ω)) and consistency checks.
* **ED Fig. 13** — Protocol-dependent CPA reconstruction near the CPA frequency: fixed-at-Ω0 coherent input versus per-frequency optimal input (singular-value ratio, optimal amplitude/phase, port-resolved outputs).

### Extended Data Figs. 14–18: electrochemical mapping on Au(111)

* **ED Fig. 14** — Polarization curves and waveguide-invariant two-mode fits on Au(111) across pH.
* **ED Fig. 15** — Modal Tafel analysis of the waveguide-invariant decomposition on Au(111).
* **ED Fig. 16** — Cai–Smith plots coloured by normalized cathodic output density on Au(111).
* **ED Fig. 17** — Feedback magnitude and retention proxy versus effective overpotential on Au(111).
* **ED Fig. 18** — Storage intensity, susceptibility, and density optimum on Au(111).

### Extended Data Figs. 19–25: electrochemical mapping on Pt(111) + regime diagram

* **ED Fig. 19** — Waveguide-invariant mapping of polarization kinetics on Pt(111) across pH.
* **ED Fig. 20** — Polarization curves and waveguide-invariant two-mode fits on Pt(111) across pH.
* **ED Fig. 21** — Modal Tafel analysis on Pt(111) across pH.
* **ED Fig. 22** — Cai–Smith plots coloured by normalized cathodic output density on Pt(111).
* **ED Fig. 23** — Feedback magnitude and retention proxy versus effective overpotential on Pt(111).
* **ED Fig. 24** — Storage intensity, susceptibility, and density optimum on Pt(111).
* **ED Fig. 25** — Dimensionless current density versus dimensionless overpotential across regimes and asymmetry conditions.

## Data availability

Electrochemical datasets used in the manuscript may originate from publisher “Source Data” files and/or prior literature; please follow the manuscript references and licensing terms to obtain the raw data. The scripts in this repo are designed to reproduce the figures once the data are placed in the expected paths.

## License (GPL-3.0)

This project is released under the **GNU General Public License v3.0**.
You are free to use, modify, and redistribute this code under GPL-3.0 terms.
If you distribute a modified version, you must also release the source code under the same license and retain the license notice.

## Citation

If you use this code in academic work, please cite the associated manuscript.
