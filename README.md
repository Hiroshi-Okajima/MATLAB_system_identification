# MATLAB_system_identification

**MATLAB codes for system identification — from classical methods to multirate and periodically time-varying systems.**

This repository accompanies the blog hub article:  
**[System Identification: From Data to Dynamical Models — A Comprehensive Guide](https://blog.control-theory.com/entry/system-identification)**

Each folder contains MATLAB scripts with explanatory comments. The codes are designed to be self-contained and educational — run them directly in MATLAB to reproduce the results discussed in the corresponding blog articles and papers.

---

## Repository Structure

| Folder | Topic | Blog Article | Paper |
|--------|-------|-------------|-------|
| `01_basic_sysid/` | Introduction to system identification | [System Identification: Obtaining Dynamical Model](https://blog.control-theory.com/entry/2024/10/03/151451) | — (educational) |
| `02_subspace_n4sid/` | Subspace identification: N4SID, MOESP, CVA | [Subspace System Identification: N4SID, MOESP, and CVA](https://blog.control-theory.com/entry/subspace-identification) | — (educational) |
| `03_lptv_cyclic_sysid/` | Cyclic reformulation for LPTV systems | [Cyclic Reformulation-Based System Identification for LPTV Systems](https://blog.control-theory.com/entry/2026/03/04/232709) | [IEEE Access 2025](https://doi.org/10.1109/ACCESS.2025.3537086) |
| `04_multirate_sysid/` | Multirate system identification | [System Identification Under Multirate Sensing Environments](https://blog.control-theory.com/entry/2026/03/04/233302) | [JRM 2025 (Open Access)](https://doi.org/10.20965/jrm.2025.p1102) |
| `05_parametric_pem/` | ARX, ARMAX, OE, BJ, and the prediction error method | [Classical Parametric System Identification: ARX, ARMAX, and the PEM](https://blog.control-theory.com/entry/parametric-identification) | — (educational) |

---

## Quick Start

### Requirements

- **MATLAB** R2020b or later
- **System Identification Toolbox** (for `n4sid`, `arx`, `armax`, `pem`, `ssest`)
- **Control System Toolbox** (for `tf`, `ss`, `bode`)
- **Robust Control Toolbox** (for LMI-related examples in `03_lptv_cyclic_sysid/`)

### Running the Examples

```matlab
% Clone the repository
% >> git clone https://github.com/Hiroshi-Okajima/MATLAB_system_identification.git

% Navigate to a folder and run the main script
cd 01_basic_sysid
run('basic_sysid_demo.m')
```

Each folder contains a main script (named `*_demo.m` or `main_*.m`) that reproduces the key results.

---

## Folder Details

### 01_basic_sysid — Introduction to System Identification

A simple educational example demonstrating the basic workflow of system identification: data generation → model estimation → validation. Uses a 2nd-order discrete-time system.

**What you will learn:**
- How to create `iddata` objects in MATLAB
- Basic ARX estimation with `arx`
- Model validation with `compare` and `resid`

### 02_subspace_n4sid — Subspace Identification Methods

Demonstrates the three major subspace identification algorithms (N4SID, MOESP, CVA) on the same dataset, comparing their results. Includes model order selection via singular value analysis.

**What you will learn:**
- Using `n4sid` with different weighting options (`'auto'`, `'MOESP'`, `'CVA'`)
- Model order selection from singular value plots
- Combining `n4sid` with `pem` refinement via `ssest`
- Comparison table of identified models

### 03_lptv_cyclic_sysid — Cyclic Reformulation for LPTV Systems

Reproduces the numerical example from the IEEE Access 2025 paper. Identifies a periodically time-varying (LPTV) system using cyclic reformulation and subspace identification.

**What you will learn:**
- Construction of cycled signals from LPTV data
- Subspace identification (N4SID) of the cycled LTI system
- State coordinate transformation to recover LPTV parameters
- Verification that the identified parameters match the true system

**Paper:** H. Okajima, Y. Fujimoto, H. Oku and H. Kondo, [Cyclic Reformulation-Based System Identification for Periodically Time-Varying Systems](https://doi.org/10.1109/ACCESS.2025.3537086), IEEE Access, 2025.

### 04_multirate_sysid — Multirate System Identification

Reproduces the numerical example from the JRM 2025 paper. Identifies a discrete-time LTI plant from multirate sensor data (different output sampling rates).

**What you will learn:**
- Multirate system formulation with periodic observation matrices
- Construction of cycled signals for multirate data
- The complete 4-step identification algorithm (Algorithm 1 in the paper)
- Transfer function recovery and comparison with the true plant

**Paper:** H. Okajima, R. Furukawa and N. Matsunaga, [System Identification Under Multirate Sensing Environments](https://doi.org/10.20965/jrm.2025.p1102), Journal of Robotics and Mechatronics, Vol. 37, No. 5, pp. 1102–1112, 2025. **(Open Access)**

**See also:** [Code Ocean capsule](https://codeocean.com/capsule/3611894/tree/v1) for the original reproducible code.

### 05_parametric_pem — Classical Parametric Methods and PEM

Demonstrates the four classical parametric model structures (ARX, ARMAX, Output-Error, Box-Jenkins) and the prediction error method (PEM) on a common example. Includes model order selection via AIC/BIC and comparison with subspace identification.

**What you will learn:**
- Estimating ARX, ARMAX, OE, and BJ models in MATLAB
- Using `pem` and `ssest` for PEM-based estimation
- Model order selection with `arxstruc`, `selstruc`, AIC, BIC
- Residual analysis with `resid`
- Comparison between parametric and subspace approaches

---

## Related Resources

### Blog Articles

- **Hub:** [System Identification: From Data to Dynamical Models](https://blog.control-theory.com/entry/system-identification)
- [System Identification: Obtaining Dynamical Model](https://blog.control-theory.com/entry/2024/10/03/151451)
- [Subspace System Identification: N4SID, MOESP, and CVA](https://blog.control-theory.com/entry/subspace-identification)
- [Classical Parametric System Identification: ARX, ARMAX, and the PEM](https://blog.control-theory.com/entry/parametric-identification)
- [Cyclic Reformulation-Based System Identification for LPTV Systems](https://blog.control-theory.com/entry/2026/03/04/232709)
- [System Identification Under Multirate Sensing Environments](https://blog.control-theory.com/entry/2026/03/04/233302)
- [State Observer and State Estimation](https://blog.control-theory.com/entry/state-observer-estimation)
- [Model Error Compensator (MEC)](https://blog.control-theory.com/entry/model-error-compensator-eng)
- [LMIs and Controller Design](https://blog.control-theory.com/entry/lmi-eng)
- [Discretization of Continuous-Time Control Systems](https://blog.control-theory.com/entry/discretization-eng)

### Research Web Pages

- [www.control-theory.com](https://www.control-theory.com/en)
- [Multi-rate System](https://www.control-theory.com/en/multi-rate-system)
- [Publications](https://www.control-theory.com/en/publications)

### Video

- [YouTube: Control Engineering Channel](https://www.youtube.com/channel/UC121T0-DD2KBuqxWx2GGRkg)
- [Video Portal](https://www.portal.control-theory.com/)

### MATLAB File Exchange / Code Ocean

- [State Estimation under Multi-Rate Sensing (File Exchange)](https://jp.mathworks.com/matlabcentral/fileexchange/182941-state-estimation-under-multi-rate-sensing-ieee-access-2023)
- [Multi-Rate System Code (Code Ocean)](https://codeocean.com/capsule/3611894/tree/v1)

### Other GitHub Repositories

- [MATLAB_state_observer](https://github.com/Hiroshi-Okajima/MATLAB_state_observer) — State observer and estimation codes (observer design uses identified models)
- [MATLAB_fandamental_control_LMI](https://github.com/Hiroshi-Okajima/MATLAB_fandamental_control_LMI) — LMI-based control design
- [Robust-control-MATLAB_MEC01](https://github.com/Hiroshi-Okajima/Robust-control-MATLAB_MEC01) — Model Error Compensator

---

## Author

**Hiroshi Okajima** — Associate Professor, Graduate School of Science and Technology, Kumamoto University, Japan. Member of SICE, ISCIE, and IEEE.

- [www.control-theory.com](https://www.control-theory.com/en)
- [Blog](https://blog.control-theory.com/)
- [YouTube](https://www.youtube.com/channel/UC121T0-DD2KBuqxWx2GGRkg)
- [GitHub](https://github.com/Hiroshi-Okajima)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
