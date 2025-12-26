# Code and Data for: Multi-timescale neural adaptation underlying long-term musculoskeletal reorganization

**Authors:** Roland Philipp, et al.
**Journal:** eLife (2025)
**Repository:** [https://github.com/animalmodel/Philipp_eLife_2025](https://github.com/animalmodel/Philipp_eLife_2025)

## Overview
This repository contains the MATLAB analysis scripts used to generate the figures for the manuscript. The code analyzes EMG activity and muscle synergies in non-human primates before and after tendon transfer surgery.

## 📂 Data Availability
The dataset (EMG and Synergy matrices) required to run these scripts will be made available via Zenodo/GitHub shortly.

## ⚙️ Installation & Setup
To run these scripts, you must maintain the following folder structure:

1.  Download the `Codes/` folder from this repository.
2.  Place the `Data/` folder (once available) in the same root directory as `Codes/`.

**Directory Structure:**
```text
/Philipp_eLife_2025
├── Codes/
│   ├── Figure5.m
│   ├── ... (other scripts)
├── Data/
│   ├── emg/
│   ├── synergy/
│   ├── behavior/
│   ├── kinematics/
│
└── README.md


📜 List of Scripts (Current Release)
Figure5.m: Behavioral and kinematic metrics
Figure6.m: Temporal EMG profiles and CC
Figure7.m: Primary Synergies (A/B) - Spatial structure & activation.
Figure8.m: Secondary Synergies (C/D) - Spatial structure & activation.
Figure9.m: Cross-correlation analysis of primary synergy activation
Figure10.m: Cross-correlation analysis of secondary synergy activation
Figure11.m: Aggregated and averaged EMG (aaEMG)
Figure12.m: Kinematic analysis of joint angles (Monkey B)
Figure13.m: Kinematic analysis (gradual refinement)
FigureS1.m: Evolution of Time Lag at Peak CC (EMG)
FigureS2.m: EMG Patterns (Pre-surgery vs. Post-surgery).
FigureS3.m: EMG Profiles for landmark days.
FigureS4.m: Cross-correlation of EMG profiles.
FigureS5.m: VAF plots
FigureS6.m: Synergy weights for all sessions.
FigureS7.m: Time-varying activation profiles.
FigureS8.m: Synergy activation comparison (Pre vs. Final).
FigureS9.m: Evolution of Time Lag at Peak CC (Synergy)

