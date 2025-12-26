Code and Data for: Multi-timescale neural adaptation underlying long-term musculoskeletal reorganization
Authors: Roland Philipp, et al.
Journal: eLife (2025)
Repository: https://github.com/animalmodel/Philipp_eLife_2025
Data DOI: https://doi.org/10.5281/zenodo.18058525

📝 Overview
This repository contains the MATLAB analysis scripts used to generate the figures for the manuscript. The code analyzes EMG activity, muscle synergies, and kinematic behavior in non-human primates before and after tendon transfer surgery.

💻 System Requirements
This code was developed and tested on Microsoft Windows 11 Pro.
MATLAB Version
MATLAB R2024a (Version 24.1) or newer is recommended.
Required Toolboxes
To run the analysis scripts without errors, the following MathWorks toolboxes are required:
Signal Processing Toolbox
Statistics and Machine Learning Toolbox
Curve Fitting Toolbox
Parallel Computing Toolbox (optional)
External Dependencies
FastICA (for Independent Component Analysis)
Tensor Toolbox (Sandia National Labs)

📂 Data Availability
The full dataset required to run these scripts (including raw int16 EMG data, processed double arrays, and synergy matrices) is archived on Zenodo.
Download the data here: [INSERT ZENODO LINK HERE]
Please download the data and extract it into a folder named Data in the root directory of this repository.

⚙️ Installation & Setup
To reproduce the figures, please maintain the following folder structure:
Clone or download this repository (Codes/ folder).
Download the data from Zenodo.
Place the Data/ folder in the same root directory as Codes/.
Important: In MATLAB, right-click the Codes folder and select "Add to Path -> Selected Folders and Subfolders".
Directory Structure:


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
Main Figures
Figure5.m: Behavioral and kinematic metrics (Success rate, Retrieval time).
Figure6.m: Temporal EMG profiles and Cross-Correlation (CC) analysis.
Figure7.m: Primary Synergies (A/B) - Spatial structure & activation patterns.
Figure8.m: Secondary Synergies (C/D) - Spatial structure & activation patterns.
Figure9.m: Cross-correlation analysis of primary synergy activation.
Figure10.m: Cross-correlation analysis of secondary synergy activation.
Figure11.m: Aggregated and averaged EMG (aaEMG) analysis.
Figure12.m: Kinematic analysis of joint angles (Monkey B).
Figure13.m: Kinematic analysis showing gradual refinement of movement.
Supplementary Figures
FigureS1.m: Evolution of Time Lag at Peak CC (EMG).
FigureS2.m: Comparison of EMG Patterns (Pre-surgery vs. Post-surgery).
FigureS3.m: EMG Profiles for specific landmark days.
FigureS4.m: Cross-correlation matrices of EMG profiles.
FigureS5.m: Variance Accounted For (VAF) plots vs. Number of Synergies.
FigureS6.m: Synergy weights visualization for all sessions.
FigureS7.m: Time-varying activation profiles.
FigureS8.m: Synergy activation comparison (Pre-surgery vs. Final recovery).
FigureS9.m: Evolution of Time Lag at Peak CC (Synergy).

🛡️ License
© 2025 Roland Philipp. All rights reserved.
Permission is hereby granted to use, reproduce, and modify this code and data for academic, educational, and other non-commercial research purposes only.
This code and data are provided "as is", without warranty of any kind, express or implied. The authors shall not be liable for any damages arising from the use of this software.
