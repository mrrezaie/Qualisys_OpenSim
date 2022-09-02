# Qualisys_OpenSim
Qualisys Track Manager Project Automation Framework with OpenSim Workflow

See how it works: [LINK TO YOUTUBE]

----
### What does this PAF do?
- Generate OpenSim motion files:
    - Static trial:
        - Detremine lab coordinate system based on marker data
        - Rotate marker data to be consistent with OpenSim coordinate system (XYZ)
        - Calculate hip joint center based on Harrington et al. 2007
        - Calculate knee and ankle joint center
        - Generate TRC file
    - Gait trials:
        - Determine lab coordinate system based on GRF data (vertical, braking and propulsion peaks)
        - Rotate marker and force data to be consistent with OpenSim coordinate system (XYZ)
        - Detect which foot strikes which force plate
        - Calculate time of feet contact
        - Low-pass filter the force data
        - Generate TRC, MOT and OpenSim External Load files
- Run OpenSim tools:
    - Scale tool:
        - Scale MaxIsometricForce based on Handsfield et al. 2014
        - Set MaxContractionVelocity to 15
    - InverseKinematics tool
    - InverseDynamics tool
    - StaticOptimization tool
    - JointReaction tool

---
### Requirements
Install miniconda/minigorge/mambaforge with the following packages:
- Numpy, Scipy (https://anaconda.org/conda-forge/scikit-learn)
- OpenSim (https://anaconda.org/opensim_admin/opensim)

---
### IMPORTANT NOTES
USE IT AT YOUR OWN RISK.

Only Gait (Walking and Running) trials is supported.

In QTM Project Options, make sure these settings are satisfied:
- Coordinate system is set to World (Lab) in Force Data
- 3D data and Force data are checked in MATLAB file Export
- Python path is defined correctly in Folder Options

Here I used MAT format instead of C3D because no C3D parser (e.g EZC3D or BTK) takes Kistler COP correction (FPCOPPOLY) prameter into account. Fortunately, QTM does it. Why don't we exclude mediums?  

If you use this PAF or any functions of the script for your project, please cite this work as:

```bibtex
@misc{Qualisys_OpenSim,
    title={Qualisys Track Manager Project Automation Framework with OpenSim Workflow},
    author={Rezaie, Mohammadreza},
    year={2022},
    url={https://github.com/mrrezaie/Qualisys_OpenSim}
}
```

