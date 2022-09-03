# Qualisys_OpenSim
OpenSim Workflow in Qualisys Track Manager Project Automation Framework

See how it works: https://www.youtube.com/watch?v=Nz8CSuOh0xg

----
### What does this PAF do?
- Generate OpenSim motion files:
    - Static trial:
        - Detremine lab coordinate system based on marker data
        - Rotate marker data to be consistent with OpenSim coordinate system (XYZ)
        - Calculate hip joint center based on Harrington et al. (2007)
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
        - Scale MaxIsometricForce based on Handsfield et al. (2014)
        - Set MaxContractionVelocity to 15 (m/s)
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

The model is described by Rajagopal et al. (2016) with calibrated passive muscle force curves and improved hip abductor muscle paths: https://simtk.org/projects/fbmodpassivecal

The Workflow is inspired by Rajagopal sample simulation: https://simtk.org/projects/full_body

The PAF is inspired by Qualisys PAF example: https://github.com/qualisys/paf-opensim-example

The markerset:
![sample](./Template/Markerset.csv)

Only Gait (Walking and Running) trials is supported.

Define subject Height and Mass accurately. They will be used in the analysis.

Personally I do not filter markers data. However, I fill the gaps and remove the spikes in QTM using 'Trajectory Editor' tool.

In QTM 'Project Options', make sure these settings are satisfied before start processing:
- Coordinate system is set to World (Lab) in 'Force Data'
- 3D data and Force data are checked in 'MATLAB file Export'
- Python path is defined correctly in 'Folder Options'

Here I used MAT format instead of C3D because no C3D parser (e.g EZC3D or BTK) takes Kistler COP correction (FPCOPPOLY) prameter into account. Fortunately, QTM does it. Why don't we exclude mediums? QTM is also able to import C3D files and change force plate parameters.

If you use this PAF or any functions of the script for your project, please cite this work as:

```bibtex
@misc{Qualisys_OpenSim,
    title={OpenSim Workflow in Qualisys Track Manager Project Automation Framework},
    author={Rezaie, Mohammadreza},
    year={2022},
    url={https://github.com/mrrezaie/Qualisys_OpenSim}
}
```

