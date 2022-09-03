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

If you use this PAF or any functions of the script for your project, please cite this work as:
```bibtex
@misc{Qualisys_OpenSim,
    title={OpenSim Workflow in Qualisys Track Manager Project Automation Framework},
    author={Rezaie, Mohammadreza},
    year={2022},
    url={https://github.com/mrrezaie/Qualisys_OpenSim}
}
```

The model is described by Rajagopal et al. (2016) with calibrated passive muscle force curves and improved hip abductor muscle paths: https://simtk.org/projects/fbmodpassivecal

The Workflow is inspired by Rajagopal sample simulation: https://simtk.org/projects/full_body

The PAF is inspired by Qualisys PAF example: https://github.com/qualisys/paf-opensim-example

Only Gait (Walking and Running) trials is supported.

Define subject Height and Mass accurately. They will be used in the analysis.

Personally I do not filter markers data. However, I fill the gaps and remove the spikes in QTM using 'Trajectory Editor' tool.

In QTM 'Project Options', make sure these settings are satisfied before start processing:
- Coordinate system is set to World (Lab) in 'Force Data'
- 3D data and Force data are checked in 'MATLAB file Export'
- Python path is defined correctly in 'Folder Options'

Here I used MAT format instead of C3D because no C3D parser (e.g EZC3D or BTK) takes Kistler COP correction (FPCOPPOLY) prameter into account. Fortunately, QTM does it. Why don't we exclude mediums? QTM is also able to import C3D files and change force plate parameters.

The markerset:
![sample](./Templates/Markerset.csv)


| Marker | Definition | Type | Static | Dynamic |
| --- | --- | --- | --- | --- |
| C7 | C7 spinous process | Anatomical | * | * |
| T10 | T10 spinous process | Anatomical | * | * |
| CLAV | Right sternoclavicular prominent | Anatomical | * | * |
| RACR | Right acromion process | Anatomical | * | * |
| LACR | Left acromion process | Anatomical | * | * |
| RLEL | Right lateral elbow (humeral lateral epicondyle) | Anatomical | * | * |
| LLEL | Left lateral elbow (humeral lateral epicondyle) | Anatomical | * | * |
| RFAradius | Right lateral wrist (radial styloid process) | Anatomical | * | * |
| LFAradius | Left lateral wrist (radial styloid process) | Anatomical | * | * |
| RASI | Right anterior superior iliac spine | Anatomical | * | * |
| LASI | Left anterior superior iliac spine | Anatomical | * | * |
| RPSI | Right posterior superior iliac spine | Anatomical | * | * |
| LPSI | Left posterior superior iliac spine | Anatomical | * | * |
| RHJC | Right hip joint center | Virtual | * |  |
| LHJC | Left hip joint center | Virtual | * |  |
| RTH1-RTH4 | 4 clusters on right thigh | Tracking |  | * |
| LTH1-LTH4 | 4 clusters on left thigh | Tracking |  | * |
| RLFC | Right lateral epicondyle of femur | Anatomical | * | * |
| LLFC | Left lateral epicondyle of femur | Anatomical | * | * |
| RMFC | Right medial epicondyle of femur | Anatomical | * |  |
| LMFC | Left medial epicondyle of femur | Anatomical | * |  |
| RKJC | Right knee joint center | Virtual | * |  |
| LKJC | Left knee joint center | Virtual | * |  |
| RTB1-RTB4 | 4 clusters on right tibia | Tracking |  | * |
| LTB1-LTB4 | 4 clusters on left tibia | Tracking |  | * |
| RLMAL | Right lateral malleolus | Anatomical | * |  |
| LLMAL | Left lateral malleolus | Anatomical | * |  |
| RMMAL | Right medial malleolus | Anatomical | * |  |
| LMMAL | Left medial malleolus | Anatomical | * |  |
| RAJC | Right ankle joint center | Virtual | * |  |
| LAJC | Left ankle joint center | Virtual | * |  |
| RCAL | Right calcaneus tuberosity | Anatomical | * | * |
| LCAL | Left calcaneus tuberosity | Anatomical | * | * |
| RMT1 | Right 1st metatarsal head | Anatomical | * | * |
| LMT1 | Left 1st metatarsal head | Anatomical | * | * |
| RMT5 | Right 5th metatarsal head | Anatomical | * | * |
| LMT5 | Left 5th metatarsal head | Anatomical | * | * |


