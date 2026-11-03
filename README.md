# Vascular phenotype of non-small cell lung cancer from H&E whole-slide images

Code outline for the accompanying manuscript. Vessels are segmented on routine H&E slides, measured
inside the tumour region and summarised as slide-level vascular features.

The code is provided as an outline of the methods; parameter values are given in the Methods and
Supplement of the article. No patient data, clinical tables, slide images or model weights are included.

## Pipeline

| Step | Folder | Description |
|---|---|---|
| 1 | `01_preprocessing` | Tissue detection with TRIDENT and resampling of each slide to a common resolution |
| 2 | `02_vessel_segmentation` | Vessel segmentation with SegFormer (mmsegmentation); vessel polygons per slide |
| 3 | `03_vessel_features` | Per-vessel morphology and slide-level vascular features |

## Requirements

See `requirements.txt`. TRIDENT: https://github.com/mahmoodlab/TRIDENT ;
mmsegmentation: https://github.com/open-mmlab/mmsegmentation

## Citation

Reference to be added on publication.
