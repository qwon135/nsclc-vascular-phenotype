"""Tissue detection and resampling of whole-slide images.

1. Tissue regions are detected with the TRIDENT tissue segmenter (https://github.com/mahmoodlab/TRIDENT).
2. Each slide is resampled to a common resolution; all later steps work in this pixel space,
   and the tissue contours are rescaled to it.
"""
from trident import load_wsi
from trident.segmentation_models.load import segmentation_model_factory

TARGET_MPP = ...   # common resolution (um/px), see Methods


def detect_tissue(slide_path, out_dir, device="cuda"):
    model = segmentation_model_factory(...)          # TRIDENT tissue segmenter
    wsi = load_wsi(slide_path=slide_path, lazy_init=False)
    wsi.segment_tissue(segmentation_model=model, job_dir=out_dir, device=device)
    # -> tissue contours (GeoJSON, level-0 coordinates)


def resample(slide, target_mpp=TARGET_MPP):
    scale = slide.mpp / target_mpp
    image = read_slide(slide)                        # read from the closest pyramid level
    return resize(image, scale), scale


def rescale_contours(contours, scale):
    return [[(x * scale, y * scale) for x, y in ring] for ring in contours]


def main(slide_list, out_dir):
    for slide_path in slide_list:
        contours = detect_tissue(slide_path, out_dir)
        image, scale = resample(open_slide(slide_path))
        save(image, rescale_contours(contours, scale))
