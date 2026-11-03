"""Vessel segmentation with SegFormer (mmsegmentation).

Training: pathologist-annotated vessel boundaries in regions of interest are tiled into patches
(classes: background, vessel), and a SegFormer model is trained with the standard mmsegmentation
training loop. Training settings are given in the Methods.

Inference: overlapping patches inside the tissue mask, blended vessel probability, threshold,
and vessel polygons (outer contours) per slide.
"""
import cv2
import numpy as np
from mmengine.config import Config
from mmengine.runner import Runner
from mmseg.apis import inference_model, init_model

PATCH = ...        # patch size (px)
OVERLAP = ...
THRESHOLD = ...
MIN_AREA = ...     # smallest polygon kept (px^2)


def train(config_path, work_dir):
    cfg = Config.fromfile(config_path)               # SegFormer, 2 classes
    cfg.work_dir = work_dir
    Runner.from_cfg(cfg).train()


def predict_slide(model, image, tissue_mask):
    prob = np.zeros(image.shape[:2], np.float32)
    weight = np.zeros_like(prob)
    window = blending_window(PATCH)
    for x, y in patch_grid(tissue_mask, PATCH, OVERLAP):
        p = inference_model(model, image[y:y + PATCH, x:x + PATCH]).seg_logits.data.softmax(0)[1]
        prob[y:y + PATCH, x:x + PATCH] += p.cpu().numpy() * window
        weight[y:y + PATCH, x:x + PATCH] += window
    mask = (prob / np.maximum(weight, 1e-6) > THRESHOLD) & tissue_mask
    return mask.astype(np.uint8)


def mask_to_polygons(mask):
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    return [c.squeeze() for c in contours if cv2.contourArea(c) >= MIN_AREA]


def main(config_path, checkpoint, slides):
    model = init_model(config_path, checkpoint)
    for image, tissue_mask in slides:
        save(mask_to_polygons(predict_slide(model, image, tissue_mask)))   # GeoJSON per slide
