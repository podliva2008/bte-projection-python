from .projection.transform.scale_projection_transform import ScaleProjectionTransform
from .projection.transform.flip_vertical_projection_transform import FlipVerticalProjectionTransform
from .projection.dymaxion.bte_dymaxion_projection import BTEDymaxionProjection

BTE_PROJECTION = ScaleProjectionTransform(FlipVerticalProjectionTransform(BTEDymaxionProjection), 7318261.522857145, 7318261.522857145)