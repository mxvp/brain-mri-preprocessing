"""Rigid registration of DWI / ADC to a per-patient T2w reference.

Same-patient, same-session, same physical anatomy — rigid registration with
mutual-information is sufficient and matches PI-CAI baseline practice.
ADC and high-b DWI are usually acquired together (and often share the same
header geometry), but in practice they drift relative to T2w because they're
acquired separately and the patient moves.
"""
from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)


def register_to_t2(
    moving_path: Path,
    fixed_path: Path,
    interpolator: str = "linear",
) -> "object":  # SimpleITK.Image
    """Rigidly register `moving` to `fixed`. Returns a SimpleITK.Image."""
    import SimpleITK as sitk

    fixed = sitk.ReadImage(str(fixed_path), sitk.sitkFloat32)
    moving = sitk.ReadImage(str(moving_path), sitk.sitkFloat32)

    # Center the moving image on the fixed image's CoM as the initial transform.
    init_tx = sitk.CenteredTransformInitializer(
        fixed, moving, sitk.Euler3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY,
    )

    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.20)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsGradientDescent(
        learningRate=1.0, numberOfIterations=200,
        convergenceMinimumValue=1e-6, convergenceWindowSize=10,
    )
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    reg.SetInitialTransform(init_tx, inPlace=False)

    final_tx = reg.Execute(fixed, moving)
    log.debug(f"Registration metric value: {reg.GetMetricValue():.4f}, "
              f"iterations: {reg.GetOptimizerIteration()}")

    interp = {"linear": sitk.sitkLinear, "bspline": sitk.sitkBSpline,
              "nearest": sitk.sitkNearestNeighbor}[interpolator]
    resampled = sitk.Resample(moving, fixed, final_tx, interp, 0.0,
                              moving.GetPixelID())
    return resampled
