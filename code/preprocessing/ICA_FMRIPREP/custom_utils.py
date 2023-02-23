# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utility workflows."""
from packaging.version import parse as parseversion, Version
from pkg_resources import resource_filename as pkgr_fn

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl, afni

from templateflow.api import get as get_template

from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.fixes import (
    FixHeaderRegistration as Registration,
    FixHeaderApplyTransforms as ApplyTransforms,
    FixN4BiasFieldCorrection as N4BiasFieldCorrection,
)
from niworkflows.interfaces.images import ValidateImage, MatchHeader
from niworkflows.interfaces.masks import SimpleShowMaskRPT
from niworkflows.interfaces.registration import EstimateReferenceImage
from niworkflows.interfaces.utils import CopyXForm
from niworkflows.utils.connections import listify
from niworkflows.utils.misc import pass_dummy_scans as _pass_dummy_scans


DEFAULT_MEMORY_MIN_GB = 0.01




def init_skullstrip_bold_wf(name="skullstrip_bold_wf"):
    """
    Apply skull-stripping to a BOLD image.

    It is intended to be used on an image that has previously been
    bias-corrected with
    :py:func:`~niworkflows.func.util.init_enhance_and_skullstrip_bold_wf`

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from niworkflows.func.util import init_skullstrip_bold_wf
            wf = init_skullstrip_bold_wf()


    Inputs
    ------
    in_file : str
        BOLD image (single volume)

    Outputs
    -------
    skull_stripped_file : str
        the ``in_file`` after skull-stripping
    mask_file : str
        mask of the skull-stripped input file
    out_report : str
        reportlet for the skull-stripping

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["in_file"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["mask_file", "skull_stripped_file", "out_report"]
        ),
        name="outputnode",
    )
    skullstrip_first_pass = pe.Node(
        fsl.BET(frac=0.2, mask=True), name="skullstrip_first_pass"
    )
    apply_mask = pe.Node(fsl.ApplyMask(), name="apply_mask")
    mask_reportlet = pe.Node(SimpleShowMaskRPT(), name="mask_reportlet")

    workflow.connect([
        (inputnode, skullstrip_first_pass, [("in_file", "in_file")]),
        (skullstrip_first_pass, outputnode, [("mask_file", "mask_file")]),
        # Masked file
        (inputnode, apply_mask, [("in_file", "in_file")]),
        (skullstrip_first_pass, apply_mask, [("mask_file", "mask_file")]),
        (apply_mask, outputnode, [("out_file", "skull_stripped_file")]),
        # Reportlet
        (inputnode, mask_reportlet, [("in_file", "background_file")]),
        (skullstrip_first_pass, mask_reportlet, [("mask_file", "mask_file")]),
        (mask_reportlet, outputnode, [("out_report", "out_report")]),
    ])
    # fmt: on

    return workflow
