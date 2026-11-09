# Licence and data notes

This directory is a reproducibility package for a manuscript. It contains code,
model identifiers and aggregate results; it does not redistribute third-party CAD
models.

- **ABC dataset**: copyright of the CAD models belongs to their creators. See the
  ABC dataset licence (Onshape Terms of Use) before using the models. The STEP files
  are fetched separately and are not included here.
- **Fusion 360 Gallery**: research dataset distributed by Autodesk; see the dataset
  licence in the Fusion 360 Gallery repository.
- **MFCAD++**: CC BY; see the Queen's University Belfast dataset page.
- **Verifier source**: `verify/ssi_verifier.py` and `../../cpp/tests/cert_verify.cpp`
  are provided for inspection and reproduction of the manuscript results. No licence
  is granted for third-party data beyond the terms of the datasets above.
- **Certificate data**: the full certificate archive is on Zenodo
  (`10.5281/zenodo.23602457`); only small example certificates are included here.

Do not upload the raw STEP files, converted case files or generated certificate
archives to public repositories unless the dataset licence explicitly permits it.
