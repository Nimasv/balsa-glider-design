# Third-Party Data and Attribution

## UIUC Experimental Airfoil Data

The experimental airfoil performance data used in this project
were produced under the UIUC Low-Speed Airfoil Test program.

The analysis uses:

- `source/aero_data/SD7037D.LFT`
- `source/aero_data/SD7037D.DRG`

Both files identify the clean SD7037 (D) specimen.
The listed specimen builders are H. Stokely / B. Williams.

Official sources:

- https://m-selig.ae.illinois.edu/uiuc_lsat.html
- https://m-selig.ae.illinois.edu/pd.html

The original supporting files are retained:

- `GPL.TXT`
- `MANIFEST.TXT`
- `readme.TXT`
- `FORMAT01.TXT`

Consult these files and the official provider's notices for
the applicable redistribution and usage conditions.
No additional restrictions are imposed by this project on
the original experimental datasets.

The bundled `readme.TXT` is labeled Volume 1. It is preserved
as received and is not treated as proof of the volume provenance
of the SD7037 (D) files.

## SD7037 Airfoil Coordinates

The geometry uses `source/sd7037.dat`, with the header
`SD7037-092-88`.

Source:
https://m-selig.ae.illinois.edu/ads/coord/sd7037.dat

The coordinate file is third-party material.
Do not assume that the experimental polar-data license also
covers the coordinate file. No ownership of the original
airfoil coordinates is claimed by this project.

## Project Calculations

Aircraft-level performance, trim, and stability results are
project calculations based on experimental section data and
documented modeling assumptions.

They are not UIUC measurements of the complete glider,
and no endorsement by UIUC is implied.