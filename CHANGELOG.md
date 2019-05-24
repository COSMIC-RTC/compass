- [COMPASS Change logs](#compass-change-logs)
  - [Release v4.2.0](#release-v420)
  - [Release v4.1.0](#release-v410)
  - [Release v4.0.0](#release-v400)
  - [Release v3.4](#release-v34)
  - [Release v3.3](#release-v33)
  - [Release v3.2](#release-v32)
  - [Release v3.0](#release-v30)
  - [Release v2.0](#release-v20)
  - [Release v1.1](#release-v11)

# COMPASS Change logs

## Release v4.2.0

* add pyramid focal plane visualization
* add md report generation in check.py
* update documentation generation using doxygen
* add natural integration of FP16 and complex in python
* add KECK aperture type
* better ELT Pupil generation
* drop BRAHMA, use CACAO interface instead
* computes PSF in every iteration
* add frame calibration
* remove KALMAN support
* add pyramid modulation weight
* add ```filter_TT``` parameter in centroider
* add modal integrator control law

for Internal developments:

* add custom pyr modulation script to handle natural guide dual-stars and resolved NGS
* add hunter integration optional

## Release v4.1.0

* Add multiple input/computation/output type for RTC module
* uniformize axis in widget display
* better ELT Pupil generation
* add unit tests of rtc module
* add DM petal
* add fake camera input
* debug ```load_config_from_file```
* debug ```wfs_init```

for Internal developments:

* add compile script

## Release v4.0.0

* change centroid computation using CUB library

for Internal developments:

* rename internal variables

## Release v3.4

* Naga anf Shesha are pur-python module
* CarmaWrap and SheshaWrap are pur-wrapper module using pybind
* minor debug

for Internal developments:

* rename internal variables

## Release v3.3

minor changes

## Release v3.2

* Re-up the database feature that allows to skip initialisation phase by re-using results of previous similar simulation runs

## Release v3.0

* Binding based on pyBind11
* N-faces pyramid WFS
* Introduce masked pixels centroider
* Introduce GUARDIANS package
* Introduce a way to check installation
* Debug pupil alignment (for non circular pupil)
* Shesha supervisor module improvement
* Remove the use of bytes (string instead)

## Release v2.0

* code restructuration
* Debug SH spots and PSF moves
* Debug Widget
* Fix build failure with cython 0.28
* Other minor debug

## Release v1.1

* update parameter files
* add pyr_misalignments
* add rtc_standalone
* add dm_standalone
* add supervisor between simulation and widget
