# lss-impute

Module to create catalogs containing imputations of missing galaxies given a full galaxy catalog containing RA, Dec for all galaxies and Z for observed galaxies.

`lssimpute` module contains several subfiles:
 * `cat.py` which contains several utility galaxy catalog management functions.
 * `dirs.py` which manages directories for io operations (Generally tries to use `PSCRATCH` on NERSC).
 * `impute.py` containing code for different imputations methods, most notably the `ImputeModel` class based around nearest neighbor models, descirbed below.
 * `plotting.py` containing many plots used for judging the behavior and quality of the imputation short of calculating 2-point statistics.

Current best method uses nearest neighbor (NN) galaxies to model how missing Zs should be drawn. Optionally, only galaxies associated with NN are imputed, in this mode, the random catalog should be downsampled relative to the number of passes.

Tested and ran on DESI Y1 first gen mock catalogs.

# Getting Started

The `bin/` directory contains several scripts for executing imputation on galaxy catalogs. Many of the python scripts support a number of command line arguments which are described using `--help`. Generally the order of operations to do an imputation is `run_imputation.py` which generates imputed galaxy samples, `make_plots.py` which makes some QA plots, `downsample_randoms.py` which uses the imputations, observed and complete catalogs to downsample randoms per `NTILE` based on the ratio given `NTILE` passes in imputation+observed vs complete catalogs and lastly `stitch_impute.py` which stitches together the catalog containing observed galaxies and the catalog containing imputed galaxies into a single catalog for 2-point scripts like pkrun.py or xirunpc.py found in desihub/LSS.

An example shell script running these steps is found in `bin/do_impute.sh` taking a command line argument for the mock number. (Modifications are needed to point it to the original mock catalogs)

# Requirements
 * astropy
 * scipy
 * numpy
 * matplotlib
 * LSS.tabulated_cosmo at (https://github.com/desihub/LSS/tree/main)

