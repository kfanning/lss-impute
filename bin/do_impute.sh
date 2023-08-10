#!/bin/bash
python run_imputation.py -t ELG -s y1mock -v 0 -i ${1} -nb -p -f -ft dynamic -rb 1 -pb 11 -e -sc -o
python run_imputation.py -t LRG -s y1mock -v 0 -i ${1} -nb -p -f -ft dynamic -rb 1 -pb 11 -e -sc -o
python make_plots.py -t ELG -s y1mock -v 0 -i ${1} -p -e
python make_plots.py -t LRG -s y1mock -v 0 -i ${1} -p -e
python downsample_randoms.py -t ELG -s y1mock -v 0 -i ${1} -o --stitched
python downsample_randoms.py -t LRG -s y1mock -v 0 -i ${1} -o --stitched
python stitch_impute.py -t ELG -s y1mock -v 0 -i ${1} -o -c
python stitch_impute.py -t LRG -s y1mock -v 0 -i ${1} -o -c

