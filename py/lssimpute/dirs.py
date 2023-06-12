import os
from pathlib import Path

#Find CSCRATCH for basedir, else use home
if os.environ['NERSC_HOST'] == 'perlmutter':
    basedir = os.environ['PSCRATCH']
else:
    basedir = os.getenv('CSCRATCH', '~/')

def get_catdir(survey, version):
    '''
    genrerally dir is basdir/catalogs/survey/version/

    for imputation survey is impute type, version is 
    '''
    version=str(version) #often passed as an int
    path = os.path.join(basedir, 'impute', 'catalogs', survey, version)
    Path(path).mkdir(parents=True, exist_ok=True)
    return path

def get_stagedir(survey, version):
    '''
    genrerally dir is basdir/staging/survey/version/

    generally used for intermediate catalogs like obs/missed
    that aren't used for statistics but are for others 
    '''
    version=str(version) #often passed as an int
    path = os.path.join(basedir, 'impute', 'staging', survey, version)
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def pk_fn(tracer, zlo, zhi, weightt=None, bint='lin', n='1', region='', impute=''):
    if region != '':
        region += '_'
    if impute != '':
        impute += '_'
    return 'pkpoles' + '_' + tracer + '_' + region + impute + str(zlo) + '_' + str(zhi) + '_' + str(weightt) + '_' + bint + str(n) + '.txt'

def pk_npy(tracer, zlo, zhi, weightt=None, bint='lin', region='', impute=''):
    if region != '':
        region += '_'
    if impute != '':
        impute += '_'
    return 'pkpoles' + '_' + tracer + '_' +region + impute + str(zlo) + '_' + str(zhi) + '_' + str(weightt) + '_' + bint + '.npy'
