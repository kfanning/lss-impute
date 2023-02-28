import os

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
    return os.path.join(basedir, 'impute', 'catalogs', survey, version)

def get_stagedir(survey, version):
    '''
    genrerally dir is basdir/staging/survey/version/

    generally used for intermediate catalogs like obs/missed
    that aren't used for statistics but are for others 
    '''
    version=str(version) #often passed as an int
    return os.path.join(basedir, 'impute', 'staging', survey, version)
