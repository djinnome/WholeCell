import os, sys

ht = '/msc/neurospora/src/WholeCell/bin/simulateHighthroughputExperiments/run_simulateHighthroughputExperiments.sh'
mcr = '/msc/neurospora/src/MATLAB/MATLAB_Compiler_Runtime/v80'
genesTusRxns = [
    'MG_006', 'TU_003', 'Tmk',
    'MG_023', 'TU_011', 'Fba',
    'MG_047', 'TU_027', 'MetK',
    'MG_111', 'TU_069', 'Pgi',
    'MG_272', 'TU_180', 'AceE',
    'MG_299', 'TU_203', 'Pta',
    'MG_330', 'TU_233', 'CmkA2',
    'MG_357', 'TU_260', 'AckA',
    'MG_407', 'TU_294', 'Eno',
    'MG_431', 'TU_307', 'TpiA',
    ]
parameterTypes = [
    'PromAffinity',
    'HalfLife',
    'RxnKcat'
    ]
parameterVals = {
    '05X': 0.5,
    '2X':  2.0
    } 
paramfiles = ['parameters_fluxRatio_HalfLife.mat', 'parameters_fluxRatio_PromAffinity.mat', 'parameters_fluxRatio_RxnKcat.mat']
inputdir = '/msc/neurospora/src/WholeCell/input'
outputdir = '/msc/neurospora/src/WholeCell/output'
logdir  = '/msc/neurospora/src/WholeCell/log'
bsub = 'bsub -o %s/sim-%%s-%%d.log -P %s -q %s' % (logdir , 'WholeCell', 'week')
for param in paramfiles:
    for i in range(8):
        run = "%s %s %s seed %d parameterValsPath %s/%s simPath output/sim-%s-%d.mat" % (bsub % (param[10:-4], i+1), ht, mcr, i+1, inputdir, param,  param[10:-4], i+1)
        print run
        os.system(run)
