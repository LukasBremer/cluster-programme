import numpy as np
from datetime import date
import time as t
from argparse import ArgumentParser


from ESN import ESN, hyperparameter
import general_methods as gm
import hdf5_ESN as he
import prepare_trainingdata as pt
import vizualisation as viz

from time import time


def _get_arguments():
    """Return script arguments as dictionary
    Args:
        - None
    Returns:
        - dictionary containing script arguments
    """
    parser = ArgumentParser()
    parser.add_argument(
        "ProgressFeedback", help="Integer Representation of bool, deciding whether ProgressFeedback is returned"
    )
    parser.add_argument(
        "ExperimentOutputFolder", help="Folderpath to store Output"
    )
    parser.add_argument(
        "ExperimentInputFolder", help="Folderpath to get Data and Hyperparameterfile"
    )
    parser.add_argument(
        "HyperparameterFolder", help="Name of Hyperparameterfolder in 'ExperimentInputFolder/HyperparameterConfigs'"
    )
    parser.add_argument(
        "iterator", help="Iterator of Hyperparameterset"
    )
    """ Potentially more arguments """
    return parser.parse_args()

def EvaluateOneESN(i:int, 
                   TrainingData:np.ndarray, t_train:np.ndarray, EvaluationData:np.ndarray, t_eval:np.ndarray,
                   PATH_Out:str,
                   hyps:hyperparameter, EP:dict,
                   ProgressFeedback:bool=False, 
                   Testdimensions=None):
    """
    The EvaluateOneESN function returns a measure of the quality of a given ESN. 
    Parameter: 
        i:      itterator of TrainingData
    """
    ################################################# Prepare TrainingData #########################################################
    IN, OUT, t = pt.TimeseriesPrediction(TrainingData[i], t_train, hyps.TrainingFrames, hyps.TrainingDensity, hyps.TransientFrames, hyps.HistoricFrames_T, hyps.HistoricDensity_T)
    
    ################################################# Initiate ESN #########################################################
    e = ESN(hyps)
    
    try:
        e.Seed = EP['RandomSeed'][str(i)]
    except KeyError:
        pass
    Test = e.initiate_reservoir(DimData_In=IN.shape[-1])

    ################################################# Store ESN and RandomSeed #########################################################
    EP['RandomSeed'][str(i)] = e.Seed
    gm.yml.create(PATH_Out+'Parameter.yml', EP)
    if i == 0:
        print('\tSaved Seeds in\t\t', PATH_Out+'Parameter.yml')

    ################################### Terminate Evaluation (return zeros) if initialization was not successfull #############################
    PredictionTimes = np.zeros(EP['Evaluation']['NrEval'])
    if Test==-1:
        if ProgressFeedback:
            gm.progress_featback.printProgressBar(i, EP['Evaluation']['NrNetworks'], decimals=1, name="\tTrained and Evaluated ({}) Networks:\t".format(i))
        return PredictionTimes
    else:
        ################################################# Train ESN #########################################################################
        e.train('mean_square_with_reg', TrainingData_In=IN, TrainingData_Out=OUT)


        ################################################# Evaluate ESN #####################################################################
        for j in range(EP['Evaluation']['NrEval']):

            ############################################# Preparing Adjust and Eval Data ######################################################
            AdjustData,_,_=pt.TimeseriesPrediction(EvaluationData[j], t_eval, 0, hyps.TrainingDensity, hyps.TransientFrames, hyps.HistoricFrames_T, hyps.HistoricDensity_T)
            [EvalData], t = pt.helper_LengthDensityCorrection([EvaluationData[j]], hyps.EvaluationFrames, hyps.TrainingDensity, hyps.TransientFrames, hyps.HistoricFrames_T, hyps.HistoricDensity_T, t_train=t_eval)

            ############################################# Adjust inner nodes and Evaluate ######################################################
            e.adjust_inner_nodes(AdjustData, hyps.TransientFrames) 
            PredictionTimes[j] = e.evaluate_esn(EvalData, EvalData, t, hyps.TransientFrames, ErrorTolerance=EP['Evaluation']['ErrorTolerance'], Lyapunov=EP['TrainingData']['LLE'], Testdimension=Testdimensions)

        if ProgressFeedback:
            gm.progress_featback.printProgressBar(i, EP['Evaluation']['NrNetworks'], decimals=1, name="\tTrained and Evaluated ({}) Networks:\t".format(i))
        return PredictionTimes

def LoadData(PATH_Data:str):
    print('\tLoading Data from \t {}.'.format(PATH_Data))
    Data = np.load(PATH_Data)
    
    TrainingData = Data['TrainingData']
    t_train = Data['t_train']
    T_train = Data['T_train']
    EvaluationData = Data['EvaluationData']
    t_eval = Data['t_eval']
    T_eval = Data['T_eval']
    dt = Data['dt']

    MaxTrainSteps = t_train.shape[0]
    MaxEvalSteps = t_eval.shape[0]
    
    return TrainingData, t_train, MaxTrainSteps, EvaluationData, t_eval, MaxEvalSteps

def LoadParameter(PATH_Config:str):
    print('\tLoading Parameters from\t', PATH_Config)
    ExperimentParameter = gm.yml.read(PATH_Config)
    ExperimentParameter['Date'] = date.today()
    try:
        test_seed = ExperimentParameter['RandomSeed']['0']    
    except KeyError:
        ExperimentParameter['RandomSeed'] = {}
    hyps = hyperparameter(ExperimentParameter)
         
    return (ExperimentParameter, hyps)

def TestCompatabillity(EP:dict, hyps:hyperparameter, MaxTrainSteps:int, MaxEvalSteps:int, Meta_Data:dict, EvalNr:int, TrainNr:int):
    if EP['Evaluation']['NrNetworks'] > EvalNr or EP['Evaluation']['NrEval'] > TrainNr:
        raise ValueError("'TrainNr' or 'EvalNr' are too large for Data. Data contains", TrainNr, "TrainingData sets and", EvalNr, "EvaluationData sets")

    if (hyps.TrainingFrames+hyps.TransientFrames) * hyps.TrainingDensity > MaxTrainSteps or (hyps.EvaluationFrames+hyps.TransientFrames)*hyps.TrainingDensity>MaxEvalSteps:
        print('\tTrainingFrames:', hyps.TrainingFrames,'\nTrainingDensity:', hyps.TrainingDensity, '\nAvailable TrainingFrames', MaxTrainSteps)
        print('\n\tEvaluationFrames:', hyps.EvaluationFrames,'\nEvaluationDensity:', hyps.TrainingDensity, '\nAvailable EvaluationFrames', MaxEvalSteps,'\n\n')
        raise ValueError("'TrainingFrames' or 'EvalFrames' are too great for Data. Data contains", MaxTrainSteps, "TrainingData sets and", MaxEvalSteps, "EvaluationData sets")
    if Meta_Data['PATH_TRAININGFILE'] != EP['PATH_TRAININGFILE']:
        raise NameError('Training Data Path of Meta-File does not agree with Training Data Path of Parameter-File.')

def ExperimentParameterTesting(EP):
    try: 
        Testdimensions = np.zeros(np.shape(EP['Evaluation']['TestDimensions']), dtype=int)
        for id, td in enumerate(EP['Evaluation']['TestDimensions']):
            bools_td = (td==np.array(EP['Evaluation']['TestDimensions']))
            if np.sum(bools_td) !=1:
                print(np.sum(bools_td), bools_td)
                raise ValueError("Test dimension {} appears multiple times in Parameter['Evaluation']['TestDimensions']. ".format(td))
            bools = (td == np.array(EP['TrainingData']['KnownDimensions']))
            if np.sum(bools) !=1:
                raise ValueError("Test dimension {} is not contained in Parameter['TrainingData']['KnownDimensions']. ".format(td))
            else:
                Testdimensions[id]=int(np.where(bools)[0])
    except KeyError:
        EP['Evaluation']['TestDimensions'] = EP['TrainingData']['KnownDimensions']
        Testdimensions = np.arange(0,len(EP['TrainingData']['KnownDimensions']),1)
        print("Warning: Parameter['Evaluation']['TestDimensions'] is not defined. Parameter['TrainingData']['KnownDimensions'] will be used instead.")
    return Testdimensions

def SaveParameterAndResults(PATH_Out:str, PredictionTimes:np.ndarray, ExperimentParameter:dict, Time:float):
    
    np.savez_compressed(PATH_Out+'PredictionTimes', PredictionTimes=PredictionTimes, ComputationTime=Time)
    gm.yml.create(PATH_Out+'Parameter.yml', ExperimentParameter)
    print('\tSaved Hypeparameter in\t', PATH_Out+'Parameter.yml', '\n\tSaved Valid Times in\t', PATH_Out+'PredictionTimes.npz')

    print('\tValid Time:\t\t', np.mean(PredictionTimes),'+/-', np.std(PredictionTimes), 'Lyapunov Times')
    print("\tEvaluation Time: \t {} minutes (with {} different ESNs and {} evaluations each).\n".format(Time, ExperimentParameter['Evaluation']['NrNetworks'], ExperimentParameter['Evaluation']['NrEval']))
    
def main():

    start = t.time()
    ########################################################## Pre Experiment ###################################################################
    args = _get_arguments()

    #Define Paths
    PATH_META = args.ExperimentInputFolder+'HyperparameterConfigs/'+args.HyperparameterFolder+'/MetaParameter.yml'
    PATH_Config = args.ExperimentInputFolder+'HyperparameterConfigs/'+args.HyperparameterFolder+'/Parameter'+args.iterator+'.yml'
    Meta_Data = gm.yml.read(PATH_META)
    PATH_Out = args.ExperimentOutputFolder
    
    #Bool Progressfeedback
    ProgressFeedback = bool(int(args.ProgressFeedback))

    #Loading Input
    EP, hyps = LoadParameter(PATH_Config)
    Testdimensions = ExperimentParameterTesting(EP)
    PATH_Data = EP['PATH_TRAININGFILE']
    TrainingData, t_train, MaxTrainSteps, EvaluationData, t_eval, MaxEvalSteps = LoadData(PATH_Data)
    if ProgressFeedback:
        print('\tData in memory.')
    
    TestCompatabillity(EP, hyps, MaxTrainSteps, MaxEvalSteps, Meta_Data, EvalNr=EvaluationData.shape[0], TrainNr=TrainingData.shape[0])    
    
    #Restrict TrainingData to Known Dimensions
    TrainingData = TrainingData[:,:,EP['TrainingData']['KnownDimensions']]
    EvaluationData = EvaluationData[:,:,EP['TrainingData']['KnownDimensions']]

    #CreatingOutputs
    PredictionTimes = np.zeros((EP['Evaluation']['NrNetworks'], EP['Evaluation']['NrEval']))
    EP['TrainingData']['PATH_TRAININGFILES'] = PATH_Data
    
    
    ########################################################## Experiment #######################################################################
    if ProgressFeedback:
        print('\tTraining, saving and evaluating ESNs.')
    for i in range(EP['Evaluation']['NrNetworks']):
        PredictionTimes[i] = EvaluateOneESN(i, TrainingData, t_train, EvaluationData, t_eval, PATH_Out, hyps, EP, ProgressFeedback, Testdimensions)
    
    ########################################################## Post Experiment ###################################################################
    if ProgressFeedback:
        print("\tSaving Results.")
    SaveParameterAndResults(PATH_Out, PredictionTimes, EP, ((t.time()-start)/60))

if __name__ =='__main__':
    main()
