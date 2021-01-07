clear
close all

% Load the model with wr=1 (change in case for another wr)
load('./globalModels/estimation_model_wr1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set inputs here
% Input data must contain a cell array featureMatrixPerPatient where lines
% are the patients and columns are the sessions of the respective patient.
% In each cell, the output of the feature extraction and the concentrations
% are stored as matrix. The lines are the measurement points (#blood tests),
% the columns are: 1: K+ Concentration, 2: Ca2+ Concentration, 3 to end are
% the outputs of the feature extraction.
% For this script, you need to give a concentration for the first
% session of a patient (calibration)!

inputData='';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
load(inputData)

% Define the index of the features to use; models were trained with them
idxFeat=[8,9,10]+2;

% remove patients with less than 3 sessions
featureMatrixPerPatient=featureMatrixPerPatient(~cellfun(@isempty,featureMatrixPerPatient(:,3)),:);

% initialize different arrays needed for the fitting
featureMatrixTraining2=reshape(featureMatrixPerPatient,numel(featureMatrixPerPatient),1);
indexFeatTrain=1:1:numel(featureMatrixPerPatient);
indexFeatTrain=indexFeatTrain(~cellfun(@isempty,featureMatrixTraining2))';
featureMatrixTraining2=featureMatrixTraining2(~cellfun(@isempty,featureMatrixTraining2),1);
[orgRow,orgCol]=ind2sub(size(featureMatrixPerPatient),indexFeatTrain);

nSamples=size(featureMatrixPerPatient,1);
nIon=1;

estimConcK=cell(nSamples,1);


for vali=1:1:nSamples
    dataValidation=cell2mat(featureMatrixTraining2(orgRow==vali & orgCol~=1));
    dataCalibration=cell2mat(featureMatrixTraining2(orgRow==vali & orgCol==1));
    % Step one: Get the patient-specific calibration factor
    calibData=dataCalibration(:,idxFeat)';
    calibData=(calibData-repmat(meanTrainData,[1,size(calibData(1,:),1)]));
    calibData=(calibData./repmat(stdTrainData,[1,size(calibData(1,:),1)]));
    
    Xcalib=[calibData'.^3,calibData'.^2,calibData',ones(size(calibData',1),1)];
    targetCalibData=dataCalibration(:,nIon)';
    yCalib=targetCalibData(1,:)';
    yestCalib=Xcalib*linearRegCoeff;
    calibFac=mean((yestCalib(1:end-1)-yCalib(1:end-1)));
    
    % Step two: Apply the patient-specific calibration factor and evaluate results
    valiData=dataValidation(:,idxFeat)';
    valiData=(valiData-repmat(meanTrainData,[1,size(valiData(1,:),1)]));
    valiData=(valiData./repmat(stdTrainData,[1,size(valiData(1,:),1)]));
    
    Xvali=[valiData.^3',valiData.^2',valiData',ones(size(valiData',1),1)];
    yestVali=Xvali*linearRegCoeff;
    yestValiCali=yestVali-calibFac';
    
    % save predictions for each patient
    estimConcK{vali,1}=yestValiCali;
    
end

