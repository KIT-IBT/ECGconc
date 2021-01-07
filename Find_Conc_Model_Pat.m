%% Preparation
clear
close all

% Add cvx
addpath(genpath('../../Thirdparty/'));
cvx_setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set inputs here
% Input data must contain a cell array featureMatrixPerPatient where lines 
% are the patients and columns are the sessions of the respective patient.
% In each cell, the output of the feature extraction and the concentrations 
% are stored as matrix. The lines are the measurement points (#blood tests),
% the columns are: 1: K+ Concentration, 2: Ca2+ Concentration, 3 to end are
% the outputs of the feature extraction.

inputData=''; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
load(inputData)

% Define the index of the features to use; e.g. 8: T amplitude, 10: T downslope
idxFeat=[8,10]+2;

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

error_mmol=cell(nSamples,1);
estimConcK=cell(nSamples,1);
realConcK=cell(nSamples,1);


% Go through all data in the dataset and estimate the model
for vali=1:1:nSamples
    
    %% Step one: fit model with sessions for training
    dataValidation=cell2mat(featureMatrixTraining2(orgRow==vali & ~(orgCol==1|orgCol==2)));
    dataTrain=cell2mat(featureMatrixTraining2(orgRow==vali & (orgCol==1|orgCol==2)));
    dataTrain=dataTrain(:,:);
    dataValidation=dataValidation(:,:);
    
    %Prepare (scale, etc.) data for training (inputs for training)
    trainData=dataTrain(:,idxFeat)';% take out the valith patient
    meanTrainData=min(trainData,[],2);
    trainData=(trainData-repmat(meanTrainData,[1,size(trainData(1,:),1)]));
    stdTrainData=max(trainData,[],2)-min(trainData,[],2);
    stdTrainData(stdTrainData==0)=1;
    trainData=(trainData./repmat(stdTrainData,[1,size(trainData(1,:),1)]));
    
    targetTrainData=dataTrain(:,nIon)';
    midTargetData=0; 
    scaTargetData=1;
    
    % calculate the coefficients
    X2=[trainData',ones(length(trainData),1)];
    y2=targetTrainData';
    lambda=0.01;
    dimN=size(X2,2);
    cvx_begin
    cvx_precision best
    
    variable x(dimN)
    minimize(norm((X2*x-y2),2) ...
        + lambda.* norm(x,2))
    cvx_end

    %% Step two: Apply the model to the remaining session
    valiData=dataValidation(:,idxFeat)';
    valiData=(valiData-repmat(meanTrainData,[1,size(valiData(1,:),1)]));
    valiData=(valiData./repmat(stdTrainData,[1,size(valiData(1,:),1)]));
    
    Xvali=[valiData',ones(size(valiData',1),1)];
    yestVali=Xvali*x;
    
    targetValidData=dataValidation(:,nIon)';
    yValid=targetValidData(1,:)';

    % save errors
    error_mmol{vali,1}=(yestVali-yValid);
    estimConcK{vali,1}=yestVali;
    realConcK{vali,1}=targetValidData';
    
end

