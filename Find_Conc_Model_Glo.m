%% Preparation
clear
close all

% Add cvx
addpath(genpath('../../Thirdparty/'));
cvx_setup

addpath(genpath('./helper/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set inputs here
% Input data must contain a cell array featureMatrixPerPatient where lines
% are the patients and columns are the sessions of the respective patient.
% In each cell, the output of the feature extraction and the concentrations
% are stored as matrix. The lines are the measurement points (#blood tests),
% the columns are: 1: K+ Concentration, 2: Ca2+ Concentration, 3 to end are
% the outputs of the feature extraction.

inputData='./inputData';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
load(inputData)

% Define the index of the features to use; e.g. 8: T amplitude, 10: T downslope
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

error_mmol=cell(nSamples,1);
estimConcK=cell(nSamples,1);
realConcK=cell(nSamples,1);

% set the weighting ratios
wr=[0,0.5,1];
for wrn=1:1:4
    linearRegTMP=[];
    nIon=1;
    
    chosenLambda=nan(nSamples,1);
    
    
    %% Get the weight matrix
    all_conc=cell2mat(featureMatrixTraining2);
    all_conc=all_conc(:,nIon);
    if wrn<4
        cSpan=max(all_conc(:,1))-min(all_conc(:,1));
        dConc=abs(diff(sort(all_conc(:,1))));
        cRes=0.1;
        concVal=sort(all_conc(:,1));
        
        figure(1);
        h=histfit(concVal,round(cSpan/cRes),'loglogistic');
        pd = fitdist(concVal,'loglogistic');
        y = pdf(pd,unique(concVal));
        y=y-min(y);
        w=1-y;
        w=w-wr(wrn)*min(w);
        w=w/max(w);
        figure(2);
        plot(unique(concVal),w,'LineWidth',2)
        hold on;
        W=nan(size(all_conc(:,1)));
        for ij=1:1:size(W,1)
            W(ij,1)=w(all_conc(ij,1)==unique(concVal));
        end
    elseif wrn==4
        W=ones(size(all_conc(:,1)));
    end
    %% Step one: fit model with sessions for training
    dataTrain=cell2mat(featureMatrixTraining2);
    
    %Prepare (scale, etc.) data for training (inputs for training)
    trainData=dataTrain(:,idxFeat)';
    meanTrainData=min(trainData,[],2);
    trainData=(trainData-repmat(meanTrainData,[1,size(trainData(1,:),1)]));
    stdTrainData=max(trainData,[],2)-min(trainData,[],2);
    trainData=(trainData./repmat(stdTrainData,[1,size(trainData(1,:),1)]));
    
    %Prepare (scale, etc.) target data for training (outputs for training)
    targetTrainData=dataTrain(:,nIon)';
    midTargetData=0;
    scaTargetData=1;
    
    % Build Vandermonde matrix
    X2=[trainData.^3',trainData.^2',trainData',ones(length(trainData),1)];
    y2=targetTrainData';
    
    % Calculate the L-curve
    lambdaVec=logspace(-1,0.5,20);
    solNorm=nan(size(lambdaVec));
    resNorm=nan(size(lambdaVec));
    for lbd=1:1:length(lambdaVec)
        dimN=size(X2,2);
        lambda=lambdaVec(lbd);
        cvx_begin quiet
        cvx_precision best
        
        variable x(dimN)
        minimize(norm(W.*(X2*x-y2),2) ...
            + lambda.* norm(x,2))
        cvx_end
        solNorm(lbd)=norm(x,2);
        resNorm(lbd)=norm(W.*(X2*x-y2),2);
    end
    % Evaluate the L-curve
    xlcurve=log10(resNorm');
    ylcurve=log10(solNorm');
    xlinterp=linspace(min(xlcurve),max(xlcurve),200);
    ylinterp=interp1(xlcurve,ylcurve,xlinterp,'spline');
    lambdaVecInterp=interp1(1:1:length(lambdaVec),lambdaVec,linspace(1,20,200),'spline');
    curv=calculateCurvature(xlinterp,ylinterp);
    d=quantile(curv,0.8);
    d=find(curv>=d);
    [~,d2]=max(diff(d));
    d=d(d2+1); %take always the last connected region
    lambda=lambdaVecInterp(d);
    
    % Build the model
    dimN=size(X2,2);
    cvx_begin quiet
    cvx_precision best
    variable x(dimN)
    minimize(norm(W.*(X2*x-y2),2) ...
        + lambda.* norm(x,2))
    cvx_end
    linearRegCoeff=x;
    
    
    % Write out the models
    if wrn==4
        save(['./globalModels/estimation_model_wrNone'],'linearRegCoeff','meanTrainData','stdTrainData','lambda')
    else
        save(['./globalModels/estimation_model_wr' num2str(wr(wrn))],'linearRegCoeff','meanTrainData','stdTrainData','lambda')
    end

end