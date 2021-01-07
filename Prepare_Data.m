clear
addpath(genpath('./dependencies/'))
addpath(genpath('./helper/'))

% set filtering parameters
fnotch=1; % Notch bandwidth
flp=70; % Lowpass cut-off
fhp=0.1; % Highpass  cut-off
CutWindow=60; % length of ECG window for template generation in sec.
nPat=1;

% Load the examples
load('./example/s0273lre_signal.mat') % real signals
load('./example/exampleMeasurementRandom.mat') % totally random values
signal=signal(:,[1,2,7:12]); % remove the redundant leads

% prepare some variables...
ECGTemplates_Meta=struct;
measurementPoints=ionicConc.Measurements(:,1);
MeasurementsNumber=size(measurementPoints,1);
measurementPoints=round(ionicConc.Measurements(:,1)*fs); % from time in samples
lengthSnippet=numel(measurementPoints(1,1)-(CutWindow/2)*fs+1:measurementPoints(1,1)+(CutWindow/2)*fs);
ECG_window=nan(lengthSnippet,8,MeasurementsNumber);
patientIDs=patientID;

ECG_Dialysis=nan(size(signal));

% Apply pre-filtering from ECGdeli
for l=1:1:size(signal,2)
    ECG_Dialysis(:,l)=Notch_Filter(signal(:,l),fs,50,fnotch);
    ECG_Dialysis(:,l)=ECG_Low_Filter(ECG_Dialysis(:,l),fs,flp);
    ECG_Dialysis(:,l)=ECG_High_Filter(ECG_Dialysis(:,l),fs,fhp);
end

% Cut the ECG around the measurement points
for i=1:1:MeasurementsNumber
    if measurementPoints(i,1)>length(ECG_Dialysis(:,1))
        disp('Invalid Measurement Point')
    else
        if measurementPoints(i,1)-(CutWindow/2)*fs+1 < 1
            ECG_window(:,:,i)=ECG_Dialysis(1:(CutWindow)*fs,:);
        elseif measurementPoints(i,1)+(CutWindow/2)*fs+1 > length(ECG_Dialysis(:,1))
            ECG_window(:,:,i)=ECG_Dialysis(measurementPoints(i,1)-(CutWindow)*fs+1:end,:);
        else
            ECG_window(:,:,i)=ECG_Dialysis(measurementPoints(i,1)-(CutWindow/2)*fs+1:measurementPoints(i,1)+(CutWindow/2)*fs,:);
        end
    end
end

ECGTemplates_Meta.ionicConc=ionicConc;
ECGTemplates=cell(MeasurementsNumber,8);
clear ECG_Dialysis

% Build the templates
for t=1:1:MeasurementsNumber
    disp(['Starting with measurement ' num2str(t) ' of ' num2str(MeasurementsNumber)])
    [FPT_Multichannel]=Annotate_ECG_Multi(ECG_window(:,:,t),fs);
    %-------------------------------Per Lead-----------------------------------
    for q=1:1:8
        ECG_f=ECG_window(:,q,t);
        [FPT]=Annotate_ECG_Multi(ECG_f,fs,'all',FPT_Multichannel);
        FPT(:,13)=ones(size(FPT(:,13)));
        FPT=Artifact_Detection(ECG_f,fs,FPT);
        if sum(FPT(:,13)==20)>0.8*length(FPT(:,13))
            Template=[];
            posRpeak=[];
            ampRpeak=[];
        else
            [Template,posRpeak,ampRpeak,Template_matrix,FPTtemplate]=Create_Template_Class(ECG_f,fs,FPT,1);
        end
        ECGTemplates{t,q}=Template;
        ECGTemplates_Meta.posRpeak{t,q}=posRpeak;
        ECGTemplates_Meta.ampRpeak{t,q}=ampRpeak;
        ECGTemplates_Meta.FPTtemplate{t,q}=FPTtemplate;
        ECGTemplates_Meta.FPT{t,q}=FPT;
        ECGTemplates_Meta.SingleWaves{t,q}=single(Template_matrix);
    end
    
end

% Save the templates
save('Templates.mat','ECGTemplates_Meta','ECGTemplates', '-v7.3')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now extract the features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
individualTransform=0; % if 0: the first template is used to calculate the transform, 1: every template is transformed for itself

% Some needed variables
cellSize=cellfun(@size,ECGTemplates,'UniformOutput',0);
d=cell2mat(cellSize(1,:)');


% Bring the templates to one length (template length is variable...)
spanRPeak=max(cell2mat(ECGTemplates_Meta.posRpeak(1,:)))-min(cell2mat(ECGTemplates_Meta.posRpeak(1,:)));
spanLength=max(d(:,1),[],1)-min(d(:,1),[],1);
if spanLength>=spanRPeak
    ECGMat=zeros(max(d(:,1),[],1),size(ECGTemplates,2));
    for ld=1:1:size(ECGTemplates,2)
        offLength=size(ECGMat,1)-d(ld,1);
        offRpeak=max(cell2mat(ECGTemplates_Meta.posRpeak(1,:)))-ECGTemplates_Meta.posRpeak{1,ld};
        ECGMat(1+offRpeak:end-(offLength-offRpeak),ld)=ECGTemplates{1,ld};
    end
elseif spanLength<spanRPeak
    ECGMat=zeros(max(d(:,1),[],1)+spanRPeak,size(ECGTemplates,2));
    for ld=1:1:size(ECGTemplates,2)
        offLength=size(ECGMat,1)-d(ld,1);
        offRpeak=max(cell2mat(ECGTemplates_Meta.posRpeak(1,:)))-ECGTemplates_Meta.posRpeak{1,ld};
        ECGMat(1+offRpeak:end-(offLength-offRpeak),ld)=ECGTemplates{1,ld};
    end
else
    disp('Please check the inputs')
end


% Get the tranform
FPTall=cell2mat(ECGTemplates_Meta.FPTtemplate(1,:)');
TEst=median(FPTall(:,11));
boundTright=round(min(median(FPTall(:,11))+round(150*10^(-3)*fs),size(ECGMat,1)));
boundTleft=round(max(median(FPTall(:,11))-round(150*10^(-3)*fs),median(FPTall(:,6))+round(0.02*fs)));
boundRright=round(min(size(ECGMat,1),median(FPTall(:,6))+round(0.02*fs)));
boundRleft=round(max(1,median(FPTall(:,6))-round(0.02*fs)));
[~,transformationT]=leadTransformMaxT(ECGMat,boundTleft,boundTright);


% Some needed variables 
neededTemplates=(1:1:size(ECGTemplates_Meta.ionicConc.Measurements,1))';
legendvec=cell(length(neededTemplates),1);
featureMatrix=cell(1,length(neededTemplates));
concMatrix=cell(1,length(neededTemplates));


% Analyze each template
for t=1:1:length(neededTemplates)
    cellSize=cellfun(@size,ECGTemplates,'UniformOutput',0);
    d=cell2mat(cellSize(t,:)');
    spanRPeak=max(cell2mat(ECGTemplates_Meta.posRpeak(t,:)))-min(cell2mat(ECGTemplates_Meta.posRpeak(t,:)));
    spanLength=max(d(:,1),[],1)-min(d(:,1),[],1);
    if spanLength>=spanRPeak
        ECGMat=zeros(max(d(:,1),[],1),size(ECGTemplates,2));
        for ld=1:1:size(ECGTemplates,2)
            offLength=size(ECGMat,1)-d(ld,1);
            offRpeak=max(cell2mat(ECGTemplates_Meta.posRpeak(t,:)))-ECGTemplates_Meta.posRpeak{t,ld};
            ECGMat(1+offRpeak:end-(offLength-offRpeak),ld)=ECGTemplates{t,ld};
        end
    elseif spanLength<spanRPeak
        ECGMat=zeros(max(d(:,1),[],1)+spanRPeak,size(ECGTemplates,2));
        for ld=1:1:size(ECGTemplates,2)
            offLength=size(ECGMat,1)-d(ld,1);
            offRpeak=max(cell2mat(ECGTemplates_Meta.posRpeak(t,:)))-ECGTemplates_Meta.posRpeak{t,ld};
            ECGMat(1+offRpeak:end-(offLength-offRpeak),ld)=ECGTemplates{t,ld};
        end
    else
        disp('Check again!')
    end
    
    if individualTransform==1
        FPTall=cell2mat(ECGTemplates_Meta.FPTtemplate(t,:)');
        TEst=median(FPTall(:,11));
        boundTright=round(min(median(FPTall(:,11))+round(150*10^(-3)*fs),size(ECGMat,1)));
        boundTleft=round(max(median(FPTall(:,11))-round(150*10^(-3)*fs),median(FPTall(:,6))+round(0.02*fs)));
        boundRright=round(min(size(ECGMat,1),median(FPTall(:,6))+round(0.02*fs)));
        boundRleft=round(max(1,median(FPTall(:,6))-round(0.02*fs)));
        [~,transformationT]=leadTransformMaxT(ECGMat,boundTleft,boundTright);
    end
    
    % Transform the signal
    signals2analyseT=(ECGMat*transformationT)';
    FPTcell{1,1}=nan(1,13);
    FPTallLeads=cell2mat(ECGTemplates_Meta.FPTtemplate(t,:)');
    
    % Get the R peak T-signal
    dRsig=diff(abs(signals2analyseT))';
    dRsig=dRsig(max(1,min(FPTallLeads(:,6))-round(0.1*fs)):min(max(FPTallLeads(:,6))+round(0.05*fs),length(dRsig)),1);
    posRcand=find(dRsig(2:end)<0 & dRsig(1:end-1)>=0)+max(1,min(FPTallLeads(:,6))-round(0.1*fs));
    [~,tmp]=max(abs(signals2analyseT(1,posRcand)));
    FPTcell{1,1}(:,6)=posRcand(tmp,1);
    
    % Get T peak T-signal
    bndLeft=min(round(quantile(FPTallLeads(:,11),0.8))-round(0.025*fs),FPTcell{1,1}(:,6)+round(0.1*fs));
    dTsig=diff(abs(signals2analyseT(1,:)))';
    dTsig=dTsig(bndLeft:min(size(signals2analyseT,2)-1,max(FPTallLeads(:,11))+round(0.05*fs)),1);
    posTcand=find(dTsig(2:end)<0 & dTsig(1:end-1)>=0)+bndLeft;
    [~,tmp]=max(abs(signals2analyseT(1,posTcand)));
    FPTcell{1,1}(:,11)=posTcand(tmp,1);
    
    % Get P peak T-signal
    invP=sign(median(signals2analyseT(1,FPTallLeads(:,2))));
    dTsig=diff(invP*signals2analyseT(1,:))';
    dTsig=dTsig(max(1,min(FPTallLeads(:,2))-round(0.02*fs)):max(FPTallLeads(:,2))+round(0.02*fs),1);
    posPcand=find(dTsig(2:end)<0 & dTsig(1:end-1)>=0)+max(1,min(FPTallLeads(:,2))-round(0.02*fs));
    [~,tmp]=max(invP*signals2analyseT(1,posPcand));
    if ~isempty(posPcand(tmp,1))
        FPTcell{1,1}(:,2)=posPcand(tmp,1);
    else
        FPTcell{1,1}(:,2)=mean(FPTallLeads(:,2));
    end
    
    % Calculate features T-signal
    if signals2analyseT(1,FPTcell{1,1}(:,11))<0
        [featureMatrix1,featureNames]=calculateFeaturesOneBeat_git(-signals2analyseT(1,:)',fs,0,FPTcell);
    else
        [featureMatrix1,featureNames]=calculateFeaturesOneBeat_git(signals2analyseT(1,:)',fs,0,FPTcell);
    end
    
    % Save the features and concentrations
    featureMatrix{1,t}={featureMatrix1; FPTcell}; %add RR, Pon, Poff, QRSon, QRS off, and Toff to feature matrix
    concMatrix{1,t}=ECGTemplates_Meta.ionicConc.Measurements(t,:);
end

% find unique patients and sort them into one huge array
uniquePatientIDs=unique(string(patientID));
patientIDs=string(patientID);
patientECGnum=nan(size(uniquePatientIDs,1),3);
featureMatrixPerPatient=cell(size(uniquePatientIDs,1),6); %Rows: Patient, Columns: Session of the patient
for id=1:1:size(uniquePatientIDs,1)
    indFound=strfind(patientIDs,uniquePatientIDs(id));
    tmpFeat=featureMatrix(indFound,:);
    tmpConc=concMatrix(indFound,:);
    for sess=1:1:size(tmpFeat,1)
        SessConc=cell2mat(tmpConc(sess,:)');
        idxEmpt=cellfun(@isempty,tmpFeat(sess,:));
        Tfeat=cellfun(@(x)x(1,1),tmpFeat(sess,~idxEmpt)');
        featureMatrixPerPatient{id,sess}=[SessConc(:,3),SessConc(:,4),cell2mat(Tfeat')'];
    end
end

save('Estimation_Input.mat','featureMatrixPerPatient')



