function [ecgMatrixTransformed,directVec] = leadTransformMaxT(ecg_matrix,tStart,tEnd,varargin)
%LEADTRANSFORMMAXT Summary of this function goes here
%   Detailed explanation goes here
signal=ecg_matrix(tStart:tEnd,:);
if nargin>3 && strcmp(varargin{1},'norm')
    signal=ecg_matrix(tStart:tEnd,:)./repmat(max(abs(ecg_matrix(tStart:tEnd,:)),[],1),tEnd-tStart+1,1);
end
[~,b]=max(sum(abs(signal),2)); %[~,b]=max(sum((signal).^2,2));
directVec=signal(b,:)';
directVec=directVec/norm(directVec,1); %directVec=directVec/norm(directVec,2);
ecgMatrixTransformed=signal*directVec;
end

