clc
close all
clear all

% LMM model to simulate 3-month LIBOR from August 2003 to May 2008
% Copyright @ Xiaowen Xu in Georgia Tech

Settle = datenum('1-May-2003');
CurveTimes = [ 1 3 6 12]';
ZeroRates = [0.01320 0.01310 0.01290 0.01357]';
CurveDates = daysadd(Settle,30*CurveTimes,1);

irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);

LMMVolFunc = @(a,t) (a(1)*t + a(2)).*exp(-a(3)*t) + a(4);
LMMVolParams = [.3 -.02 .7 .14];

numRates = 40
VolFunc(1:numRates-1) = {@(t) LMMVolFunc(LMMVolParams,t)};

Beta = .08;
CorrFunc = @(i,j,Beta) exp(-Beta*abs(i-j));
Correlation = CorrFunc(meshgrid(1:numRates-1)',meshgrid(1:numRates-1),Beta);

LMM = LiborMarketModel(irdc,VolFunc,Correlation,'Period',1)
[ZeroRates, ForwardRates] = LMM.simTermStructs(20,'nTrials',1,'Tenor',3);



