% Interpolate the yield curve of the term structure on May 1 2003
% Copyright @ Guoqiang Jin in Georgia Tech
sigma = 1.4832/100;
Settle = datenum('01-May-2003');

CurveTimes = [1/12 3/12 6/12 1:3 5 7 10 20 30]';
%Assume 30 yr yield 0.055
ZeroRates = [.0109 .0110 .0113 .0121 .0150 .0193 .0282 .0337 0.0388 0.0480 0.055]';
CurveDates = daysadd(Settle,360*CurveTimes,1);

DeltaTime = 0.5;
NewCurveTimes = (0.5:DeltaTime:30)'; 
NewCurveDates = daysadd(Settle,360*NewCurveTimes,1);
%irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates,'InterpMethod','spline');
%RateSpec = intenvset('Rates',ZeroRates,'EndDates',CurveDates,'StartDate',Settle);
%get zero curve 
Termstructure = irdc.getZeroRates(NewCurveDates); 

BondPrice = 100*exp(-Termstructure.*NewCurveTimes); 

RateTree = zeros(length(BondPrice)+1, length(BondPrice)+1);
thetaBDT = zeros(1, length(BondPrice)+1); % theta of BDT tree
RateTree(1,1) = 0.0113;

% BDT Tree construction and compute theta
% Copyright @ Guoqiang Jin and Yuliang Li in Georgia Tech

% for i = 2:length(BondPrice)+1
%     syms tet;
%     BP = sym('BP', [i i]);
%     for j = 1:i-1
%         temp = 100*exp(-RateTree(j,i-1)*exp(tet*DeltaTime+sigma*sqrt(DeltaTime))*DeltaTime);
%         BP(j,i) = temp;
%     end
%     temp = 100*exp(-RateTree(i-1,i-1)*exp(tet*DeltaTime-sigma*sqrt(DeltaTime))*DeltaTime);
%     BP(i,i) = temp;
%     for j = i-1:-1:1
%         for k = 1:j
%             temp = 0.5*(BP(k,j+1)+BP(k+1,j+1)) * exp(-RateTree(k,j)*DeltaTime);
%             BP(k,j) = temp;
%         end
%     end
%     thetaBDT(i) = double(solve(BP(1,1)-BondPrice(i),'tet'));
%     for j = 1:i-1
%         RateTree(j,i) = RateTree(j,i-1)*exp(thetaBDT(i)*DeltaTime+sigma*sqrt(DeltaTime));
%     end
%     RateTree(i,i) = RateTree(i-1,i-1)*exp(thetaBDT(i)*DeltaTime-sigma*sqrt(DeltaTime));
% end
tic
for i = 2:length(BondPrice)+1
    syms tet;
    BP = sym('BP', [i i]);
    BP(1:i-1,i) = 100*exp(-RateTree(1:i-1,i-1)*exp(tet*DeltaTime+sigma*sqrt(DeltaTime))*DeltaTime);

    BP(i,i) = 100*exp(-RateTree(i-1,i-1)*exp(tet*DeltaTime-sigma*sqrt(DeltaTime))*DeltaTime);
    
    for j = i-1:-1:1

        BP(1:j,j) = 0.5*(BP(1:j,j+1)+BP(2:j+1,j+1)) .* exp(-RateTree(1:j,j)*DeltaTime);
    end
    thetaBDT(i) = double(solve(BP(1,1)-BondPrice(i),'tet'));
    RateTree(1:i-1,i) = RateTree(1:i-1,i-1)*exp(thetaBDT(i)*DeltaTime+sigma*sqrt(DeltaTime));
    RateTree(i,i) = RateTree(i-1,i-1)*exp(thetaBDT(i)*DeltaTime-sigma*sqrt(DeltaTime));
end
toc
%%
thetaBDT_Use = thetaBDT(2:end-1);
thetaBDT_Use = [thetaBDT_Use, thetaBDT_Use(59)];
thetaBDT_Use = repmat(thetaBDT_Use,252/2,1);
[m,n] = size(thetaBDT_Use);
thetaBDT_Use = reshape(thetaBDT_Use, 1, m*n);
thetaBDT_Use = thetaBDT_Use(1:end-1);


