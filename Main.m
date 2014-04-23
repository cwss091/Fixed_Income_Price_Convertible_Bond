clear all;
clc;
close all;

%Parameters

DofY = 252; % trading days of one year;
Maturity = 30; % Maturity is 30 years
mu = 0.0117; % drift term of stock price process
kappa = 10.2568; % rate of reversion in Heston's mean reverting SDE
xi = 1.3283; % volatility of Heston's mean reverting process
sigma = 0.3087; % mean variance of stock price
IniS = 48.38; % today's stock price
IniV = 0.0935; % today's estimated volatility
rho = -0.5275; %correlation between stock process and volatility process

%% Heston Model of the stock price (Monte Carlo)
% Copyright @ Yuliang Li in Georgia Tech
NumberOfOuterLoop = 10;
NumberOfSamplePath = 2000;

dt = 1/DofY;
Time = dt:dt:Maturity;

VolatilityProcess = nan(length(Time), NumberOfSamplePath); % Volatility Process

temp = sqrt(dt) * (randn(length(Time), 2*NumberOfSamplePath));
dW_1 = temp(:, 1:NumberOfSamplePath);
dW_2 = temp(:, (NumberOfSamplePath+1):(2*NumberOfSamplePath));
clear temp;
dW_s = dW_1;
dW_v = rho*dW_1 + sqrt(1-rho^2)*dW_2;

clear dW_1;
clear dW_2;

VolatilityProcess  = [IniV*ones(1,NumberOfSamplePath); VolatilityProcess];

% Matrix-volatility process
for t = 2:length(Time)+1
    VolatilityProcess(t,:) = kappa*(sigma-VolatilityProcess(t-1,:))*dt + xi*sqrt(VolatilityProcess(t-1,:)).*dW_v(t-1,:);
    VolatilityProcess(t,:) = abs(VolatilityProcess(t-1,:)+VolatilityProcess(t,:));
end

VolatilityProcess = VolatilityProcess(1:end-1,:);

%clear storage

%Stock Process
dLogS = (mu-0.5*VolatilityProcess)*dt + sqrt(VolatilityProcess).*dW_s;
LogS = cumsum([log(IniS)*ones(1,NumberOfSamplePath);dLogS]);
StockProcess = exp(LogS);
%StockProcess = StockProcess(2:end,:);

%clear storage
clear dLogS;
clear LogS;
clear VolatilityProcess
clear dW_s;
clear dW_v;

%% Remarketing reset date
% Copyright @ Yuliang Li in Georgia Tech

RreDay = DofY*(5:5:25);

RemarkResetDate = zeros(5, NumberOfSamplePath);
RreThres = 100; %if average price of five days is below it, rre happen
for i=1:5
    AveFiveBeforeRreD = mean(StockProcess(RreDay(i)-5:RreDay(i)-1,:),1);
    AveFiveBeforeRreD(AveFiveBeforeRreD >= RreThres) = 0;
    AveFiveBeforeRreD(AveFiveBeforeRreD ~= 0) = 1;
    RemarkResetDate(i,:) = RreDay(i)*AveFiveBeforeRreD;
end
RemarkResetDate(RemarkResetDate == 0) = DofY*30;
% Remarketing Reser date
RemarkResetFirstDate = min(RemarkResetDate, [], 1);

clear RemarkResetDate;
clear AveFiveBeforeRreD;
clear RreDay;

    
%% Conversion Right goes into effect day
% Copyright @ Yuliang Li in Georgia Tech

ConvThres = 120; %if 20 days' prices are over this threshold, holders will have conversion right
ConvThresDay = 20;

Quarter = DofY/4:DofY/4:DofY*30;
ConvertDay = zeros(length(Quarter), NumberOfSamplePath);
for i = 1:length(Quarter)
    DaysOverThres = sum(StockProcess(Quarter(i)-30:Quarter(i)-1,:) > ConvThres);
    DaysOverThres(DaysOverThres < 20) = 0;
    DaysOverThres(DaysOverThres ~= 0) = 1;
    ConvertDay(i, :) = Quarter(i)*DaysOverThres;
end
ConvertDay(ConvertDay == 0) = DofY*30+1;
ConvertFirstDay = min(ConvertDay, [], 1);
ConvertFirstDay(ConvertFirstDay > RemarkResetFirstDate) = 0;

clear ConvertDay;
clear Quarter;

%% Remption goes into effect day
% Copyright @ Yanjie Ji in Georgia Tech 

RedeemDay = (DofY*5+2) * ones(1, NumberOfSamplePath); % After May 5, 2008. WFC can redeem its bonds
RedeemDay(RedeemDay > RemarkResetFirstDate) = 0;

%% Conversion Rate
% Conversion Rate before May 1, 2008
% Copyright @ Yanjie Ji in Georgia Tech 

ConverRate = StockProcess(1:5*DofY,:);
BaseConvPrice = 100;
ConverRate(ConverRate <= BaseConvPrice) = 10;
ConverRate(ConverRate > BaseConvPrice) = min(21.0748, 10+(ConverRate(ConverRate>BaseConvPrice)-100)./ConverRate(ConverRate>BaseConvPrice)*33.5);
Tail = repmat(ConverRate(1260-8,:),DofY*25,1);
ConvertRate = [ConverRate;Tail];

%% Interest Rate Tree
% Copyright @ Yanjie Ji in Georgia Tech 

IR_0 = 0.0113;
theta = csvread('Theta.csv');
sigma_IR = 1.4832/100;
Term_Struct_Tree = Gen_BDT_Tree(IR_0, theta, sigma_IR, Maturity, dt);

% Libor
Libor = [0.016261015	0.022066437	0.016086439	0.02005502	0.021246918	0.019194328	0.029828332	0.029876223	0.045554909	0.047633558	0.052850133	0.04334467	0.039510032	0.040086999	0.044223169	0.042316112	0.067642062	0.033759052	0.042973632	0.045498213];

%% Pricing
% Copyright @ Yanjie Ji in Georgia Tech 

tic;
Aggregate = nan(1,NumberOfSamplePath);
AAA = nan(3,NumberOfSamplePath);
for Path_index = 1: NumberOfSamplePath
    [Bond_Price, Cash_Flow] = Gen_Bond_Price_Tree(Path_index,RemarkResetFirstDate,ConvertRate,Term_Struct_Tree,Libor,dt,Maturity,1000,0.25);
    Convert_Price = Convert_Pricing(StockProcess(:,Path_index), Bond_Price, ConvertFirstDay(Path_index),ConvertRate(:,Path_index), Term_Struct_Tree,dt);
    Redeem_Price = Redeem_Pricing(Cash_Flow,Bond_Price,RedeemDay, Term_Struct_Tree,dt);
    AAA(:,Path_index)=[Bond_Price(1,1);Convert_Price;Redeem_Price];
    Aggregate(Path_index) = Bond_Price(1,1)+Convert_Price-Redeem_Price;
end
toc;
Price = mean(Aggregate)
PartPrice = mean(AAA,2)

