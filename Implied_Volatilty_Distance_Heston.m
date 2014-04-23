function Gaps = Implied_Volatilty_Distance_Heston( parameters_to_be_estimated , parameters_plugged , Strike_Price , Time_to_Maturity, Put_Price_Observed)
% Heston model parameters calibration
% Copyright @ Yuting Wu in Georgia Tech

r = parameters_plugged(1);       % risk free return
S_0 = parameters_plugged(2);     % today's price is normalized to be 100
v_0 = parameters_plugged(3);


sigma = 0.26313245 + 0.26313245*exp(parameters_to_be_estimated(1))/(1+exp(parameters_to_be_estimated(1)));
xi = 0.885508883 + 0.885508883*exp(parameters_to_be_estimated(2))/(1+exp(parameters_to_be_estimated(2)));
kappa = (xi^2)/(2*(sigma^2)) + exp(parameters_to_be_estimated(3));
rho = -0.1 - 0.9*exp(parameters_to_be_estimated(4))/(1+exp(parameters_to_be_estimated(4))); 

T = Time_to_Maturity(end); 
Number_of_Periods_Over_UnitTime = 1000;
dt = (452/252)/Number_of_Periods_Over_UnitTime;
Time = (dt:dt:T);

Number_of_SamplePaths = 5000;
Heston_Price = nan( length(Strike_Price),length(Time_to_Maturity) );

% simulating independent Brownian
temp = sqrt(dt) * (-1+2*(rand( length(Time) , 2*Number_of_SamplePaths) > 0.5));
dW_1 = temp(:,1:Number_of_SamplePaths);
dW_2 = temp(:,(Number_of_SamplePaths+1):2*Number_of_SamplePaths);

dW_S = dW_1;
dW_V = rho*dW_1 + sqrt(1-rho^2)*dW_2;

% Volatility Process Simulation
Volatility_Process = nan( length(Time)+1 , Number_of_SamplePaths );
Volatility_Process(1,:) = v_0; % initializing the volatility

for t = 1:length(Time)

    v_t = Volatility_Process(t,:);
    dv_t = kappa*( sigma^2 - v_t)*dt + xi*sqrt(v_t).*dW_V(t,:);
    Volatility_Process(t+1,:) = v_t+dv_t;

end

V = Volatility_Process(1:end-1,:); 

dLog_S = (r - 0.5*V)*dt + sqrt(V).*dW_S;

% Converting dLogS into the Stock Price Process
LogS = cumsum( [log(S_0)*ones(1,Number_of_SamplePaths);dLog_S] );
Stock_Price_Process = exp(LogS);
Stock_Price_Process = Stock_Price_Process(2:end,:); % first row is time dt, not time zero

% Compute the Heston Price and Convert the price into Implied Volatility
Temp_Keep_Heston_Price = nan(length(Strike_Price),length(Time_to_Maturity));

for k = 1:length(Strike_Price)

    for t = 1:length(Time_to_Maturity)

        T = Time_to_Maturity(t);
        K = Strike_Price(k);
        Time_Location = round(T/dt);

        S_T = Stock_Price_Process(Time_Location,:);
        Put_Payoff = max(K - S_T,0);
        Temp_Keep_Heston_Price(k,t) = exp(-r*T)*mean(Put_Payoff,2);

    end

end

Heston_Price(:,:) = Temp_Keep_Heston_Price;

Gaps = sum(sum(abs(Heston_Price - Put_Price_Observed)));
Gaps
end

