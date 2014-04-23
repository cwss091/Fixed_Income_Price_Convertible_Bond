function BDT_Tree = Gen_BDT_Tree (r_0, theta, sigma, T, dt)
% Generate the daily interest rate binomial tree of BDT model
% Copyright @ Yanjie Ji in Georgia Tech
Interval_Number = T/dt;
BDT_Tree = nan(Interval_Number, Interval_Number);
BDT_Tree(1,1) = r_0;

for t = 1:Interval_Number-1
    
   BDT_Tree(1:t,t+1) = BDT_Tree(1:t,t)*exp( theta(t)*dt + sigma*sqrt(dt) ); 
   BDT_Tree(t+1,t+1) =  BDT_Tree(t,t)*exp( theta(t)*dt - sigma*sqrt(dt) );
    
end
end
