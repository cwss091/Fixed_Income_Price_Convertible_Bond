function [Bond_Price_Tree,Cash_Flow] = Gen_Bond_Price_Tree(Path_Number,RemarkResetFirstDate,Conv_R,BDT_Tree,LIBOR,dt,T,F,q)
% Generate the bond price binomial tree
% Copyright @ Yanjie Ji in Georgia Tech
Payment_Factor = 0.44*ones(1,25/q+1);
for i = 2:25/q+1
    if mod(i,3) == 1
        Payment_Factor(i) = Payment_Factor(i-1)+0.02;
    else 
        Payment_Factor(i) = Payment_Factor(i-1);
    end
end

if RemarkResetFirstDate(Path_Number) == T/dt
    T_Maturity = T;
else
    T_Maturity = RemarkResetFirstDate(Path_Number)*dt;
end
Interval_Number = T_Maturity/dt;
Bond_Price_Tree = nan(Interval_Number+1, Interval_Number+1);
Cash_Flow_Tree = zeros(Interval_Number+1, Interval_Number+1);
Cash_Flow_Tree(:,end) = ones(Interval_Number+1,1)*F;
for i= T_Maturity:-q:q
    temp = i/dt+1;
    if i >=20*q
        Cash_Flow_Tree(1:temp,temp) = ones(temp,1)*Conv_R(end,Path_Number)*Payment_Factor(i/q-19);
    else 
        Cash_Flow_Tree(1:temp,temp) = ones(temp,1)*F*(LIBOR(i/q)-0.0025)/4;
    end
end
Cash_Flow_Tree(:,end) = Cash_Flow_Tree(:,end)+ones(Interval_Number+1,1)*F;
Cash_Flow = Cash_Flow_Tree(1,:);
Bond_Price_Tree(:,end) = Cash_Flow_Tree(:,end);

for t = Interval_Number:-1:1   
    Bond_Price_Tree(1:t,t) = exp(-BDT_Tree(1:t,t)*dt).*...
        ( .5* Bond_Price_Tree(1:t,t+1) +  .5* Bond_Price_Tree(2:t+1,t+1) )+...
        Cash_Flow_Tree(1:t,t);
end
end
