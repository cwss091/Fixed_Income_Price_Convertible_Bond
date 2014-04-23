function Redeem_Price = Redeem_Pricing(Cash_Flow,Bond_Price_Tree,Redeem_Right, BDT_Tree,dt)
% Pricing the redemption as a call option
% Copyright @ Yanjie Ji in Georgia Tech
if Redeem_Right == 0
    Redeem_Price = 0;
else
    
n = size(Bond_Price_Tree,1)-1;
Redeem_Value = nan(n,n); 
Coupon_Flow = Cash_Flow;
Coupon_Flow(end) = Coupon_Flow(end)-1000;
Last_Coupon = Coupon_Flow(end);
Redeem_Value(:,end) = max(  Bond_Price_Tree(1:n,n)-(1000+Last_Coupon*((mod(n,63)-1)/63)) , 0 );

for t = (n-1):-1:1
    if t >= 252*5+2
        Coupon_After_t = Coupon_Flow(t:end);
        Next_Coupon =  Coupon_After_t(find(Coupon_After_t,1));
        temp = Bond_Price_Tree(1:t,t)- (1000+Next_Coupon*mod(n,63)/63);
    else
        temp = 0;
    end
    Redeem_Value(1:t,t) = max (temp , ... 
        exp(-BDT_Tree(1:t,t)*dt).*( .5* Redeem_Value(1:t,t+1) +  .5* Redeem_Value(2:t+1,t+1) ) );
end

Redeem_Price = Redeem_Value(1,1);
end
end
