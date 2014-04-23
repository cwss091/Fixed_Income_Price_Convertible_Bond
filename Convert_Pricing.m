function Convert_Price = Convert_Pricing(One_Stock_Path,Bond_Price_Tree, Conv_Right,Conv_Rate, BDT_Tree,dt)
% Pricing the conversion right as a put option on bonds
% Copyright @ Yanjie Ji in Georgia Tech

if Conv_Right == 0
    Convert_Price = 0;
else
    n = Conv_Right+63+1;
    Convert_Value = nan(n,n); 
    Convert_Value(:,end) = max( One_Stock_Path(n)*Conv_Rate(n-1) - Bond_Price_Tree(1:n,n) , 0 );

    for t = (n-1):-1:1
        if t>= Conv_Right
            temp = One_Stock_Path(t)*Conv_Rate(t-1) - Bond_Price_Tree(1:t,t);
        else 
            temp = 0;
        end
    Convert_Value(1:t,t) = max (temp , ... 
        exp(-BDT_Tree(1:t,t)*dt).*( .5* Convert_Value(1:t,t+1) +  .5* Convert_Value(2:t+1,t+1) ));    
    end
Convert_Price = Convert_Value(1,1);
end
end
