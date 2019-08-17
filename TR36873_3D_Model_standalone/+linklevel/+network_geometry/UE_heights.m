function heights = UE_heights(nUE, UE_is_indoor)
% heights = UE_heights(nUE, UE_is_indoor) generate the heights for nUE
% users with location (indoor/outdoor) according to UE_is_indoor
heights = zeros(1, nUE);
for uu = 1:nUE
    if UE_is_indoor(uu)
        N_fl = 3 + unidrnd(5);
        n_fl = unidrnd(N_fl);
        heights(uu) = 3*(n_fl-1) + 1.5;
    else    
        heights(uu) = 1.5;
    end
end