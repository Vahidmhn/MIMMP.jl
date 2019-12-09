function MinimalT1(S)
    Min = S[1].LB;
    Ind = 1;
    for i = 2:length(S)
        if S[i].LB <= Min
            Min = S[i].LB;
            Ind = i;
        end
    end
    return Ind,Min;
end
