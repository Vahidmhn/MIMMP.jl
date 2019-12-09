function MinimalT2(S)
    Min = S[1].T;
    Ind = 1;
    for i = 2:length(S)
        if S[i].T <= Min
            Min = S[i].T;
            Ind = i;
        end
    end
    return Ind,Min;
end
