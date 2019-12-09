function SelectNode(Tree;VariableDef = 0, MU = 0, RO1 = 0, RO2 = 0, Type = 0, Type2Flag = 0)

    const  EPS = 0.000001
    
    # Type = (Type_of_selection, Key_Having_Feasible_Sol=1)
    if Type == 0 || isempty(Type)
        return length(Tree)
    elseif Type == 1
        return MinimalT(Tree)[1]
    elseif Type == 2
        return SelectNode(Tree,Type = Type2Flag)
    elseif Type == 3
        ExpectedImprove = zeros(length(Tree))
        for i = 1:length(Tree)
            delta1 = (ceil.(Tree[i].X[VariableDef[1]+1:end])-Tree[i].X[VariableDef[1]+1:end]).*RO1./(EPS+n1)
            delta2 = (Tree[i].X[VariableDef[1]+1:end]-floor.(Tree[i].X[VariableDef[1]+1:end])).*RO2./(EPS+n2)
            NonIntegers = find(.!isinteger.(Tree[i].X[VariableDef[1]+1:end]))
            S = MU*min.(delta1,delta2) + (1-MU)*max.(delta1,delta2)
            ind = findmax(S[NonIntegers])[2]
            ExpectedImprove[i] = max.(delta1[ind]+Tree[i].LB, delta2[ind]+Tree[i].LB)
        end
        ind = findmax(ExpectedImprove)[2]
        return ind
    elseif Type == 4
        ExpectedImprove = zeros(length(Tree));
        for i = 1:length(Tree)
            delta1 = (ceil.(Tree[i].X[VariableDef[1]+1:end])-Tree[i].X[VariableDef[1]+1:end]).*RO1./(EPS+n1)
            delta2 = (Tree[i].X[VariableDef[1]+1:end]-floor.(Tree[i].X[VariableDef[1]+1:end])).*RO2./(EPS+n2)
            delta = min.(delta1, delta2);
            NonIntegers = find(.!isinteger.(Tree[i].X[VariableDef[1]+1:end]))
            ExpectedImprove[i] = sum(delta[NonIntegers])+Tree[i].LB
        end
        ind = findmax(ExpectedImprove)[2]
        return ind
    end
end
