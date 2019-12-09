function SelectVariable(VariableDef,X_Tilde, delta1 ,delta2, MU ;
                  n1 = [], n2 = [], TAV = 0,ZZ = 0,xlb = [],xub = [],Type=1)

  if Type == 1 || isempty(Type) # Randomly Chosen
    NonIntegers = find(.!isinteger.(X_Tilde[VariableDef[1]+1:end]))
    return VariableDef[1]+NonIntegers[rand(1:length(NonIntegers))]
  elseif Type == 2 # Randomly among those which have the most value violation.
    NonIntegers = find(.!isinteger.(X_Tilde[VariableDef[1]+1:end]))
    MaxInd = findmax( min.(X_Tilde[VariableDef[1]+NonIntegers]-floor.(X_Tilde[VariableDef[1]+NonIntegers]),
                      ceil.(X_Tilde[VariableDef[1]+NonIntegers])-X_Tilde[VariableDef[1]+NonIntegers]) )[2];
    return VariableDef[1] + MaxInd ;
  elseif Type == 3
    NonIntegers = find(.!isinteger.(X_Tilde[VariableDef[1]+1:end]))
    S = MU*min.(delta1,delta2) + (1-MU)*max.(delta1,delta2)


    ind = findmax(S[NonIntegers])[2]
    # println(ind)
    return NonIntegers[ind]+VariableDef[1]

  elseif Type == 4
    NonIntegers = find(.!isinteger.(X_Tilde[VariableDef[1]+1:end]))
    S = zeros(length(NonIntegers))
    for i = 1:length(NonIntegers)
      NewXUB = copy(xub);
      NewXUB[VariableDef[1]+NonIntegers[i]] = floor(X_Tilde[VariableDef[1]+NonIntegers[i] ])
      Z2 = RSolver(A,b,C,d,xlb,NewXUB)[1];

      NewXLB = copy(xlb);
      NewXLB[VariableDef[1]+NonIntegers[i]] = ceil(X_Tilde[VariableDef[1]+NonIntegers[i]])
      Z1 = RSolver(A,b,C,d,NewXLB,xub)[1];

      if n1[NonIntegers[i]] < TAV && !isempty(Z1)
        Delta1 = Z1 - ZZ;
      else
        Delta1 = delta1[i]
      end
      if n2[NonIntegers[i]] < TAV && !isempty(Z2)
        Delta2 = Z2 - ZZ;
      else
        Delta2 = delta2[i]
      end
      S[i] = MU*min(Delta1, Delta2) + (1-MU)*max(Delta1,Delta2)
      ind = findmax(S)[2]

      return NonIntegers[ind]+VariableDef[1]

    end
  end

end
