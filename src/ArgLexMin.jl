function ArgLexMin(A,b,C,d,XLB,XUB,VariableDef,ObjOrd; TimeLimit = 1e+20)

    const EPS = 1e-10
    m,n = size(A);
    ObjNum = size(C)[1];
    MyModel = Model();
    ObjVal = zeros(ObjNum)
    @variable(MyModel,XLB[i] <= X[i=1:n] <= XUB[i],category = VariableType(VariableDef,i))
    for i = 1:m
        @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end

    for i in ObjOrd
        @objective(MyModel, Min, sum(C[i,j]*X[j] for j = 1:n )+d[i] )
        setsolver(MyModel, CplexSolver(CPXPARAM_ScreenOutput = 0,CPXPARAM_Threads = 1,CPXPARAM_TimeLimit = TimeLimit))
        status = solve(MyModel)
        if (i == ObjOrd[1]) && (status != :Optimal)
            return [],[]
        elseif (status != :Optimal)
            error("In ArgLexMin the first objective was optimal but we have some infeasible in the rest.")
        end
        ObjVal[i] = getobjectivevalue(MyModel)
        @constraint(MyModel,sum(C[i,j]*X[j] for j = 1:n )+d[i] <= ObjVal[i]+EPS )
    end
    XVal = getvalue(X)
    return ObjVal, XVal
end
