function WSM(A,b,C,d,XLB,XUB,VariableDef,Lambda;MyModel = [],TimeLim = 999999)
    if MyModel == []
        MyModel = WSM(A,b,C,d,XLB,XUB,VariableDef)
    end
    X = getindex(MyModel,:X);
    ObjNum = size(C)[1];
    m,n = size(A);
    @objective(MyModel, Min, sum(Lambda[i]*(sum(C[i,j]*X[j] for j = 1:n )+d[i])  for i = 1:ObjNum))

    setsolver(MyModel, CplexSolver(CPXPARAM_Threads = 1,CPXPARAM_ScreenOutput = 0,CPXPARAM_TimeLimit = TimeLim ))
    status = solve(MyModel)
    if status != :Optimal
        return [], [], [];
    end

    Var = getvalue(X)
    Y = zeros(ObjNum)
    for i =1:ObjNum
        Y[i] = sum(C[i,:].*Var)+d[i]
    end
    return Y, Var,MyModel;
end

function WSM(A,b,C,d,XLB,XUB,VariableDef)
    m,n = size(A);
    MyModel = Model();
    @variable(MyModel,XLB[i] <= X[i=1:n] <= XUB[i],category = VariableType(VariableDef,i))#
    for i = 1:m
        @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end

    return MyModel
end
