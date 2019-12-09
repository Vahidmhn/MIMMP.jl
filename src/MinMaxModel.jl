function MinMaxModel(A,b,C,d,XLB,XUB,VariableDef,YR)

    m,n = size(A);
    nObj = size(C)[1];
    MyModel = Model();
    @variable(MyModel,XLB[i] <= X[i=1:n] <= XUB[i],category = VariableType(VariableDef,i))#
    @variable(MyModel,Z)

    @objective(MyModel,Min,Z)
    # @constraintref Ucons[1:m];
    # @constraintref LAMBDAcons[1:nObj];
    for i = 1:m
         @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end
    for j = 1:nObj
         @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j] - Z <= YR[j] )
    end
    for j = 1:nObj
         @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j] >= YR[j])
    end
    return MyModel;

end


function MinMaxModel(A,b,C,d,XLB,XUB,VariableDef,YR; Cuts = [])

    m,n = size(A);
    nObj = size(C)[1];
    MyModel = Model();
    @variable(MyModel,XLB[i] <= X[i=1:n] <= XUB[i],category = VariableType(VariableDef,i))#
    @variable(MyModel,Z)

    @objective(MyModel,Min,Z)
    # @constraintref Ucons[1:m];
    # @constraintref LAMBDAcons[1:nObj];
    for i = 1:m
         @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end
    for j = 1:nObj
         @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j] - Z <= YR[j] )
    end
    for j = 1:nObj
         @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j] >= YR[j])
    end
    if !isempty(Cuts)
        for i in Cuts
            ZeroVals = setdiff(VariableDef[1]+VariableDef[2]+1:n,i.OneVals)
            @constraint(MyModel, sum(1-X[Int(j+VariableDef[1]+VariableDef[2])] for j in i.OneVals)
            + sum(X[Int(j)] for j = ZeroVals) >= 1)
        end
    end
    return MyModel;

end
