function WeakRemoverModel(A,b,C,d,XLB,XUB,VariableDef,Y_bar1,YR)

    const EPS = 1e-5;

    m,n = size(A);
    nObj = size(C)[1];
    MyModel = Model();
    @variable(MyModel,XLB[i] <= X[i=1:n] <= XUB[i],category = VariableType(VariableDef,i))#
    @variable(MyModel,Z)

    @objective(MyModel,Min, sum(sum(C[j,i]*X[i] for i = 1:n) + d[j] for j =1:nObj) )
    for i = 1:m
         @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end
    for j = 1:nObj
         @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j]  <= Y_bar1[j] + EPS)
    end

    for j = 1:nObj
         @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j] >= YR[j]  )
    end
    return MyModel;

end
