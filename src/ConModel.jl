function ConModel(A,b,C,d,XLB,XUB,y;Cuts = [])
    m,n = size(A);
    nObj = size(C)[1];
    MyModel = Model();
    @variable(MyModel, XLB[i] <= X[i=1:n] <= XUB[i])
    @variable(MyModel, Z[i=1:nObj]>=0)
    @objective(MyModel,Min,sum(1*(sum(C[j,i]*X[i] for i = 1:n) + d[j]) for j=1:nObj))
    @constraintref Ucons[1:m];
    for i = 1:m
        Ucons[i] = @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end
    for j = 1:nObj
        @constraint(MyModel, Z[j] == sum(C[j,i]*X[i] for i = 1:n) + d[j] )
    end
    for j = 2:nObj
        @constraint(MyModel, y[j]*(sum(C[1,i]*X[i] for i = 1:n) + d[1]) == y[1]*(sum(C[j,i]*X[i] for i = 1:n) + d[j]) )
        # @constraint(MyModel, y[j]*Z[1] == y[1]*Z[j] )
    end
    for c = 1:length(Cuts)
        @constraint(MyModel,sum(Cuts[c].Coeff[i]*(sum(C[i,j]*X[j] for j = 1:n )+d[i])  for i = 1:nObj) >= Cuts[c].RHS)
    end
    return MyModel;

end
