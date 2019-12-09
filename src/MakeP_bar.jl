function MakeP_bar(A,b,C,d,XLB,XUB,y;Cuts = [])
    m,n = size(A);
    nObj = size(C)[1];
    MyModel = Model();
    @variable(MyModel, XLB[i] <= X[i=1:n] <= XUB[i])
    @variable(MyModel,Z)
    @objective(MyModel,Min,Z)
    @constraintref Ucons[1:m];
    @constraintref LAMBDAcons[1:nObj];
    for i = 1:m
        Ucons[i] = @constraint(MyModel,sum(A[i,j]*X[j] for j = 1:n) >= b[i] )
    end
    for j = 1:nObj
        LAMBDAcons[j] = @constraint(MyModel,sum(C[j,i]*X[i] for i = 1:n) + d[j] - Z <= y[j] )
    end
    for c = 1:length(Cuts)
        @constraint(MyModel,sum(Cuts[c].Coeff[i]*(sum(C[i,j]*X[j] for j = 1:n )+d[i])  for i = 1:nObj) >= Cuts[c].RHS)
    end
    return MyModel,Ucons,LAMBDAcons;

end
