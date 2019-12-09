function Preprocess2(A,b,C,d,XLB,XUB,VariableDef;TimeLim = 999999)
    t0 = time_ns();
    GUB = Inf
    nObj = size(C)[1];
    nVar = size(A)[2];

    Cuts = [];Xs = zeros(nVar); Ys = zeros(nObj)
    while true
        Lambda = abs.(randn(nObj))
        Lambda = Lambda./sum(Lambda)
        # Time commands
        # ---------------------------------------------
        Time = (time_ns() - t0)/1e9
        if TimeLim-Time < 0.001;
            return Xs, Ys, GUB, Cuts
        end
        # ---------------------------------------------
        Y,X, Model = WSM(A,b,C,d,XLB,XUB,VariableDef,Lambda,TimeLim = TimeLim-Time)
        if Y == []
            return Xs, Ys, GUB, Cuts
        end

        if GUB > prod(Y)
            GUB = prod(Y); Xs = copy(X); Ys = copy(Y);
        end
        push!(Cuts,CutType(Lambda,sum(Lambda.*Y)))
        # Time commands
        # ---------------------------------------------
        Time = (time_ns() - t0)/1e9
        if TimeLim-Time < 0.001;
            return Xs, Ys, GUB, Cuts
        end
        # ---------------------------------------------

    end
end
