function WeightedExplore(A,b,C,d,XLB,XUB,VariableDef,YR,TimeLim = 1e10)
        const EPS = 1e-5;
        t0 = time_ns();
        nObj = length(YR)

        MMModel= WeightedMinMaxModel(A,b,C,d,XLB,XUB,VariableDef,YR)

        # Time commands ---------------------------------------------------
        if (TimeLim-(time_ns() - t0)/1e9) < 0.001;
            return [],[],[],false
        end
        # -----------------------------------------------------------------

        setsolver(MMModel,CplexSolver(CPXPARAM_Threads = 1,CPXPARAM_ScreenOutput = 0,CPXPARAM_TimeLimit = TimeLim-(time_ns() - t0)/1e9));
        status = solve(MMModel);

        # Time commands ---------------------------------------------------
        if (status != :Optimal) | (TimeLim-(time_ns() - t0)/1e9 < 0.001)
            return [],[],[],false
        end
        # -----------------------------------------------------------------

        Var = getvalue(getVar(MMModel,:X))
        Y_bar1 = zeros(size(YR))
        for j = 1:nObj
            Y_bar1[j] = sum(Var.*C[j,:]) + d[j];
        end


        # println(Y_bar2)

        Y = copy(Y_bar1)
        if any(abs(YR-Y_bar1) .< EPS)
            Y = YR + maximum(abs(Y_bar1-YR))
        end

        return Y_bar1,Var,Y,true
end
