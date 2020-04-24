function CBMIMMPSolver(A,b,C,d,VariableDef;XLB=[],XUB=[], rEPS = 1e-6,aEPS = 1e-6,TimeLim = 2e+10)
    t0 = time_ns();
    # Initialization
    const EPS = 0.000001;
    const inf = 9999999;

    nObj = size(C)[1];
    nVar = size(A)[2];
    nCons = length(b);

    if XLB == []
        XLB = -1e9*ones(nVar)
    end
    if XUB== []
        XUB = 1e9*ones(nVar)
    end

    NodeCounter = 1;
    GUB = Inf;
    GLB = 0;
    # Node 0
    OptX = zeros(nVar)
    YI = zeros(nObj);
    Var = zeros(nVar);
    PayoffTable = zeros(nObj,nObj);
    for i = 1:nObj
        ObjOrd = [j for j = i:i+nObj-1]
        ObjOrd[ObjOrd .> nObj] = ObjOrd[ObjOrd .> nObj]-nObj
        ObjVal, Var = ArgLexMin(A,b,C,d,XLB,XUB,VariableDef,ObjOrd; TimeLimit = TimeLim-(time_ns() - t0)/1e9)

        if (TimeLim-(time_ns() - t0)/1e9) < 0.001;
            println("Time Limit reached! The solution is not a proved optimal")
            return GUB, OptX;
        end
        YI[i] = ObjVal[i]
        PayoffTable[i,:] = ObjVal
        if GUB >= prod(PayoffTable[i,:])
            GUB = prod(PayoffTable[i,:])
            OptX = copy(Var)
        end
    end
    GLB = prod(YI);
    S = [NodeType2(YI,GLB,GUB)]
    Y = zeros(nObj)
    # Main Loop
    s = MinimalT1(S)[1];
    while (length(S) > 0) && (abs(GUB-GLB) >= aEPS) && ((GUB-GLB)/GUB >= rEPS)

        # Time commands ---------------------------------------------------
        if (TimeLim-(time_ns() - t0)/1e9) < 0.001;
            println("Time Limit reached! The solution is not a proved optimal")
            return GUB, OptX;
        end
        # -----------------------------------------------------------------
        # Explore node s
        Y_bar,Var,Y,IsExplored = Explore(A,b,C,d,XLB,XUB,VariableDef,S[s].Point)

        NodeCounter += 1

        if IsExplored
            if (prod(Y_bar) < GUB)
                GUB = prod(Y_bar)
                OptX = copy(Var)
            end

            # Branch
            for i = 1:nObj
                Point = copy(S[s].Point)
                Point[i] = Y[i];
                if (GUB - prod(Point) > aEPS) & ((GUB-prod(Point))/GUB >= rEPS) & !(Dominated(Point,S,s))
                    push!(S,NodeType2(Point,prod(Point),GUB))
                end
            end
        end
        deleteat!(S,s);
        if ~((length(S) > 0) && (abs(GUB-GLB) >= aEPS) && ((GUB-GLB)/GUB >= rEPS)); break; end
        s,GLB = MinimalT1(S)

    end

    return GUB, OptX;
end
