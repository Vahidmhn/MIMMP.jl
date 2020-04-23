function DBMIMMPSolver(A,b,C,d,VariableDef;XLB=[],XUB=[],MU = .9,EnhancedUB = 0,EnhancedCut = 0,SelectVariableMode = 1,SelectNodeMode = 0, TimeLim = 2e+10,
                rEPS = .000001, aEPS = .000001)
    ### Initialization
    t0 = time_ns();
    nVar = size(A)[2];
    if XLB == []
        XLB = -1e9*ones(nVar)
    end
    if XUB== []
        XUB = 1e9*ones(nVar)
    end

    # Algorithm's setting ##################################################

    # Initializing the Type 2 node selection
    Type2Iter = 0;
    Type2LastLB = 0;
    Type2Flag = 0;
    Type2MaxIter = 1000;
    Type2MaxImprovement = 0.05;

    # Initializing the Type 3 or 4 variable selection
    if (SelectVariableMode == 3) || (SelectVariableMode == 4)
        TempSelectVariableMode = SelectVariableMode
        SelectVariableMode = 1
    end
    # Initializing the Type 3 or 4 node selection
    if (SelectNodeMode == 3) || (SelectNodeMode == 4)
        TempSelectNodeMode = SelectNodeMode
        SelectNodeMode = 0
    end
    ########################################################################

    const  inf = 9999999;
    const  EPS = 0.000001
    nObj = size(C)[1];
    nVar = size(A)[2];
    nCons = length(b);

    #  Initial Global Bounds
    GLB = 0;
    GUB = inf;

    # Initialization of Temp Variables
    NodeCounter = 1;
    RO1 = zeros(sum(VariableDef[2:3]));
    n1 = zeros(sum(VariableDef[2:3]));
    RO2 = zeros(sum(VariableDef[2:3]));
    n2 = zeros(sum(VariableDef[2:3]));



    ## Main Algorithm
    Tree = []
    ## Preprocessing
    Cuts = [];
    if EnhancedUB == 1 | EnhancedCut == 1
        nVars = size(C,2)
        Xs, Ys, GUB, Cuts = Preprocess2(A,b,C,d,XLB,XUB,VariableDef,TimeLim = .01*nVars)
        if Xs != []
            OptX = Xs;
        end
        if EnhancedUB == 0
            UB = inf;
        end
        if EnhancedCut == 0
            Cut = []
        end
    else
        OptX = [];
    end
    # ---------------------------------------------
    Time = (time_ns() - t0)/1e9
    if TimeLim-Time < 0.001;
        println("Time Limit reached! The solution is not a proved optimal")
        return GUB, OptX
    end
    # ---------------------------------------------

    GLB, X_Tilde,NewUB,NewX = ShaoRSolver(A,b,C,d,XLB,XUB,VariableDef,CUTs = Cuts,TimeLim = TimeLim - Time);

    if GLB == -2
        Time = (time_ns() - t0)/1e9
        println("Something went wrong, the problem seems to be infeasible or has numerical issues!")
        return -2, []
    end
    iter = 1;
    if isempty(X_Tilde)
        if GUB > NewUB
            GUB = NewUB
            OptX = NewX
            GLB = GUB
            return GUB, OptX
        else
            println("The problem is infeasible!")
            return GUB, []
        end
    elseif VariableDef[1] == nVar
        OptX = X_Tilde;
        Time = (time_ns() - t0)/1e9
    elseif all(isinteger.(X_Tilde[VariableDef[1]+1:end]))
        # We get the optimal solution
        println("The relaxed solution was optimal")
        OptX = X_Tilde;
        GUB = GLB
        Time = (time_ns() - t0)/1e9
    else
        delta1 = (ceil.(X_Tilde[VariableDef[1]+1:end])-X_Tilde[VariableDef[1]+1:end]).*RO1./(EPS+n1)
        delta2 = (X_Tilde[VariableDef[1]+1:end]-floor.(X_Tilde[VariableDef[1]+1:end])).*RO2./(EPS+n2)
        SelectedVariable = SelectVariable(VariableDef,X_Tilde,delta1,delta2,MU,
        n1 = n1, n2 = n2, TAV =100, ZZ = GLB ,xlb = XLB,xub = XUB, Type = SelectVariableMode)

        # Node 1
        NewXUB = copy(XUB);
        NewXUB[SelectedVariable] = floor(X_Tilde[SelectedVariable])
        # Time commands
        # ---------------------------------------------
        Time = (time_ns() - t0)/1e9
        if TimeLim-Time < 0.001;
            println("Time Limit reached! The solution is not a proved optimal")
            return GUB, OptX
        end
        # ---------------------------------------------
        Z, X_Tilde1,NewUB = ShaoRSolver(A,b,C,d,XLB,NewXUB,VariableDef,CUTs = Cuts,TimeLim = TimeLim - Time);
        if GUB > NewUB
            GUB = NewUB
        end
        # Time commands
        # ---------------------------------------------
        Time = (time_ns() - t0)/1e9
        if TimeLim-Time < 0.001;
            println("Time Limit reached! The solution is not a proved optimal")
            return GUB, OptX
        end
        # ---------------------------------------------
        if Z == -2
            println("Something went wrong, the problem seems to be infeasible or has numerical issues!")
            return -2, []
        end
        if !isempty(Z) && Z < GUB
            if all(isinteger.(X_Tilde1[VariableDef[1]+1:end]))
                GUB = Z;
                OptX = copy(X_Tilde1);
            else
                Tree = [Node(copy(X_Tilde1),copy(XLB),NewXUB,Z,GUB)];
                RO2[SelectedVariable-VariableDef[1]] = RO2[SelectedVariable-VariableDef[1]] + (Z-GLB)/(X_Tilde[SelectedVariable] - NewXUB[SelectedVariable])
                n2[SelectedVariable-VariableDef[1]] = n2[SelectedVariable-VariableDef[1]] + 1
            end
        end
        NodeCounter = NodeCounter + 1;

        # Node 2

        NewXLB = copy(XLB);
        NewXLB[SelectedVariable] = ceil(X_Tilde[SelectedVariable])
        # Time commands
        # ---------------------------------------------
        Time = (time_ns() - t0)/1e9
        if TimeLim-Time < 0.001;
            println("Time Limit reached! The solution is not a proved optimal")
            return GUB, OptX
        end
        # ---------------------------------------------
        Z, X_Tilde2,NewUB = ShaoRSolver(A,b,C,d,NewXLB,XUB,VariableDef,CUTs = Cuts,TimeLim = TimeLim - Time);
        if GUB > NewUB
            GUB = NewUB
        end
        # Time commands
        # ---------------------------------------------
        Time = (time_ns() - t0)/1e9
        if TimeLim-Time < 0.001;
            println("Time Limit reached! The solution is not a proved optimal")
            return GUB, OptX
        end
        # ---------------------------------------------
        if Z == -2
            println("Something went wrong, the problem seems to be infeasible or has numerical issues!")
            return -2, []
        end
        if !isempty(Z) && Z < GUB
            if all(isinteger.(X_Tilde2[VariableDef[1]+1:end]))
                GUB = Z;
                OptX = copy(X_Tilde2);
            else
                push!(Tree,Node(copy(X_Tilde2),NewXLB,copy(XUB),Z,GUB))
                RO1[SelectedVariable-VariableDef[1]] = RO1[SelectedVariable-VariableDef[1]] + (Z-GLB)/(NewXLB[SelectedVariable] - X_Tilde[SelectedVariable] )
                n1[SelectedVariable-VariableDef[1]] = n1[SelectedVariable-VariableDef[1]] + 1
            end
        end
        NodeCounter = NodeCounter + 1;
        if !isempty(Tree)
            GLB = GLBUpdater(Tree)
        end
        # Main Loop
        while (length(Tree) > 0) && (abs(GUB-GLB) >= aEPS) && ((GUB-GLB)/GUB >= rEPS)
            # Type 2 Node Selection modification
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            Type2Iter = Type2Iter + 1;
            if Type2Iter > Type2MaxIter
                Type2Iter = 0;
                Type2LastLB = GLB;
                Type2Flag = 1 - Type2Flag;
            elseif Type2Flag == 1 && (abs(Type2LastLB-GLB)/Type2LastLB >Type2MaxImprovement)
                Type2Flag = 0;
                Type2Iter = 0;
            end
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Type 3 or 4 Variable Selection modification
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (NodeCounter > 10) && ((SelectVariableMode == 3) || (SelectVariableMode == 4))
                SelectVariableMode = TempSelectVariableMode;
            end
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Type 3 or 4 Variable Selection modification
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (NodeCounter > 10) && ((SelectNodeMode == 3) || (SelectNodeMode == 4))
                SelectNodeMode = TempSelectNodeMode
            end

            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Select node to branch
            SelectedNode = SelectNode(Tree,VariableDef = VariableDef, MU = MU, RO1 = RO1,
                                     RO2 = RO2, Type = SelectNodeMode, Type2Flag = Type2Flag);
            node = deepcopy(Tree[SelectedNode])
            deleteat!(Tree,SelectedNode);

            # Select variable to branch
            delta1 = (ceil.(node.X[VariableDef[1]+1:end])-node.X[VariableDef[1]+1:end]).*RO1./(EPS+n1)
            delta2 = (node.X[VariableDef[1]+1:end]-floor.(node.X[VariableDef[1]+1:end])).*RO2./(EPS+n2)
            SelectedVariable = SelectVariable(VariableDef,node.X,delta1,delta2,MU,n1 = n1, n2 = n2,
                                    TAV =100,ZZ = node.LB, xlb = node.XLB, xub = node.XUB, Type = SelectVariableMode)

            # Node 1,TimeLim = TimeLim - Time
            NewXUB = copy(node.XUB);
            NewXUB[SelectedVariable] = floor(node.X[SelectedVariable])
            # Time commands
            # ---------------------------------------------
            Time = (time_ns() - t0)/1e9
            if TimeLim-Time < 0.001;
                println("Time Limit reached! The solution is not a proved optimal")
                return GUB, OptX
            end
            # ---------------------------------------------
            Z, X_Tilde,NewUB = ShaoRSolver(A,b,C,d,copy(node.XLB),NewXUB,VariableDef,CUTs = Cuts,TimeLim = TimeLim - Time);
            if GUB > NewUB
                GUB = NewUB
            end
            # Time commands
            # ---------------------------------------------
            Time = (time_ns() - t0)/1e9
            if TimeLim-Time < 0.001;
                println("Time Limit reached! The solution is not a proved optimal")
                return GUB, OptX
            end
            # ---------------------------------------------
            if Z == -2
                println("Something went wrong, the problem seems to be infeasible or has numerical issues!")
                return -2, []
            end
            if !isempty(Z) && Z < GUB
                if all(isinteger.(X_Tilde[VariableDef[1]+1:end]))
                    GUB = Z;
                    OptX = copy(X_Tilde);
                    # Type 2 Node Selection modification
                    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    if Type2Flag == 0
                        Type2Flag = 1;
                        Type2Iter = 0;
                        Type2LastLB = GLB;
                    end
                    # Type 2 Node Selection modification
                    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                else
                    push!(Tree,Node(copy(X_Tilde),copy(node.XLB),NewXUB,Z,GUB));
                    RO2[SelectedVariable-VariableDef[1]] = RO2[SelectedVariable-VariableDef[1]] +
                        (Z-node.LB)/(node.X[SelectedVariable] - NewXUB[SelectedVariable])
                    n2[SelectedVariable-VariableDef[1]] = n2[SelectedVariable-VariableDef[1]] + 1
                end
            end
            NodeCounter = NodeCounter + 1;
            # Node 2
            NewXLB = copy(node.XLB);
            NewXLB[SelectedVariable] = ceil(node.X[SelectedVariable])
            # Time commands
            # ---------------------------------------------
            Time = (time_ns() - t0)/1e9
            if TimeLim-Time < 0.001;
                println("Time Limit reached! The solution is not a proved optimal")
                return GUB, OptX
            end
            # ---------------------------------------------
            Z, X_Tilde,NewUB = ShaoRSolver(A,b,C,d,NewXLB,copy(node.XUB),VariableDef,CUTs = Cuts,TimeLim = TimeLim - Time);
            if GUB > NewUB
                GUB = NewUB
            end
            # Time commands
            # ---------------------------------------------
            Time = (time_ns() - t0)/1e9
            if TimeLim-Time < 0.001;
                println("Time Limit reached! The solution is not a proved optimal")
                return GUB, OptX
            end
            # ---------------------------------------------
            if Z == -2
                println("Something went wrong, the problem seems to be infeasible or has numerical issues!")
                return -2, []
            end
            if !isempty(Z) && Z < GUB
                if all(isinteger.(X_Tilde[VariableDef[1]+1:end]))
                    GUB = Z;
                    OptX = copy(X_Tilde);
                    # Type 2 Node Selection modification
                    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    if Type2Flag == 0
                        Type2Flag = 1;
                        Type2Iter = 0;
                        Type2LastLB = GLB;
                    end
                    # Type 2 Node Selection modification
                    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                else
                    push!(Tree,Node(copy(X_Tilde),NewXLB,copy(node.XUB),Z,GUB))
                    RO1[SelectedVariable-VariableDef[1]] = RO1[SelectedVariable-VariableDef[1]] +
                        (Z-node.LB)/(NewXLB[SelectedVariable] - node.X[SelectedVariable] )
                    n1[SelectedVariable-VariableDef[1]] = n1[SelectedVariable-VariableDef[1]] + 1
                end
            end
            NodeCounter = NodeCounter + 1;

            if !isempty(Tree)
                GLB = GLBUpdater(Tree)
            end
            # println("Iteration = $iter  |Time = $Time| Node = $NodeCounter  |  GLB = $GLB  |   GUB =  $GUB" )
            iter = iter + 1;

        end

    end

    Time = (time_ns() - t0)/1e9
    return GUB, OptX
end
