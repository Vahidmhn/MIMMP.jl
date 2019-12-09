function Vert(PP :: Array{VerticeType,1}, NewCut,CutNum)
    const EPS = 0;
    V_Plus = [];
    v_Plus = 0;
    V_Minus = [];
    v_Minus = 0;
    for i = 1:length(PP)
        if (sum(NewCut[1:end-1].*PP[i].Point) + NewCut[end] > EPS)
            push!(V_Plus, i)
            v_Plus = v_Plus + 1;
        else
            push!(V_Minus, i)
            v_Minus = v_Minus + 1;
        end
    end
    # println(V_Plus)
    # println(V_Minus)
# ---------------------------------------------------------------------------
    # Step I
    if v_Minus == 0
        error("The cut makes it infeasible!!!!")
        return []
    elseif v_Plus == 0
        println("The cut is redundant!!!!")
##############################################################
        return []
##########################################################3##
        return PP
    end

    if v_Plus <= v_Minus
        V1 = V_Plus;
        V2 = V_Minus;
    else
        V2 = V_Plus;
        V1 = V_Minus;
    end
# ---------------------------------------------------------------------------
    # Step II

    VS = [];
    counter = 1;
    for u in V1
        for v in PP[u].adj
            if v in V2
                W = FindIntersectingPoint(PP[u].Point,PP[v].Point,NewCut)
                NewPoint = VerticeType([],[],[],-1)
                NewPoint.Point = W;
                if v_Plus <= v_Minus
                    PP[v].adj = filter(i -> (i != u), PP[v].adj);
                    push!(PP[v].adj,length(PP)+counter);
                    counter = counter + 1;
                    push!(NewPoint.adj,v);
                else
                    PP[u].adj = filter(i -> (i != v), PP[u].adj);
                    push!(PP[u].adj,length(PP)+counter)
                    counter = counter + 1;
                    push!(NewPoint.adj,u)
                end
                temp = intersect(PP[u].J,PP[v].J)
                for ii in temp
                    push!(NewPoint.J,ii)
                end
                push!(NewPoint.J,CutNum)
                if isempty(VS)
                    VS = [NewPoint];
                else
                    push!(VS,NewPoint);
                end
            end
        end
    end
# ---------------------------------------------------------------------------
    IndTrans = Dict();
    for i = 1:v_Minus
        IndTrans[V_Minus[i]]=i;
        for j = 1:length(VS)
            VS[j].adj[VS[j].adj .== V_Minus[i]] =i;
        end
    end

    for i in V_Minus
        for j = 1:length(PP[i].adj)
            if PP[i].adj[j] > length(PP)
                PP[i].adj[j] = PP[i].adj[j]- v_Plus;
            else
                PP[i].adj[j] = IndTrans[PP[i].adj[j]];
            end
        end
    end
    PP = [PP[i] for i in V_Minus];

    # # Step III
    for u = 1:length(VS), v = 1:length(VS)
        if length(intersect(VS[u].J,VS[v].J)) == length(PP[1].Point)-1
            if !(length(V_Minus)+v in VS[u].adj)
                push!(VS[u].adj,length(V_Minus)+v)
            end
            if !(length(V_Minus)+u in VS[v].adj)
                push!(VS[v].adj,length(V_Minus)+u)
            end
        end
    end

    for i = 1:length(VS)
        VS[i].T = prod(VS[i].Point);
        push!(PP,VS[i])
    end
    return PP;

end
