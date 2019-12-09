function Dominated(Point, S,s)
    for i = 1:length(S)
        if i == s; continue; end
        if all(S[i].Point .<= Point)
            return true
        end
    end
    return false
end
