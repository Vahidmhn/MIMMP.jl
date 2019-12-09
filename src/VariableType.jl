function VariableType(Def,i)

    if i <= Def[1]
        return :Cont
    elseif i > Def[1] & i <= Def[2] +Def[1]
        return :Int
    else
        return :Bin
    end
end
