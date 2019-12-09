function GLBUpdater(Tree)

    GLB = Tree[1].LB;
    for i = 2:length(Tree)
        if Tree[i].LB < GLB
            GLB = Tree[i].LB
        end
    end
    return GLB
end
