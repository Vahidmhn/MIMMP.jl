mutable struct NodeType2
    Point :: Array{Float64,1}
    LB :: Float64
    UB :: Float64
end
mutable struct Node
    X :: Array{Float64,1}
    XLB :: Array{Float64}
    XUB :: Array{Float64}
    LB :: Float64
    UB :: Float64
end
mutable struct CutType
    Coeff :: Array{Float64,1}
    RHS :: Float64
end
mutable struct VerticeType
    Point :: Array{Float64,1}
    adj ::  Array{Int64,1}
    J :: Array{Int64,1}
    T :: Float64
end
