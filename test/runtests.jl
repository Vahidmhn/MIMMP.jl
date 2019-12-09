using MIMMP
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here

C = trunc.(Int,readdlm("C.txt"));
A = trunc.(Int,readdlm("A.txt"));
b = trunc.(Int,readdlm("b.txt"));
d = trunc.(Int,readdlm("d.txt"));
VariableDef = [0 0 400];


@test MIMMPSolver(A,b,C,d,VariableDef,  XLB = zeros(400), XUB = ones(400), Solver = "CB")[1] == 6468
@test MIMMPSolver(A,b,C,d,VariableDef,  XLB = zeros(400), XUB = ones(400), Solver = "DB")[1] == 6468
