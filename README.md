# MIMMP
This package need JuMP and CPLEX.
This Package solves minimum mixed integer bi-linear programming problems. This problem is formulated as follows:

$  min_{x,y} {y_1*y_2; y = Cx+d , Ax>=b, y>=0, x in {R^n_R, Z^n_Z, B^n_B} }  $,
where n = n_R + n_Z + n_B is the number of variables in different types. For more information refer to [1].

This package prvides a function called MIMMPSolver which solves above problem for given matrices and vectors A, b, C, d, VariableDef. VariableDef is a 3 element vector whose entries are the number of variables from each type real, integer, and binary, respectively. Meaning, the order of the variables must follow this order in the problem as wel. 

MIMMPSolver(A,b,C,d,VariableDef)

It also accepts the bounds of the x variable as key arguments in case the user wants to define them distinct from A matrix. (In this problem, variable y must be non-negative.)

MIMMPSolver(A,b,C,d,VariableDef,XLB=zeros(n),XUB=ones(n))

Other key arguments for this function are as follows:


			Explanation                 				default value

rEPS 			relative optimality gap (termination condition)		1e-6
-----------------------------------------------------------------------------------------------
aEPS 			absolute optimality gap (termination condition)		1e-6
-----------------------------------------------------------------------------------------------
TimeLim			time limit (termination condition)	 		2e+10
-----------------------------------------------------------------------------------------------
Solver		determines if the solver uses decision based algorithm 
				or criterion based algorithm			DB
			

