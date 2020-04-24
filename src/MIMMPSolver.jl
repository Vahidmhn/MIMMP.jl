function MIMMPSolver(A,b,C,d,VariableDef;Solver = "CB", XLB=[],XUB=[], rEPS = 1e-6,aEPS = 1e-6,TimeLim = 2e+10)
	VariableDef = Int.(VariableDef)
	if Solver == "CB"
		GUB, OptX = CBMIMMPSolver(A,b,C,d,VariableDef,XLB=XLB,XUB=XUB, rEPS = rEPS,aEPS = aEPS,TimeLim = TimeLim)
	else
		NODE,VARIABLE,ENHANCED_UB,MU,ENHANCED_CUT = 4,	3,	1,	1,1;
		GUB, OptX = DBMIMMPSolver(A,b,C,d,VariableDef,XLB=XLB,XUB=XUB,MU = MU,EnhancedUB = ENHANCED_UB,EnhancedCut = ENHANCED_CUT ,SelectVariableMode = VARIABLE,SelectNodeMode = NODE,TimeLim=TimeLim,rEPS = rEPS, aEPS = aEPS)
	end
	return GUB, OptX;
end
