*! version 1.0.0 Joe Long 05oct2018
pr car_sfe, eclass
	version 12.0
	syntax varlist [if] [in], Strata(varlist) [COVars(varlist) table]
	
	************
	*ROBUSTNESS*
	************
	marksample touse
	gettoken Y rest: varlist
	gettoken A leftover: rest
	
	*Check there aren't extra variables
	if "`leftover'" != "" {
		di as err "extraneous independent variables"
		exit 198
	}
	*Check strata is one variable 
	unab slist: `strata'
	if `:word count `slist'' != 1 {
		di as err "strata must be as one variable"
		exit 198
	}
	*Check variables are numeric
	cap confirm numeric variable `strata'
	if _rc != 0 {
		di as err "strata must be numeric"
		exit 198
	}
	cap confirm numeric variable `A'
	if _rc != 0 {
		di as err "treatment must be numeric"
		exit 198
	}
	qui cou if `A' == 0 & `touse'
	if r(N) == 0 {
		di as err "must have a control equal to zero"
		exit 198
	}
	
	*******************************************
	*GENERATE VARIABLES/LOCALS FOR CALCULATION*
	*******************************************
	qui levelsof `A' if `touse', l(treats)
	qui levelsof `strata' if `touse', l(groups)
	loc s = 0
	foreach strat of loc groups {
		*First check if the strata has any empty cells
		loc check_strat = 0 
		foreach a of loc treats {
			qui cou if `strata' == `strat' & `A' == `a' & `touse'
			if `r(N)' == 0 {
				loc ++check_strat
			}
		}
		*If strata has no empty cells, then save all indicator variables in strata to use
		if `check_strat' == 0 {
			loc ++s 
			foreach a of loc treats {
				if `a' != 0 {
					tempvar I_`a'_`s'
					gen `I_`a'_`s'' = `strata' == `strat' & `A' == `a' & `touse'
				}
			}
			tempvar I_`s'
			gen `I_`s'' = `strata' == `strat' & `touse'
			loc Scal `Scal' `I_`s''
		}
		*Otherwise, mark the strata as defective, ignore strata
		else {
			replace `touse' = 0 if `strata' == `strat' & `touse'
		}
	}
	/*
	foreach a of loc treats {
		if `a' != 0 {
			loc s = 0 
			foreach strat of loc groups {
				loc ++s
				loc interaction `interaction' `I_`a'_`s''
			}
			loc column `column' "`A'_`a'"
			tempvar A_`a'
			gen `A_`a'' = `A' == `a' & `touse'
			loc T `T' `A_`a''
		}
	}
	*/
	loc s = 0
	foreach strat of loc groups {
		loc ++s
		foreach a of loc treats {
			if `a' != 0 {
				loc interaction `interaction' `I_`a'_`s''
			}
		}
	}
	foreach a of loc treats {
		if `a' != 0 {
			loc column `column' "`A'_`a'"
			tempvar A_`a'
			gen `A_`a'' = `A' == `a' & `touse'
			loc T `T' `A_`a''
		}
	}
	
	*********
	*REGRESS*
	*********
	if "`table'" != "" {
		reg `Y' i.`A' i.`strata' `covars' if `touse', robust
	}
	tempname b V N df b2 V2 V_H V_hc mult
	
	*Check for perfect collinearity
	noi mata: collinear("`Scal'", "`interaction'", "`covars'", "`touse'", "`mult'")
	if `mult' != 0 {
			di as err "issue with collinear regressors"
			exit 198
	}
	noi mata: myregress("`Y'", "`Scal'", "`interaction'", "`touse'", "`b'", "`V'", "`N'", "`df'", "`V_H'", "`V_hc'", "`covars'", "`T'")

	********
	*OUTPUT*
	********
	mat coln `b' = `column'
	mat rown `b' = "`A'"
	mat coln `V' = `column'
	mat rown `V' = `column'
	
	ereturn post `b' `V', esample(`touse') dof(`=`df'')
	eret display
	eret scalar N = `N'
	eret matrix V_H = `V_H'
	eret matrix V_hc = `V_hc'
	eret local cmd "car_sfe"
	
end

*****************
*MATA REGRESSION*
*****************
mata:
	mata clear
	void myregress(string scalar depvar_s, string scalar strata_s, string scalar interact_s, ///
		string scalar touse_s, string scalar b_s, string scalar v_s, ///
		string scalar n_s, string scalar df_s, string scalar v_h, ///
		string scalar v_y, string scalar cov, string scalar treatment)
	{
		real vector y, Xpy, beta, beta_sfe, e2, N_s, T, I_S
		real matrix X, X_sfe, XpXi, XpeepX, B, V_H, R, vc, S_cal, V_hc
		real scalar k, n, S, A, as, c
		
		y 		= st_data(., depvar_s, touse_s)
		X 		= st_data(., (tokens(strata_s), tokens(interact_s), tokens(cov)), touse_s)
		S_cal 	= st_data(., strata_s, touse_s)
		n 		= rows(y)
		S 		= cols(S_cal)
		N_s		= colsum(S_cal)/n
		I_S 	= J(1, S, 1)
		k		= cols(X)
		as 		= cols(st_data(., (tokens(strata_s), tokens(interact_s)), touse_s))
		c 		= k - as
		
		Xpy 	= quadcross(X, y)
		XpXi 	= cholinv(quadcross(X, X))
		beta 	= XpXi * Xpy
		e2 		= (y - X*beta):^2
		XpeepX 	= quadcross(X, e2, X)
		vc 		= n * quadcross(quadcross(XpXi', XpeepX)', XpXi) * n/(n-k)

		beta 	= beta[1..as, 1]
		B 		= colshape(beta, S)[|2,1\.,.|]
		A 		= rows(B)
		B		= colshape(colshape(B, A)', S)
		T 		= B * N_s'
		V_H 	= quadcross((B - T*I_S)', N_s, (B - T*I_S)')
		R 		= J(A, S, 0)
		for (i=1; i<=S; i++) {
			R 	= R, diag(J(A, 1, N_s[i]))
		}
		R		= R, J(A, c, 0)
		V_hc 	= quadcross(quadcross(R', vc)', R')

		X_sfe 	= st_data(., (tokens(treatment), tokens(strata_s), tokens(cov)), touse_s)
		beta_sfe = cholinv(quadcross(X_sfe, X_sfe))*quadcross(X_sfe, y)

		st_matrix(b_s, beta_sfe[1..A,1]')
		st_matrix(v_s, (V_H + V_hc)/n)
		st_matrix(v_h, V_H/n)
		st_matrix(v_y, V_hc/n)
		st_numscalar(n_s, n)
		st_numscalar(df_s, n-k) 
	}	
end

***********************
*PERFECT COLLINEARITY?*
***********************

mata:
	mata clear
	void collinear(string scalar strata_s, string scalar interact_s, string scalar cov, string scalar touse_s, string scalar mult)
	{
		real matrix X
		real scalar rc
		
		X 		= st_data(., (tokens(strata_s), tokens(interact_s), tokens(cov)), touse_s)
		rc	 	= rank(X) < cols(X)
		
		st_numscalar(mult, rc)
	}
end
