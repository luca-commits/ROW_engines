#ifndef BIGSTAB_HPP
#define BIGSTAB_HPP

#include <iostream> 
#include <Eigen/Dense>
#include <Eigen/Sparse>


//Matrix im CRS format, Basis 1 --> double* Matrix, int* ia, int* ja
//Vorkonditionierer: PardisoCockpit<double>* Precond
//Rechte seite: double* const b
//Lˆsung: double* const x
//Genauigkeit: double& eps
//Max iterations: unsigned& nsteps
//Minimal Anzahl iterations: unsigned min_no_steps


template <typename Preconditioner>
unsigned BiCGstab(Eigen::SparseMatrix<double>  Matrix, unsigned N,  Preconditioner& precond, const Eigen::VectorXd& b, Eigen::VectorXd & x, double eps, unsigned nsteps, unsigned min_no_steps)
{
	double D_PREC = 1e-16;
	bool verbose = 0;
	double bestRes = 1.0e100;
	double EuklidRes = 1.0e100;
	double bestEuklidRes = 1.0e100;
	bool   Euklid_Output = 0;

	Eigen::VectorXd v, p, phat, s, shat, t, r, r0;
	Eigen::VectorXd bestAbsSol(N);
	bool copyBestSolution = true;



	if(verbose) std::cout << "Start BiCGstab..." << std::endl;
	//const unsigned min_no_steps=5;

	Eigen::VectorXd precond_rhs; 

	double resid;

	// normb = |b|
	double nrmb = b.norm();
	if (nrmb < D_PREC) nrmb = 1.;
	r0 = b - Matrix * x;

	r = r0;

	resid = r.norm() / nrmb;

	if (resid < eps && min_no_steps == 0) {
		eps = resid;
		nsteps = 0;
		if (verbose) std::cout << "Leave BicG with error 0." << std::endl;
		if (verbose) std::cout << "I return the solution with the best iterative residual = " << eps << std::endl;
		return 0;
	}

	double alpha = 0, omega = 1, rho2 = 1;

	for (unsigned l = 1; l <= nsteps; ++l) {

		// rho1 = r0 * r
		//cout << l << " of " << nsteps << " Steps." << endl;
		const double rho1 = r0.dot(r);
		if (abs(rho1) < D_PREC) {

			if (l > 1 && min_no_steps > 0 && resid < 1.0e-15)
			{
				eps =r.norm() / nrmb;
				if (verbose) std::cout << (resid);
				if (verbose)  std::cout << "Assume that the BiCGStab methods converged to a precise enough solution." <<  std::endl;
				return 0;
			}


			if (b.norm() < 1.0e-17 && x.norm() < 1.0e-17)
			{
				eps = r.norm() / nrmb;
				if (verbose) std::cout << (b.norm());
				if (verbose)  std::cout << "Found zero RHS and zero initial guess. Converged in 1 Step." <<  std::endl;
				return 0;
			}

		}

		if (l == 1)
			p = r;                // p = r
		else {
			p = -1 * omega * v + p;       // p = (p-omega v)*beta+r  ------------------------- doubt 
			const double beta = rho1 * alpha / (rho2 * omega);
			p = p * beta;
			p = r + p;
		}
		
		// phat = C p
		phat = precond.solve(p);

		// v = A phat
		v = Matrix * phat;
		alpha = rho1 / r0.dot(v);

		// s = r - alpha v
		s = r - alpha * v;

		resid = s.norm() / nrmb;

		if (verbose) std::cout << "Step " << l << ", resid=" << resid <<  std::endl;

		//Catch divergence if it occurs
		if ((l > 20 && resid > 10.0) || (l > 40 && resid > 1.0) || (l > 60 && resid > 0.1) || (l > 80 && resid > 0.01) && l > min_no_steps)
		{
			// x += alpha phat
			x = x + alpha * phat; 
			eps = resid;
			nsteps = l;
			if (Euklid_Output)
			{
				//Compute EuklidRes
				r = b - Matrix * x;  //r=-(-b)-Ax=b-Ax
				EuklidRes =r.norm() / nrmb;  //EuklidRes = r/b
				if (verbose) std::cout << (EuklidRes);
				if (verbose)  std::cout << "Exit Point 0 " <<  std::endl;
			}
			if (resid > bestRes && copyBestSolution)
			{
				x = bestAbsSol;
				eps = bestRes;
			}
			if (verbose) std::cout << "No hope for further convergence. Leave BicG with error 1." << std::endl;
			if (verbose) std::cout << "I return the solution with the best iterative residual =" << eps << std::endl;
			return 1;
		}


		if (resid<eps && l>min_no_steps) {
			// x += alpha phat
			x += alpha * phat;
			eps = resid;
			nsteps = l;
			if (Euklid_Output)
			{
				//Compute EuklidRes
				r = b - Matrix * x;
				EuklidRes = r.norm() / nrmb;  //EuklidRes = r/b
				if (verbose) std::cout << EuklidRes << std::endl;
				if (verbose)  std::cout << "Exit Point 1 " <<  std::endl;
			}
			if (resid > bestRes&& copyBestSolution)
			{
				x = bestAbsSol;
				eps = bestRes;
			}
			if (verbose) std::cout << "Leave BicG with error 0." << std::endl;
			if (verbose) std::cout << "I return the solution with the best iterative residual =" << eps << std::endl;
			return 0;
		}

		shat = precond.solve(s);

		t = Matrix * shat;

		// omega = t*s / t*t
		omega = t.dot(s) / t.dot(t);

		// x += alpha phat + omega shat
		x += alpha * phat + omega * shat;

		// r = s - omega t
		r = s - omega * t;

		rho2 = rho1;

		resid = r.norm() / nrmb;
		if (resid < bestRes)
		{
			bestRes = resid;
			bestAbsSol = x;

		}
		if (Euklid_Output)
		{
			//Compute EuklidRes
			r = b - Matrix * x;
			EuklidRes =  r.norm() / nrmb;  //EuklidRes = r/b
			if (verbose) std::cout << (EuklidRes);
			if (verbose)  std::cout <<  std::endl;
		}

		if (resid<eps && l>min_no_steps) {
			eps = resid;
			nsteps = l;
			if (Euklid_Output)
			{
				//Compute EuklidRes
				r = b - Matrix * x;
				EuklidRes = r.norm() / nrmb;  //EuklidRes = r/b
				if (verbose) std::cout << (EuklidRes)  << std::endl;
				if (verbose)  std::cout << "Exit Point 2 " <<  std::endl;
			}
			if (resid > bestRes&& copyBestSolution)
			{
				x = bestAbsSol;
				eps = bestRes;
			}

			if (verbose)  std::cout << "Step " << l << ", resid=" << resid <<  std::endl;
			if (verbose) std::cout << "Leave BicG with error 0." << std::endl;
			if (verbose) std::cout << "I return the solution with the best iterative residual =" << eps << std::endl;
			return 0;
		}

		if (std::abs(omega) < D_PREC) {
			eps = resid;
			if (Euklid_Output)
			{
				//Compute EuklidRes
				r = b - Matrix * x;
				EuklidRes = r.norm() / nrmb;  //EuklidRes = r/b
				std::cout << (EuklidRes) << std::endl;
				std::cout << "Exit Point 3 " << std::endl;
			}
			if (resid > bestRes&& copyBestSolution)
			{
				x = bestAbsSol;
				eps = bestRes;
			}
			if (verbose) std::cout << "Leave BicG with error 3." << std::endl;
			if (verbose) std::cout << "I return the solution with the best iterative residual =" << eps << std::endl;
			return 3;
		}
	}//end l loop
	eps = resid;
	if (Euklid_Output && verbose)
	{
		//Compute EuklidRes
		r = b - Matrix * x;
		EuklidRes = r.norm() / nrmb;  //EuklidRes = r/b
		if (verbose) std::cout << (EuklidRes);
		if (verbose)  std::cout << "Exit Point 4 " <<  std::endl;
	}
	if (resid > bestRes&& copyBestSolution)
	{
		x = bestAbsSol;
		eps = bestRes;
	}
	if (verbose) std::cout << "Leave BiCG with error 1." << std::endl;
	if (verbose) std::cout << "I return the solution with the best iterative residual =" << eps << std::endl;
	return 1;
}

#endif