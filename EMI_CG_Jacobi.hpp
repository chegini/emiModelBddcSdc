#ifndef INTEGRATE_CG_JACOBI_HH
#define INTEGRATE_CG_JACOBI_HH

template <class Grid, class Functional, class VariableSet, class Spaces, class elementType, class Vector, class Options>
typename VariableSet::VariableSet semiImplicit_CG_Jacobi(	GridManager<Grid>& gridManager,
        	                                                Functional& F,
        	                                                VariableSet const& variableSet, 
        	                                                Spaces const& spaces,
        	                                                Grid const& grid, 
                                                          Options const& options,
        	                                                std::string out,
        	                                                bool cg_semi, 
        	                                                bool direct,
        	                                                typename VariableSet::VariableSet u,
        	                                                elementType & uAll,
                                                          Vector & sol_semi)
{

  double dt = options.dt;
  int maxSteps = std::floor(options.T/options.dt);
  int order  = options.order;
  std::cerr << "semiImplicitCGJacobi: dt = " << dt << " maxSteps " << maxSteps << "\n";

  using namespace boost::fusion;
  typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Functional> >  SemiLinearization;
	typedef VariationalFunctionalAssembler<SemiLinearization> Assembler;
  typedef typename Assembler::field_type field_type;
  typedef Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<field_type,1,1>> Matrix;

  int number_cells = gridManager.grid().size(0);

	auto& timer = Timings::instance();
	Assembler assembler(spaces);

	constexpr int nvars = Functional::AnsatzVars::noOfVariables;
	constexpr int neq = Functional::TestVars::noOfVariables;

  typedef typename Functional::OriginVars::template CoefficientVectorRepresentation<0,neq>::type LinearSpace;

  size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
  size_t  size = variableSet.degreesOfFreedom(0,nvars);
  std::cout << " size  "<< size <<"   nnz  " << nnz << std::endl;

  uAll = component<0>(u);

  if(options.plot) writeVTK(uAll,out+"/emiCGInitial"+paddedString(0,2),
                    IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

  SemiImplicitEulerStep<Functional>  eq(&F,dt);
  double const tau = dt;
  eq.setTau(tau);
  
  size_t nDofs = variableSet.degreesOfFreedom(0,nvars);

	auto du(u);
  auto step_test(u);
	du *= 0;

	Vector u_semi(nDofs);

  assembler.assemble(SemiLinearization(eq,u,u,du),options.assemblyThreads);
  AssembledGalerkinOperator<Assembler> A(assembler);
  Matrix LHS = assembler.template get<Matrix>(false); 
  for (int time_step=0; time_step<maxSteps; ++time_step) 
  {

    std::cout << "CG in semi-Implicit "<<time_step << std::endl;
    assembler.assemble(SemiLinearization(eq,u,u,du),Assembler::RHS,options.assemblyThreads);

    auto step = variableSet.zeroCoefficientVector();
    auto rhs = assembler.rhs();
 
    Vector rhs_temp(nDofs), step_temp(nDofs); 
    rhs.write(rhs_temp.begin());
    step.write(step_temp.begin());

    Vector rhs_vec(nDofs);
    rhs.write(rhs_vec.begin());

    // //---------------------------------------------------------------
    // // apcg
    // //---------------------------------------------------------------
    double cgTol = options.cgTol;
    double maxCGIter = options.maxCGIter;
    auto solver = [&cgTol,&maxCGIter](Matrix const& J, Vector& du, Vector& r) {
      if (maxCGIter>0) {

        DefaultDualPairing<Vector,Vector> dp;
        Dune::MatrixAdapter<Matrix,Vector,Vector> Jop(J);
        JacobiPreconditioner<Dune::MatrixAdapter<Matrix,Vector,Vector>> preco(Jop,1.0);
        PCGEnergyErrorTerminationCriterion<double> term(cgTol,maxCGIter);
        Pcg<Vector,Vector> pcg(Jop,preco,dp,term,0);
        Dune::InverseOperatorResult res;
        pcg.apply(du,r,res);

        // A*du = r
        // the solution(du) need to be shifted when using CG
        // shift s = argmin_s 1/2 (du + sI)^T A (du + sI) - (du + sI)^T r
        // then update du = du + sI
        // f(s) = argmin_s 1/2 I^T A I s^2 + du^T A I s - I^T r s + ... constant
        // f'(s)  = I^T A I s + du^T A I - I^T r = 0
        // s = (I^T r - du^T A I)/ (I^T A I)
        size_t size = du.size();
        using BlockVectorX = Dune::BlockVector<Dune::FieldVector<double,1> >;
        BlockVectorX I(size);
        BlockVectorX v(size);

        v *=0;
        double w =0;
        double p =0;
        double c =0;
        double s =0;
        I = 1;

        J.umv(I,v);

        w = I.dot(r);//dot(I,r);
        p = du.dot(v);//dot(du,v);
        c = I.dot(v);//dot(I,v);
  
        s = (w-p);
        s *= (1/c);
        I *=s;
        du+=I;
      } 
      else {
        try {
          DirectSolver<Vector,Vector> solver(J);
          solver.apply(du,r);
        }
        catch(...) {
          std::cerr << "solver failed on " << J.N() << "x" << J.M() << " matrix\n";
          std::cerr << J;
          std::cerr.flush();
          abort();
        }
      }//else 
    };

    solver(LHS,step_temp,rhs_vec);

    for (int j=0; j<nDofs; ++j) {
       at_c<0>(step_test.data).coefficients()[j] = step_temp[j];
    }
    
    // u = u + step
    component<0>(u) += component<0>(step_test);
    
    uAll = component<0>(u);
    
    if(options.plot) 
      writeVTK(uAll,out+"/emiCG"+paddedString(time_step,2),
                  IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

    writeVTK(uAll,out+"/emiCG"+paddedString(time_step,2),
                  IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
    sol_semi *= 0;
    step_test.write(sol_semi.begin());

  }

  writeVTK(uAll,out+"/emiCGLast",
                  IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

  writeVectorTofile(u,"matlab_dir/emiCGLast");
	return u;
}
#endif