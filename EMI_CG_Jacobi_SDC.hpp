/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/en/numerik/software/kaskade-7.html               */
/*                                                                           */
/*  Copyright (C) 2011-2014 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef INTEGRATE_CG_JACOBI_SDC_HH
#define INTEGRATE_CG_JACOBI_SDC_HH

// ./emiModel --refine 5 --stol 1e-4 --T 0.02 --sdc_contraction 0.1

using namespace Kaskade;

template <class Matrix, class Vectors, class ReactionDerivatives, class Solver>
typename Matrix::field_type sdcIterationStepJacobi(SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix const& Shat, Solver const& solve, 
                                              Matrix const& M, Matrix const& Stiff,
                                              Vectors const& rUi, ReactionDerivatives const& rDu, Vectors const& Mdu, Vectors& du)
{
  auto const& pts = grid.points();
  int const n = pts.size()-1;      // number of subintervals

  assert(Mdu.size()>=n); 
  assert(du.size()>=n+1);  // including start point
  

  size_t const dofs = Mdu[0].size();

  // compute exact integration matrix
  auto const& S = grid.integrationMatrix();


  // initialize correction at starting point to zero
  du[0] = 0.0; 

  typedef typename Vectors::value_type Vector;
  Vector rhs(dofs), tmp(dofs);  // declare here to prevent frequent reallocation
  
  // perform n Euler steps
  typename Matrix::field_type norm = 0;
  Matrix J = M;
  for (int i=1; i<=n; i++)
  {
    // matrix J = M - Shat_i-1,i*(A+f_u)
    for (size_t row=0; row<J.N(); ++row)
    {
      auto colJ = J[row].begin(); 
      auto end = J[row].end();
      auto colM = M[row].begin();
      auto colA = Stiff[row].begin();
      auto colR = rDu[i][row].begin();
      auto endR = rDu[i][row].end();
      
      while (colJ != end)
      {
        *colJ = *colM - Shat[i-1][i] * *colA;
        
        if (colR != endR && colJ.index() == colR.index()) // f_u can have subset of sparsity pattern
        {
          *colJ -= std::min(0.5* *colM, Shat[i-1][i] * *colR); // guarantee M - Shat*fu is nonnegative -- reduce by at most 50%
          ++colR;
        }
        ++colJ; ++colM; ++colA;           
      }
    }
    

    //  right-hand side for linear system

    // M * ( u_i^{k} - u_{i+1}^k + du_i)
    rhs = Mdu[i-1];
    M.umv(du[i-1],rhs);

    // add sum_j S_ij r_j to right hand side
    for (int j=0; j<=n; ++j)
      rhs.axpy(S[i-1][j],rUi[j]);
       
    // add sum_j Shat_ij r'(u_j) du_j with r' = A + f_u
    tmp = 0;
    for (int j=0; j<i; ++j) // TODO: start at 1 instead of 0? du[0] is zero anyway...
    {
      rDu[j].usmv(Shat[i-1][j],du[j],rhs);
      tmp.axpy(Shat[i-1][j],du[j]);
    }
    Stiff.umv(tmp,rhs);
      
    // solve linear system
    du[i] = du[i-1]; // previous increment is probably a good starting value
    solve(J,du[i],rhs);
    
    // evaluate norm of correction
    norm += (pts[i]-pts[i-1]) * (du[i]*rhs);
  }   // end i - loop
  
  return std::sqrt(norm/(pts[n]-pts[0]));
} 

template <class State, class StateUe, class TimeGrid, class Vector, class Eq,class EqSemi, class Assem, class elementType, class CellFilter, class Options>
void computeRHS(int step, 
                State const& x,
                CellFilter & Cellfltr,
                int number_cells,
                std::vector<StateUe> const& collocationUe,
                size_t size_e, 
                TimeGrid const& grid, 
                int sweep, double t, std::vector<size_t> const& expandedIndices,
                Eq& eq, EqSemi& eqSemi, Assem& assembler, Options const& options, std::vector<Vector>& ru,
                elementType & uAll, double dt,
                std::chrono::duration<double, std::milli> & ms_double_ass)
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  using namespace boost::fusion;
  typedef typename Eq::OriginVars::template CoefficientVectorRepresentation<0,1>::type CoefficientVectorsU;  
  typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Eq> >  SemiLinearization;
  size_t const nrDofs = expandedIndices.size();
  
  State stateTmp(x);
  State stateTmp0(x);
  State dstateTmp(x);

  dstateTmp = 0;
  eqSemi.setTau(dt);
  eq.Mass_stiff(0);
  auto const& pts = grid.points();
  for (int ii=0; ii<pts.N(); ii++)
  {
    // In the first sweep, the iterates and hence the right hand sides coincide. Exceptions are 
    // - nonautonomous rhs (but this is rare here and will be otherwise corrected in the next sweeps)
    if (sweep==0 && ii>0) 
    {                    
      ru[ii] = ru[ii-1]; 
    }
    else 
    {
      eq.time(t+pts[ii]-pts[0]);
      at_c<0>(stateTmp.data) = collocationUe[ii];

      auto recordTime_ass1 = high_resolution_clock::now();
      assembler.template assemble<AssemblyDetail::TakeAllBlocks,CellFilter>(SemiLinearization(eqSemi,stateTmp,stateTmp,dstateTmp),Cellfltr,Assembler::RHS,options.assemblyThreads);
      auto recordTime_ass2 = high_resolution_clock::now();

      duration<double, std::milli> ms_double_ass_temp = recordTime_ass2 - recordTime_ass1;
      ms_double_ass +=ms_double_ass_temp;

      // extract restricted subgrid rhs values
      CoefficientVectorsU r(assembler.rhs()); 
      ru[ii] = Vector(nrDofs);
      for (size_t j=0; j<nrDofs; ++j){
        size_t ej = expandedIndices[j];
        ru[ii][j] = at_c<0>(r.data)[ej]*(1/dt);
      }
    }
  }
}


template <class Grid, class Equation, class VariableSet, class Spaces, class elementType, class CellFilter, class Options, class OptionStatistics>
typename VariableSet::VariableSet semiImplicit_CG_Jacobi_SDC( GridManager<Grid>& gridManager,
                                                              Equation& eq,
                                                              CellFilter & Cellfltr,  VariableSet const& variableSet, Spaces const& spaces,
                                                              Grid const& grid,
                                                              typename VariableSet::VariableSet x,
                                                              std::vector<std::set<int>> index2Cells_new,
                                                              Options const& options, OptionStatistics& statistics,
                                                              std::string output,
                                                              elementType & uAll,
                                                              std::vector<std::set<int>> index2IndexsSet
                                                              )
{
  std::vector<int> cellsMarked;
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  using namespace boost::fusion;
  std::setprecision(16);

  // --------------------------------------------------------------------------------------------
  // timers
  // --------------------------------------------------------------------------------------------
  boost::timer::cpu_timer assemblyTimer;
  assemblyTimer.stop();

  boost::timer::cpu_timer assemblyRhsTimer;
  assemblyRhsTimer.stop();
  boost::timer::cpu_timer assemblyReactionTimer;
  assemblyReactionTimer.stop();
  
  boost::timer::cpu_timer refineTimer;
  refineTimer.stop();
  boost::timer::cpu_timer coarsenTimer;
  coarsenTimer.stop();
  boost::timer::cpu_timer sdcTimer;
  sdcTimer.stop();
  boost::timer::cpu_timer sdcOdeTimer;
  sdcOdeTimer.stop();
  boost::timer::cpu_timer outputTimer;
  outputTimer.stop();
  boost::timer::cpu_timer MaxStepsTimer;
  MaxStepsTimer.stop();

  if(options.plot) printuAll(x,uAll,options.order, output+"/initialValue-"+paddedString(0),"potential");
  int number_cells = gridManager.grid().size(0);
  // --------------------------------------------------------------------------------------------
  // typedef: reference to its elements: states & CoefficientVectors
  // --------------------------------------------------------------------------------------------
  typedef typename Equation::OriginVars::VariableSet State;
  typedef typename boost::fusion::result_of::value_at_c<typename State::Sequence,0>::type StateUe;
  
  typedef typename Equation::OriginVars::template CoefficientVectorRepresentation<0,1>::type CoefficientVectorsU;
  // --------------------------------------------------------------------------------------------
  // assembler
  // --------------------------------------------------------------------------------------------

  typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Equation> >  SemiLinearization;
  typedef VariationalFunctionalAssembler<SemiLinearization> Assembler;
  using BlockVectorX = Dune::BlockVector<Dune::FieldVector<double,1> >;

  int const dim = Grid::dimension;
    
  Assembler assembler(spaces);

  // --------------------------------------------------------------------------------------------
  // tolerance X
  // --------------------------------------------------------------------------------------------
  std::vector<std::pair<double,double> > tolX(variableSet.noOfVariables);

  // adaptivity for transmembrane voltage only -> large tolerances for the remaining variables
  for (int i=0; i<tolX.size(); ++i) {
    tolX[i] = std::make_pair(10000,10000);
  }
  tolX[0] = std::make_pair(0,0);

  auto tolXC = tolX;
  tolXC[0] = std::make_pair(options.aTol,0);
  for (int i=1; i<tolXC.size(); ++i) {
    tolXC[i].first *= 0.7;
    tolXC[i].second *= 0.7;
  }
    
  // --------------------------------------------------------------------------------------------
  // define T, dt
  // --------------------------------------------------------------------------------------------
  double T = options.T;
  double dt = options.dt;
  int maxSteps = std::floor(T/options.dt);  
  std::cerr << "dt = " << dt <<" maxSteps " << maxSteps << "\n";
  
  // --------------------------------------------------------------------------------------------
  // initilize statistics: maxLevel, totalSweeps, what about sdcTime? avgDofs?
  // --------------------------------------------------------------------------------------------
  size_t Total_dofCount = 0;
  size_t dofCount = 0;
  statistics.maxLevel = 0;
  statistics.totalSweeps = 0;
 
  // --------------------------------------------------------------------------------------------
  // define types and sparseMatrix & vector
  // --------------------------------------------------------------------------------------------
  typedef typename Assembler::field_type field_type;
  typedef Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<field_type,1,1>> SparseMatrix;
  typedef typename Dune::BlockVector<Dune::FieldVector<double,1>> Vector;
  
  SparseMatrix  matMu;
  SparseMatrix  matStiffu;
  std::vector<size_t> expandedIndices, compressedIndex;
  std::vector<size_t> expandedIndices_pre;


  // --------------------------------------------------------------------------------------------
  // time stepping loop
  // --------------------------------------------------------------------------------------------
  MaxStepsTimer.resume();
  int steps; // time step number
  bool reassemble = true;
  bool done = false;
  double end_T = 0;
  eq.time(0);
  duration<double, std::milli> ms_double_ass;
  auto recordTime1 = high_resolution_clock::now();

  for (steps=0; !done && steps<maxSteps; ++steps) // maxSteps
  {

    std::cout << "****************************************************************************************** " <<std::endl; 
    std::cout << "\t\t\t\t step: " << steps << "\t t: "  << eq.time() <<std::endl; 
    std::cout << "****************************************************************************************** " <<std::endl; 
    std::cout << " \t\t\t\t\t steps "<< steps <<std::endl;
    eq.Mass_stiff(0);
    // --------------------------------------------------------------------------------------------
    // last time step
    // --------------------------------------------------------------------------------------------
    if (eq.time()>T) 
      done = true;

    if(eq.time()+dt < T){
      end_T = eq.time()+dt;
    }else{
       end_T = T;
       dt = T - eq.time();
    }
    
    SemiImplicitEulerStep<Equation> eqSemi(&eq,dt);
    eqSemi.setTau(1);
    // --------------------------------------------------------------------------------------------
    // spectral time grid for defect correction methods with Radau points
    // --------------------------------------------------------------------------------------------
    RadauTimeGrid grid(options.nCollocUstart,eq.time(),end_T); // start sweep with 
    // std::cerr << "time points are: " << grid.points() << '\n';
    // --------------------------------------------------------------------------------------------
    // # size of all varialbles & size of variable u
    // --------------------------------------------------------------------------------------------
    int const nvars = Equation::OriginVars::noOfVariables;
    size_t     size = variableSet.degreesOfFreedom(0,nvars);
    size_t     size_e = variableSet.degreesOfFreedom(0,1);

    size_t     size_adaptivity = size;
    size_t     size_e_adaptivity= size_e;
    
    // --------------------------------------------------------------------------------------------
    // subgrid selection on new time step: everything (new game)
    // --------------------------------------------------------------------------------------------
    int sweep = -1; // sweep step number 
    expandedIndices.resize(size); std::iota(expandedIndices.begin(),expandedIndices.end(),0);
    compressedIndex.resize(size); std::iota(compressedIndex.begin(),compressedIndex.end(),0);
    expandedIndices_pre.resize(size); std::iota(expandedIndices.begin(),expandedIndices.end(),0);
    // --------------------------------------------------------------------------------------------
    // current time step
    // --------------------------------------------------------------------------------------------
    double const t = eq.time();

    // --------------------------------------------------------------------------------------------
    // mesh adaptation loop. Just once if i>0.
    // --------------------------------------------------------------------------------------------
    gridManager.setVerbosity(options.verbosity);
    // --------------------------------------------------------------------------------------------
    // Initialize values at all collocation points to the value at the starting point t0
    // --------------------------------------------------------------------------------------------
    std::vector<StateUe> collocationUe(grid.points().N(),component<0>(x)); // the (possibly non-polynomial) initial guess    
    // --------------------------------------------------------------------------------------------
    // right hand sides at the collocation points
    // --------------------------------------------------------------------------------------------   
    std::vector<Vector> ru(grid.points().N());  // right hand sides at the collocation points

    // --------------------------------------------------------------------------------------------
    // do SDC sweeps until accuracy and the number of sweeps is reached 
    // --------------------------------------------------------------------------------------------   
    bool accurate = false; 
    double normU2, sdcContraction = options.sdc_contraction;
    std::vector<double> sweepNorm;
    bool debug = false;

    // std::cerr <<"sweep\t"<<"ndof\t" <<"||du||\t\t" << "||u||\t\t" <<"sdcContraction\t\t" <<"number of cells"<<"\n";  
    std::cerr <<"sweep\t"<<"||du||\t\t" << "||u||\t\t" <<"sdcContraction\t\t"<<"\n"; 

    std::set<int> s_temp;
    for (int i = 0; i < gridManager.grid().size(0); ++i)
    {
      s_temp.insert(i);
    }

    do
    {
      // --------------------------------------------------------------------------------------------
      // increase sweep
      // --------------------------------------------------------------------------------------------  
      ++statistics.totalSweeps;
      sweep++; 
      if (sweep==0)   
      { 
        Cellfltr.set_cells(s_temp);
      }     
      // --------------------------------------------------------------------------------------------
      // set the time step
      // --------------------------------------------------------------------------------------------  
      eq.time(t);
      // --------------------------------------------------------------------------------------------
      // Evaluate mass and stiffness matrices at the initial time point (only if FE grid has changed, or deformation changed)
      // --------------------------------------------------------------------------------------------
      assemblyTimer.resume();
      if (reassemble)
      {
        // --------------------------------------------------------------------------------------------
        // perhaps we can use semi-implict equation here
        // --------------------------------------------------------------------------------------------
        auto dx(x);
        dx *= 0;
      
        // mass matrix M 
        eq.Mass_stiff(1);
        eqSemi.setTau(0);
        assembler.template assemble<AssemblyDetail::TakeAllBlocks,CellFilter>(SemiLinearization(eqSemi,x,x,dx),Cellfltr,Assembler::MATRIX|Assembler::RHS,options.assemblyThreads);   
        matMu = assembler.template get<SparseMatrix>(false);

        // stiffness matrix A 
        dx *= 0;
        eq.Mass_stiff(0);
        eqSemi.setTau(1);
        assembler.template assemble<AssemblyDetail::TakeAllBlocks,CellFilter>(SemiLinearization(eqSemi,x,x,dx),Cellfltr,Assembler::MATRIX|Assembler::RHS,options.assemblyThreads);
        matStiffu = assembler.template get<SparseMatrix>(false);

        // --------------------------------------------------------------------------------------------
        // subgrid selection on reassembly: everything (anyways, the old indices are no longer relevant)
        // --------------------------------------------------------------------------------------------
        expandedIndices.resize(matMu.N()); std::iota(expandedIndices.begin(),expandedIndices.end(),0);
        compressedIndex.resize(matMu.N()); std::iota(compressedIndex.begin(),compressedIndex.end(),0);
        expandedIndices_pre.resize(matMu.N()); std::iota(expandedIndices.begin(),expandedIndices.end(),0);
        reassemble = false;

        int const nvars = Equation::OriginVars::noOfVariables;
        size = variableSet.degreesOfFreedom(0,nvars);
        size_e = variableSet.degreesOfFreedom(0,1);
 
        size_adaptivity = size;
        size_e_adaptivity = variableSet.degreesOfFreedom(0,1);
      }
      assemblyTimer.stop();
      if(sweep ==0){        
        size_e_adaptivity = variableSet.degreesOfFreedom(0,1);
      }

      // --------------------------------------------------------------------------------------------  
      // number of restricted degrees of freedom
      // --------------------------------------------------------------------------------------------  
      size_t const nrDofs = expandedIndices.size();
      // --------------------------------------------------------------------------------------------  
      // Loop over collocation nodes (excluding the initial point) and compute the right hand sides in ru
      // --------------------------------------------------------------------------------------------  
      assemblyRhsTimer.resume();
      computeRHS(steps,x,
                 Cellfltr,
                 number_cells,
                 collocationUe,
                 size_e,
                 grid,sweep,t,
                 expandedIndices,eq,eqSemi,
                 assembler,options,ru, uAll, dt, ms_double_ass);
      assemblyRhsTimer.stop();
      // -------------------------------------------------------------------------------------------- 
      // integration integration matrix: lu or euler -> S_hat
      // -------------------------------------------------------------------------------------------- 
      sdcTimer.resume();
      SDCTimeGrid::RealMatrix Shat;
      switch (options.sweepType)
      {
        case 0: eulerIntegrationMatrix(grid, Shat); break;
        case 1: luIntegrationMatrix(grid, Shat); break;
        default: abort();
      }
     
      // -------------------------------------------------------------------------------------------- 
      // perform one SDC sweep for transmembrane voltage, gating variables, and active stress independently
      // -------------------------------------------------------------------------------------------- 

      // -------------------------------------------------------------------------------------------- 
      // * this is the solution routine used for the transmembrane voltage sweep
      // -------------------------------------------------------------------------------------------- 
      // std::string debug_system = ; 
      auto solver = [&options](SparseMatrix const& J, Vector& du, Vector& r) {
        if (options.maxCGIter>0) {
          int verbose = 0;
          if(verbose>0) std::cout <<"\n";
          DefaultDualPairing<Vector,Vector> dp;
          Dune::MatrixAdapter<SparseMatrix,Vector,Vector> Jop(J);
          JacobiPreconditioner<Dune::MatrixAdapter<SparseMatrix,Vector,Vector>> preco(Jop,1.0);
          PCGEnergyErrorTerminationCriterion<double> term(options.cgTol,options.maxCGIter);
          Pcg<Vector,Vector> pcg(Jop,preco,dp,term,0);
          Dune::InverseOperatorResult res;
          pcg.apply(du,r,res);
          // if(options.verbosity) std::cout << "cg iterations" << res.iterations << std::endl;
      

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

        } else {
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
        } 
      };
      
      // -------------------------------------------------------------------------------------------- 
      // compute restricted rhs contributions
      // -------------------------------------------------------------------------------------------- 
      std::vector<Vector> matMdiffu(grid.points().N(),Vector(nrDofs)), duVec(matMdiffu);//collocationU.size()
      // compute restricted rhs contributions. Note that we need the WHOLE residual, hence working with the restricted
      // matMuu would not work. Hence we implement the matrix-vector multiplication (with all columns but a row subset)
      // on our own: M*(u_i-u_{i+1})

      for (int j=0; j<expandedIndices.size(); ++j)
      {
        for (int i=0; i<collocationUe.size()-1; ++i){
          matMdiffu[i][j] = 0;
        }

        auto row = matMu[expandedIndices[j]];
        for (auto ci=row.begin(); ci!=row.end(); ++ci)
        {
          for (int i=0; i<collocationUe.size()-1; ++i){

            if(ci.index()<size_e){
                matMdiffu[i][j] += *ci * ((collocationUe[i].coefficients()[ci.index()])       
                                         -(collocationUe[i+1].coefficients()[ci.index()]));  
            }
          }
        }
      }
      // -------------------------------------------------------------------------------------------- 
      // compute restricted matrices
      // -------------------------------------------------------------------------------------------- 
      SparseMatrix matMuu(expandedIndices,compressedIndex,matMu);
      SparseMatrix matStiffuu(expandedIndices,compressedIndex,matStiffu);
      // std::cout <<"matMuu " << matMuu.N() << ","<< matMuu.M() <<std::endl;
      // std::cout <<"matStiffuu " << matStiffuu.N() << ","<< matStiffuu.M() <<std::endl;
      // -------------------------------------------------------------------------------------------- 
      // Reaction matrix f_u. Only the restricted subgrid dofs are considered
      // -------------------------------------------------------------------------------------------- 
      assemblyReactionTimer.resume();
      typedef Dune::BDMatrix<typename SparseMatrix::block_type> DiagonalMatrix;
      std::vector<DiagonalMatrix> allMatFu(grid.points().N(),DiagonalMatrix(nrDofs));
      assemblyReactionTimer.stop();
      
      // -------------------------------------------------------------------------------------------- 
      // perform SDC sweep
      // -------------------------------------------------------------------------------------------- 
      // writeVTKFile(x,output + "/x-step-"+paddedString(steps));
      // printuAll(x,uAll,options, output + "/x-step-"+paddedString(steps),"u_pre");
      bool shift = true;
      if(steps==maxSteps-1)
        shift = true;
      std::string name = "test_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep);
      sweepNorm.push_back( sdcIterationStepJacobi(grid,Shat,solver,matMuu,matStiffuu,ru,allMatFu,matMdiffu,duVec));
      // -------------------------------------------------------------------------------------------- 
      // update solution: add Newton correction
      // --------------------------------------------------------------------------------------------
      for (int ii=1; ii<collocationUe.size(); ii++) // start loop at 1 as initial values is not changed
      {
        State du(x);
        State du_mark(x);
        du=0;
        du_mark=0;
        for (int j=0; j<expandedIndices.size(); ++j)
        {
          size_t ej = expandedIndices[j];
          collocationUe[ii].coefficients()[ej] += duVec[ii][j];
          at_c<0>(du.data).coefficients()[ej] = duVec[ii][j];
          if(duVec[ii][j]!=0) at_c<0>(du_mark.data).coefficients()[ej] = 1.0;    
        }
        
        // if(steps==std::floor(maxSteps/2)) printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
        // if(steps==std::floor(maxSteps/2)) printuAll(du_mark,uAll,options.order, output + "/du_mark_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du_mark");
        if(options.plot)
          printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
      }
      // --------------------------------------------------------------------------------------------
      // sdcContraction
      // --------------------------------------------------------------------------------------------      
      size_t nDofs = variableSet.degreesOfFreedom(0,nvars);

      State sol_tmp(x); 
      component<0>(sol_tmp) = collocationUe.back();
      if(options.plot) printuAll(sol_tmp,uAll,options.order, output + "/newton_update_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep),"u_update");

      if(steps==std::floor(maxSteps/2) and options.plot) printuAll(sol_tmp,uAll,options.order, output + "/newton_update_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep),"u_update");

      BlockVectorX tmp(nDofs), sol_tmp_copy(nDofs);

      tmp *= 0; 
      sol_tmp_copy *= 0;
      sol_tmp.write(sol_tmp_copy.begin());
      matMu.mv(sol_tmp_copy,tmp);
      normU2 = sol_tmp_copy * tmp;

      // std::cerr << "selectedDOF= "<< size_e_adaptivity <<
      //              "  Muu = " << matMuu.N() << 
      //              "  Kuu = " << matStiffuu.N() << 
      //              "  Mdiffu = " << matMdiffu.size() << 
      //             "  Mdiffu[0] = " << matMdiffu[0].size() << 
      //              "  allMatFu = " << allMatFu.size() << 
      //              "  allMatFu[0] = " << allMatFu[0].N() << 
      //              "  Shat = " << Shat.N() << 
      //              "  Shat[0] = " << Shat[0].N() << 
      //              "  ru = " << ru[0].size() << 
      //              "  ||du|| = " << sweepNorm.back() << 
      //              "  ||u|| = " << std::sqrt(normU2) << "\n";  


      if (sweepNorm.size()>1)
      {
        double c = sweepNorm.back()/sweepNorm[sweepNorm.size()-2];
        sdcContraction = std::sqrt(c*sdcContraction);
        // std::cerr << "=================================================================================\n";
        // std::cerr << "SDC contraction: " << c << " estim: " << sdcContraction << "  long range: " << std::pow(sweepNorm.back()/sweepNorm[0],1.0/(sweepNorm.size()-1)) << "\n";
        // std::cout << "||du|| " << sweepNorm.back()   
        //           << "\nsdcContraction/(1-sdcContraction):   " << sdcContraction/(1-sdcContraction)
        //           << "\nsweepNorm.back()*sdcContraction/(1-sdcContraction):  " << sweepNorm.back()*sdcContraction/(1-sdcContraction) << std::endl;
      }

      //std::cerr << sweep <<"\t"<< expandedIndices.size()  <<"\t"<< sweepNorm.back() << "\t"<< std::sqrt(normU2) << "\t"<<sdcContraction << "\t\t" << Cellfltr.get_size()<<"\n";  
      std::cerr << sweep <<"\t"<< sweepNorm.back() << "\t" <<std::sqrt(normU2) << "\t"<<sdcContraction <<"\n"; 

      // // --------------------------------------------------------------------------------------------
      // // plot du 
      // // --------------------------------------------------------------------------------------------  
      // if (true || sdcContraction<=1 && sweepNorm.back()*sdcContraction/(1-sdcContraction)<=options.aTol) 
      // {
      //   StateU du = collocationU[0];
      //   for (int ii=1; ii<collocationU.size(); ii++) // start loop at 1 as initial values is not changed
      //   {
      //     du = 0;
      //     for (int j=0; j<expandedIndices.size(); ++j) {
      //       du.coefficients()[expandedIndices[j]] += duVec[ii][j];
      //     }
      //     //writeVTK(du,output+"/aliev-du-sweep="+paddedString(sweep)+"-col="+paddedString(ii)+"-step="+paddedString(steps),IoOptions().setOrder(2),"du");
      //   }
      // }

      // --------------------------------------------------------------------------------------------  
      // select degrees of freedom to take into account in the next sweep. This is a 
      // further reduction of the up to now used dofs
      // --------------------------------------------------------------------------------------------  
      int count = 0;
      int discount = 0;
      if (options.tolSelect > 0)
      {
        cellsMarked.clear();
        std::map<int, int> map_expanded_indices;

        State markSelectedDOF(x);
        markSelectedDOF*=0;

        std::vector<size_t> newExpandedIndices;
        std::vector<size_t> newCompressedIndex;
        std::set<size_t> set_ExpandedIndices;
      
        expandedIndices_pre.assign(expandedIndices.begin(),expandedIndices.end());

        newCompressedIndex.resize(compressedIndex.size(),compressedIndex.size());
        compressedIndex.assign(newCompressedIndex.begin(),newCompressedIndex.end());

        for (int ii=0; ii<duVec.back().size(); ++ii)
        {
          double duMax = 0;
          for (auto const& duj: duVec)
            duMax = std::max(duMax,std::abs(duj[ii]));

          //if(sdcContraction>=1 || sdcContraction*duMax/(1-sdcContraction) > options.tolSelect)
          if(sdcContraction>1 || sdcContraction*duMax/(1-sdcContraction) > options.tolSelect)
          {
            size_t ej = expandedIndices[ii];
            std::set<int> s = index2IndexsSet[ej];
            std::set<int>::iterator it;
            for (it = s.begin(); it != s.end(); ++it) {
              set_ExpandedIndices.insert(*it);
            }
          }
        }
        
        std::set<size_t>::iterator it;
        for (it=set_ExpandedIndices.begin(); it!=set_ExpandedIndices.end(); ++it){
          newExpandedIndices.push_back(*it);
        }


        size_e_adaptivity = 0;

        for (int i = 0; i < newExpandedIndices.size(); ++i)
        {          
          size_t ej_next = newExpandedIndices[i];
        
          std::set<int> cell_set = index2Cells_new[ej_next];
          std::set<int> selected_cell_idx;
          selected_cell_idx.insert(cell_set.begin(), cell_set.end());

          std::set<int>::iterator it_cell;
          for (it_cell = cell_set.begin(); it_cell != cell_set.end(); ++it_cell) {
            cellsMarked.push_back(*it_cell);
          }
          
          at_c<0>(markSelectedDOF.data).coefficients()[ej_next] = 1.0;
          size_e_adaptivity++;
        }

        std::set<int> s_temp_al(cellsMarked.begin(), cellsMarked.end()); 

        Cellfltr.set_cells(s_temp_al);
        // Cellfltr.get_cells();
        expandedIndices.assign(newExpandedIndices.begin(),newExpandedIndices.end());
     
        for (int i = 0; i < expandedIndices.size(); ++i)
        {
          compressedIndex[expandedIndices[i]] = i;
        }

        count = expandedIndices.size();
        discount = compressedIndex.size()-expandedIndices.size();
      
        if(options.plot and compressedIndex.size()!=expandedIndices.size()) printuAll(markSelectedDOF,uAll,options.order, output + "/Selected-steps_"+paddedString(steps)+"-sweep_"+paddedString(sweep),"SelectedDOF");

        expandedIndices_pre.assign(expandedIndices.begin(),expandedIndices.end());
      }
      eq.time(t+dt);

      // // --------------------------------------------------------------------------------------------  
      // // estimate spatial error. We impose an absolute tolerance that is on the order of the SDC iteration error.
      // // --------------------------------------------------------------------------------------------  
      // if (options.adapt && !(options.rosenbrockRefinementStyle&&accurate)) 
      // {
      //   if (options.rosenbrockRefinementStyle)
      //     tolX[0] = std::make_pair(options.aTol,0);
      //   else
      //     tolX[0] = std::make_pair(std::max(options.aTol,sweepNorm.back()),0);
        
      //   State spatialError(x), xnext(x); 
      //   at_c<0>(spatialError.data) = collocationUe.back();

      //   projectHierarchically(variableSet,spatialError);
      //   at_c<0>(spatialError.data) -= collocationUe.back();
      //   at_c<0>(xnext.data) = collocationUe.back();
        
      //   // perform mesh adaptation
      //   accurate = embeddedErrorEstimator(variableSet,spatialError,xnext,IdentityScaling(),tolX,gridManager,0);
      //   if (!accurate) {
      //     size = variableSet.degreesOfFreedom(0,nvars);
      //     //std::cout << "\t\t\t  accurate = " << accurate << ",   dofs after mesh refinement: " << size << std::endl;
      //     sweepNorm.clear(); // old SDC sweep norm on coarser grid cannot be compared to sweeps on finer grid.
      //     reassemble = true;
          
      //     if (options.rosenbrockRefinementStyle)
      //       for (int i=1; i<collocationUe.size();   ++i)
      //       {
      //         collocationUe[i] = collocationUe[0];
      //       }
      //   }
      //   else if (options.rosenbrockRefinementStyle && options.writeVTK>0)
      //   {
      //     // if(options.plot) writeVTKFile(xnext,output+"/rosenbrock-euler-"+paddedString(steps),IoOptions(),2);
      //   }
      // } else
       accurate = true;

      // --------------------------------------------------------------------------------------------  
      // Refine time grid for next sweep if nominal value not reached. Perform interpolation of coefficients.
      // --------------------------------------------------------------------------------------------  
      if (grid.points().N() < options.nCollocU+1)
      {
        SDCTimeGrid::RealMatrix p;
        grid.refine(p);
        std::vector<StateUe> newUe(collocationUe.size()+1,collocationUe[0]);
        for (int i=0; i<newUe.size(); ++i)
        {
          newUe[i] = 0;
          for (int j=0; j<collocationUe.size(); ++j)
          {
            newUe[i].axpy(p[i][j],collocationUe[j]);
          }
        }
        collocationUe.swap(newUe);
        ru.push_back(ru.front()); // extend those as well
        // std::cerr << "time points in ladder method " << grid.points() << '\n';
      }

      if(sweepNorm.back() < options.SDC_TOL){
        std::cout << "options.SDC_TOL"<<std::endl;
        // std::cout << "break sweep loop: " << sweepNorm.back()   << " sweepNorm.back()*sdcContraction/(1-sdcContraction):  " << sweepNorm.back()*sdcContraction/(1-sdcContraction) << " sdcContraction/(1-sdcContraction):   " << sdcContraction/(1-sdcContraction) << std::endl;
        break;
      }
      //std::cout << "matMdiffu next!!! collocationUe0 grid.points().N() " << grid.points().N()  <<" nrDofs "<< nrDofs << " expandedIndices.size() "<< expandedIndices.size() << std::endl;
    

      if(count==0 and discount!=0){
        std::cout << "count==0"<<std::endl;
        break;
      }   
    }
    //while ( sweep+1<options.maxSweeps && (sweep+1<options.minSweeps ||  sdcContraction>1 || !accurate) );
    while ( sweep+1<options.maxSweeps && (sweep+1<options.minSweeps ||  sdcContraction>1 || sweepNorm.back()*sdcContraction/(1-sdcContraction)>options.aTol || !accurate) );
    

    // --------------------------------------------------------------------------------------------
    // extract final-time value
    // --------------------------------------------------------------------------------------------   
    component<0>(x) = collocationUe.back();

    // // --------------------------------------------------------------------------------------------
    // // perform mesh coarsening
    // // --------------------------------------------------------------------------------------------
    // if (options.adapt) {
    //   std::vector<bool> dummy(0);
    //   coarsening(variableSet,x,IdentityScaling(),tolXC,gridManager,dummy,options.verbosity,options.minRefLevel);
    //   reassemble = true;
    // }
    
    // --------------------------------------------------------------------------------------------
    // t & dt 
    // --------------------------------------------------------------------------------------------
    std::cout.flush();
    // std::cerr << "t= " << eq.time() << " dt = " << dt << '\n';
    
    // dt = std::min(dt,dtMax);
    if(options.plot) printuAll(x,uAll,options.order, output + "/emiSDC"+paddedString(steps),"u");
  }
  MaxStepsTimer.stop();
  auto recordTime2 = high_resolution_clock::now();
  MaxStepsTimer.stop();
  std::cout << " ****************************************************************************************** " <<std::endl;
  duration<double, std::milli> ms_double = recordTime2 - recordTime1;
  std::cout <<"ms_double:\t" << ms_double.count() << " ms\n";
  ms_double = ms_double - ms_double_ass;
  
  std::cout <<"ms_double:after\t" << ms_double.count() << " ms\n";

  if (!done)
    std::cerr << "\n*** maxSteps reached ***\n\n";
  
  std::cout << "assembly matrix time:   " << assemblyTimer.format() 
            << "assembly rhs time:      " << assemblyRhsTimer.format() 
            << "assembly reaction time: " << assemblyReactionTimer.format() 
            << "sdc time:               " << sdcTimer.format() 
            << "sdc ODE time:           " << sdcOdeTimer.format()
            << "refinement time:        " << refineTimer.format() 
            << "coarsening time:        " << coarsenTimer.format() 
            << "output time:            " << outputTimer.format()
            << "forloop t=1:maxsteps time:            " << MaxStepsTimer.format() << '\n';
  
  statistics.avgDofs = dofCount / (double)steps;
  statistics.sdcTime = sdcTimer.elapsed().wall;
  std::cout<< "statistics.sdcTime:  "<< statistics.sdcTime << std::endl;
  printuAll(x,uAll,options.order, output+"/emiSDCJACOBILast","u");
  writeVectorTofile(x,"matlab_dir/emiSDCJACOBILast");
  // std::cout <<" full domain: " << Total_dofCount <<std::endl;

  // std::string s_target;
  // if(options.tolSelect==0){
  //   s_target = "E_minSW_sol_stol_"+std::to_string(options.tolSelect)+"_dt_"+std::to_string(options.dt)+"_Sweep_"+std::to_string(options.maxSweeps)+"_coll_"+std::to_string(options.nCollocU)+"_TOL_"+std::to_string(options.SDC_TOL)+"_Stps_"+ std::to_string(maxSteps);
  //   writeDoubleTofile(Total_dofCount,s_target+"FULLDOFs");
  // }
  // else{
  //   s_target = "E_minSW_sol_stol_"+std::to_string(std::abs(log10(options.tolSelect)))+"_dt_"+std::to_string(options.dt)+"_Sweep_"+std::to_string(options.maxSweeps)+"_coll_"+std::to_string(options.nCollocU)+"_TOL_"+std::to_string(options.SDC_TOL)+"_Stps_"+ std::to_string(maxSteps);
  //   writeDoubleTofile(Total_dofCount,s_target+"DROPDOFs");
  // }
    
 return x;
}

#endif
