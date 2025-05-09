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

#ifndef INTEGRATE_CG_BDDC_SDC_SPLIT_HH
#define INTEGRATE_CG_BDDC_SDC_SPLIT_HH

using namespace Kaskade;

template <class Interfaces, class Matrix, class Vectors, class Vector, class ReactionDerivatives>
typename Matrix::field_type sdcIterationStepBDDCSPLIT(bool BDDC_SDC_with_initial, bool BDDC_verbose, bool BDDC_SDC_verbose, int sweep, bool cg_solver, int iter_cg_with_bddc, double tol, std::string matlab_dir, Interfaces const& interfaces,int n_subdomains ,SDCTimeGrid const& grid, SDCTimeGrid::RealMatrix const& Shat,
                                                std::map<int,std::set<int>> map_II, 
                                                std::map<int,std::set<int>> map_GammaGamma_noDuplicate, 
                                                std::vector<std::vector<LocalDof>> sharedDofsKaskade, 
                                                std::vector<Vector> weights,
                                                std::map<int, int> map_indices,
                                                Matrix const& A,
                                                double dt,
                                                std::vector<size_t> const& expandedIndices,
                                                std::map<int, int> map_index_to_subdomain,
                                                std::vector<int> sequenceOfTags,
                                                std::map<int,std::unordered_map<int, int>> local2Global,
                                                std::map<int,std::unordered_map<int, int>> global2Local,
                                                std::map<int,std::vector<int>> IG_seq,
                                                std::vector<Matrix> const& M_bddc,                //matMuu_bddc
                                                std::vector<Matrix> const& Stiff_bddc,            //matStiffuu_bddc
                                                std::vector<Vectors> & rUi_bddc,                  //ru_bddc
                                                std::vector<ReactionDerivatives> const& rDu_bddc, //allMatFu_bddc
                                                std::vector<Vectors> const& Mdu_bddc,             //matMdiffu_bddc
                                                std::vector<Vectors> & du_bddc,
                                                std::vector<Vectors> & du_bddc_initial)                   //duVec_bddc

{
  auto const& pts = grid.points();
  int const n = pts.size()-1;      // number of subintervals

  assert(Mdu_bddc[0].size()>=n); 
  assert(du_bddc[0].size()>=n+1);  // including start point

  // compute exact integration matrix
  auto const& S = grid.integrationMatrix();

  // typedef typename Vectors::value_type Vector;
  // declare here to prevent frequent reallocation
  std::vector<Vector> rhs_bddc(n_subdomains);
  std::vector<Vector> tmp_bddc(n_subdomains);
  for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
  {
    size_t const dofs = Mdu_bddc[subIndx][0].size(); 
    Vector rhs(dofs), tmp(dofs); 
    rhs_bddc[subIndx] = rhs;
    tmp_bddc[subIndx] = tmp;
  }

  std::vector<int> activeIds(n_subdomains);
  std::iota(activeIds.begin(), activeIds.end(), 0);
  
  // std::vector<int> activeIds;
  // activeIds.resize(n_subdomains);
  // parallelFor(0,n_subdomains,[&](int subIdx)
  // {
  //   activeIds[subIdx] = subIdx;
  // });


  Vector initial(A.N()), initial_temp(A.N());
  initial *= 0;
  initial_temp *=0;  

  // perform n Euler steps
  typename Matrix::field_type norm = 0;

  std::vector<Matrix> JJ(n_subdomains);
  for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
  {
    JJ[subIndx] = M_bddc[subIndx];
    // initialize correction at starting point to zero
    du_bddc[subIndx][0] = 0.0; // for each subdomain
    du_bddc_initial[subIndx][0] = 0.0;
  }

  if(BDDC_SDC_verbose){
  std::cout << "\n\n";
  std::cout << "\t\t\t\t\t\t\t\t\t\t========================================================\n";
  std::cout << "\t\t\t\t\t\t\t\t\t\titr \t" << "step length\t" << "||du||_A\t"<< "p[0]\n";
  std::cout << "\t\t\t\t\t\t\t\t\t\t========================================================\n";
  }

  // using TransmissionScalar = double;
  // using BddcSubdomain = Subdomain<1,double,double,SpaceTransfer<1,double,TransmissionScalar>>;

  using TransmissionScalar = uint32_t;//int16_t;
  using BddcSubdomain = Subdomain<1,double,double,SpaceTransferDataCompression<1,double,TransmissionScalar>>;
  // std::cout <<"sizeof(TransmissionScalar): " << sizeof(TransmissionScalar) << " sizeof(double): "<< sizeof(double) << std::endl;

  // for each collocation point
  for (int i=1; i<=n; i++)
  { 
    if(BDDC_SDC_verbose){
    std::cout << "\t\t\t\t\t\t\t\t\t\tcol \t" << i <<"\n";
    std::cout << "\t\t\t\t\t\t\t\t\t\t========================================================" <<std::endl;
    }
    std::vector<std::unique_ptr<BddcSubdomain>> subsptr(n_subdomains);

    // -----------------------------------------------------------------------
    // set initial guess from previpus collocation point
    // -----------------------------------------------------------------------
    if(BDDC_SDC_with_initial)
    {
      for (int gIdx=0; gIdx<A.N(); ++gIdx)
      {
        size_t ej = expandedIndices[gIdx];
        int subIndx = map_index_to_subdomain[ej];
        int tag =  sequenceOfTags[subIndx];
        initial_temp[ej] = du_bddc[subIndx][i-1][global2Local[tag][ej]]; 
      }
      initial *= 0;
      //  initial = initial+A * initial_temp.
      A.umv(initial_temp,initial);
      
      Vector rhs_petsc_test(A.N());
      petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, initial, rhs_petsc_test);

      for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
      {
        int tag = sequenceOfTags[subIdx];
        std::vector<int> IG = IG_seq[subIdx];
        int counter_kasakde = IG.size();
        Vector Fs_subIdx(counter_kasakde);  
        { 
          for (int c = 0; c < counter_kasakde; ++c)
          {
            int index = IG[c];
            float coef = weights[subIdx][index];
            Fs_subIdx[c] = coef*rhs_petsc_test[index]*Shat[i-1][i];
          }
        }
        du_bddc_initial[subIdx][i] = Fs_subIdx; // for each collocation
      } 
    }

    // -----------------------------------------------------------------------
    for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
    {
      // ---------------------------------------------------------------------
      // LHS: left-hand side for linear system
      // ---------------------------------------------------------------------
      // matrix J = M - Shat_i-1,i*(A+f_u)
      // ---------------------------------------------------------------------
      for (size_t row=0; row<JJ[subIndx].N(); ++row)
      {
        auto colJ = JJ[subIndx][row].begin(); 
        auto end = JJ[subIndx][row].end();
        auto colM = M_bddc[subIndx][row].begin();
        auto colA = Stiff_bddc[subIndx][row].begin();
        auto colR = rDu_bddc[subIndx][i][row].begin();
        auto endR = rDu_bddc[subIndx][i][row].end();
         
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

      // ---------------------------------------------------------------------
      // RHS:  right-hand side for linear system
      // ---------------------------------------------------------------------
      // M * ( u_i^{k} - u_{i+1}^k + du_i) = Mdiffu + Mdu
      // ---------------------------------------------------------------------
      rhs_bddc[subIndx] = Mdu_bddc[subIndx][i-1];
      M_bddc[subIndx].umv(du_bddc[subIndx][i-1],rhs_bddc[subIndx]);


      // // add sum_j S_ij r_j to right hand side
      for (int j=0; j<=n; ++j){
        rhs_bddc[subIndx].axpy(S[i-1][j],rUi_bddc[subIndx][j]);
      }
         
      // add sum_j Shat_ij r'(u_j) du_j with r' = A + f_u
      tmp_bddc[subIndx] = 0;
      for (int j=0; j<i; ++j) // TODO: start at 1 instead of 0? du[0] is zero anyway...
      {
        rDu_bddc[subIndx][j].usmv(Shat[i-1][j],du_bddc[subIndx][j],rhs_bddc[subIndx]);
        tmp_bddc[subIndx].axpy(Shat[i-1][j],du_bddc[subIndx][j]);
      }
      Stiff_bddc[subIndx].umv(tmp_bddc[subIndx],rhs_bddc[subIndx]);

  
    }// end of subdomain   
    // ------------------------------------------------------------------------------------------------------------
    



    // ------------------------------------------------------------------------------------------------------------
    // for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
    // {
    parallelFor(0,n_subdomains,[&](int subIndx)
    {
      if(BDDC_SDC_with_initial) rhs_bddc[subIndx]-=du_bddc_initial[subIndx][i]; // previous increment is probably a good starting value
      subsptr[subIndx] = std::make_unique<BddcSubdomain>(subIndx,JJ[subIndx],interfaces);
    });

    std::vector<BddcSubdomain> subs;
    for (auto& sp: subsptr)
      subs.push_back(*sp);

    BDDCSolver<BddcSubdomain> bddcSolver(subs,interfaces.coarseConstraints(),activeIds,cg_solver, BDDC_SDC_verbose);
    bddcSolver.setRhs(rhs_bddc);

    std::vector<double> resNorm;
    for (int k=0; k<iter_cg_with_bddc; ++k)
    {
      resNorm.push_back(bddcSolver.solve());
      if(resNorm.back()<tol)
        break;
    }

    // update the solution
    for (int subIndx=0; subIndx<n_subdomains; ++subIndx)
    {
        auto ui = subs[subIndx].getSolution();
        du_bddc[subIndx][i] = ui; // fix me!
        if(BDDC_SDC_with_initial) du_bddc[subIndx][i]+=du_bddc_initial[subIndx][i];
    } 
    double norm_sub = 0;
    for (int subIndx=0; subIndx<n_subdomains; ++subIndx)
    {
      norm_sub += du_bddc[subIndx][i]*rhs_bddc[subIndx];
    }

    norm += (pts[i]-pts[i-1])*norm_sub;

    if(BDDC_SDC_verbose){
    std::cout << "\t\t\t\t\t\t\t\t\t\t========================================================" <<std::endl;
    }
  }   // end i - loop
  
  return std::sqrt(norm/(pts[n]-pts[0]));
} 

template <class State, class StateUe, class TimeGrid, class Vector, class Eq,class EqSemi, class Assem, class elementType, class CellFilter, class Options>
void computeRHS_BDDC_split(int step, 
                State const& x,
                CellFilter & Cellfltr,
                int number_cells,
                std::vector<StateUe> const& collocationU,
                size_t size_e, 
                TimeGrid const& grid, 
                int sweep, double t, 
                std::vector<size_t> const& expandedIndices,
                Eq& eq, 
                EqSemi& eqSemi, 
                Assem& assembler, 
                Options const& options, 
                std::vector<std::vector<Vector>> & ru_bddc,
                std::vector<std::vector<size_t>> const& expandedIndices_bddc,
                elementType & uAll, double dt,
                std::map<int,std::set<int>> map_II, 
                std::map<int,std::set<int>> map_GammaGamma_noDuplicate, 
                int n_subdomains, 
                std::vector<std::vector<LocalDof>> sharedDofsKaskade, 
                std::vector<Vector> weights,
                std::vector<int> sequenceOfTags,
                std::map<int,std::vector<int>> IG_seq,
                int nDofs,
                std::chrono::duration<double, std::milli> & ms_double_ass,
                std::map<int,std::unordered_map<int, int>> local2Global,
                std::map<int,std::unordered_map<int, int>> global2Local,
                std::map<int, int> map_indices)
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  using namespace boost::fusion;
  typedef typename Eq::OriginVars::template CoefficientVector<0,1> CoefficientVectorsU;  
  typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Eq> >  SemiLinearization;
  size_t const nrSelectedDofs = expandedIndices.size();
  
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
      // for each subdomain we need to copy the next ru ....
      for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
      {
        ru_bddc[subIndx][ii] = ru_bddc[subIndx][ii-1];
      }

    }
    else 
    {
      eq.time(t+pts[ii]-pts[0]);
      at_c<0>(stateTmp.data) = collocationU[ii];

      auto recordTime_ass1 = high_resolution_clock::now();
      assembler.template assemble<AssemblyDetail::TakeAllBlocks,CellFilter>(SemiLinearization(eqSemi,stateTmp,stateTmp,dstateTmp),Cellfltr,Assembler::RHS,options.assemblyThreads);
      auto recordTime_ass2 = high_resolution_clock::now();

      duration<double, std::milli> ms_double_ass_temp = recordTime_ass2 - recordTime_ass1;
      ms_double_ass +=ms_double_ass_temp;

      // ------------------------------------------------
      // BDDC ---> START
      // ------------------------------------------------      
      Vector rhs_vec_test(nDofs);
      auto rhs = assembler.rhs();
      rhs.write(rhs_vec_test.begin());

      Vector rhs_petsc_test(nDofs);
      rhs.write(rhs_petsc_test.begin());

      petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, rhs_vec_test, rhs_petsc_test);

      // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      // fill the sub matrices for kaskade format
      // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      std::vector<Vector> Fs(n_subdomains);
      for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
      {
        int tag = sequenceOfTags[subIdx];
        std::vector<int> IG = IG_seq[subIdx];
        int counter_kasakde = IG.size();
        Vector Fs_subIdx(counter_kasakde);  
        { 
          // std::cout << subIdx << " : ";
          for (int i = 0; i < counter_kasakde; ++i)
          {
            int index = IG[i];
            float coef = weights[subIdx][index];
            Fs_subIdx[i] = coef*rhs_petsc_test[index];
            // std::cout << Fs_subIdx[i] << ", ";
          }
          // std::cout << "\n";
        }
        Fs[subIdx] = Fs_subIdx;
      }

      
      for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
      {
        std::vector<size_t> expandedIndices_sub = expandedIndices_bddc[subIndx];
        // CoefficientVectorsU rs(Fs[subIndx]); 
        size_t nDofs_sub = expandedIndices_sub.size();

        ru_bddc[subIndx][ii] = Vector(nDofs_sub);

        for (size_t j=0; j<expandedIndices_sub.size(); ++j)
        {
          size_t ej = expandedIndices_sub[j];// for each subdoamin the local indices are saved in expandedIndices_sub
          ru_bddc[subIndx][ii][j] = Fs[subIndx][ej]*(1/dt); 
        }
      }
    }
  }
}

template <class Grid, class Equation, class VariableSet, class Spaces, class elementType, class CellFilter, class Options, class OptionStatistics, class Matrix, class Vector>
typename VariableSet::VariableSet semiImplicit_CG_BDDC_SDC_SPLIT( GridManager<Grid>& gridManager,
                                                            Equation& eq,
                                                            CellFilter & Cellfltr,  
                                                            VariableSet const& variableSet, 
                                                            Spaces const& spaces,
                                                            Grid const& grid,
                                                            typename VariableSet::VariableSet x,
                                                            std::vector<std::set<int>> index2Cells_new,
                                                            Options const& options, 
                                                            OptionStatistics& statistics,
                                                            std::string output,
                                                            elementType & uAll,
                                                            std::vector<std::set<int>> index2IndexsSet,
                                                            bool cg_semi, 
                                                            bool direct,
                                                            std::string matlab_dir,
                                                            Vector & sol_BDDC_SDC,
                                                            std::vector<std::vector<LocalDof>> sharedDofsKaskade,
                                                            int interfaceTypes,
                                                            int n_subdomains,
                                                            Matrix A_,
                                                            Matrix M_,
                                                            Matrix K_,
                                                            std::vector<Matrix> As,
                                                            std::vector<Matrix> Ms,
                                                            std::vector<Matrix> Ks,
                                                            std::map<int,std::vector<int>> IG_seq,
                                                            std::map<int,std::set<int>> IGamma,
                                                            std::vector<int> sequenceOfTags, 
                                                            std::map<int,std::set<int>> map_II, 
                                                            std::map<int,std::set<int>> map_GammaGamma_noDuplicate, 
                                                            std::vector<Vector> weights,
                                                            std::vector<Vector> Fs,
                                                            bool cg_solver,
                                                            int iter_cg_with_bddc,
                                                            std::map<int,std::unordered_map<int, int>> local2Global,
                                                            std::map<int,std::unordered_map<int, int>> global2Local,
                                                            std::map<int, int> map_index_to_subdomain,
                                                            double tol,
                                                            // std::map<int, int> sub_length_var,
                                                            std::map<int, int> map_t2l,
                                                            std::map<int, int> map_indices,
                                                            bool BDDC_SDC_with_initial,
                                                            bool BDDC_verbose,
                                                            bool BDDC_SDC_verbose)
{
  // for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
  // {
  //   std::unordered_map<int, int> global2Local_sub = global2Local[subIndx];
  //   std::unordered_map<int, int>::iterator it;
  //   for (it = global2Local_sub.begin(); it != global2Local_sub.end(); it++)
  //   { 
  //     std::cout <<"("<<it->first << "," << it->second << ")\n";
  //   }
  //   std::cout << "=====================================================\n";
  // }

  // for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  // {
  //   std::cout << "BDDC_SDC: As[subIdx].N()-> "<< As[subIdx].N() << " Ms[subIdx].N()-> " << Ms[subIdx].N()  << " Ks[subIdx].N()-> "  << Ks[subIdx].N()  << std::endl;
  // }

  using namespace boost::fusion;
  std::setprecision(16);

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

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
  
  typedef typename Equation::OriginVars::template CoefficientVector<0,1> CoefficientVectorsU;
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
  // typedef typename Dune::BlockVector<Dune::FieldVector<double,1>> Vector;
  
  SparseMatrix  matMu;
  SparseMatrix  matStiffu;
  std::vector<size_t> expandedIndices, compressedIndex;
  std::vector<size_t> expandedIndices_pre;

  std::vector<std::vector<size_t>> expandedIndices_bddc, compressedIndex_bddc, expandedIndices_pre_bddc;
  expandedIndices_bddc.resize(n_subdomains);
  compressedIndex_bddc.resize(n_subdomains);
  expandedIndices_pre_bddc.resize(n_subdomains);


  std::vector<int> s_vec;
  s_vec.resize(gridManager.grid().size(0));std::iota(s_vec.begin(),s_vec.end(),0);
  std::set<int> elementIdx(s_vec.begin(),s_vec.end());
  std::cout << "elementIdx.size(): "<< elementIdx.size() <<std::endl;
  // --------------------------------------------------------------------------------------------
  // BDDC
  // --------------------------------------------------------------------------------------------
  std::vector<int> subdomSize(n_subdomains);
  for (int subIndx=0; subIndx<n_subdomains; ++subIndx){
    int tag =  sequenceOfTags[subIndx];
    subdomSize[subIndx] = Ms[subIndx].N();
    // std::cout << "Ms[subIndx].N() "<< Ms[subIndx].N()  << std::endl;
  }

  InterfaceAverages<1,int> ifa(sharedDofsKaskade,subdomSize,interfaceTypes);
  // --------------------------------------------------------------------------------------------
  // time stepping loop
  // --------------------------------------------------------------------------------------------
  MaxStepsTimer.resume();
  int steps; // time step number
  bool reassemble = true;
  bool done = false;
  double end_T = 0;
  eq.time(0);
  
  size_t nDofs = variableSet.degreesOfFreedom(0,1);
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
    std::cerr << "time points are: " << grid.points() << '\n';
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
    // -> to keep the selected global indices
    // --------------------------------------------------------------------------------------------
    int sweep = -1; // sweep step number 
    expandedIndices.resize(size); std::iota(expandedIndices.begin(),expandedIndices.end(),0);
    compressedIndex.resize(size); std::iota(compressedIndex.begin(),compressedIndex.end(),0);
    expandedIndices_pre.resize(size); std::iota(expandedIndices.begin(),expandedIndices.end(),0);

    // --------------------------------------------------------------------------------------------
    // subgrid selection on new time step: everything (new game) BDDC, 
    // -> to keepp the selected local indices for each subdomain
    // --------------------------------------------------------------------------------------------
    for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
    {
      int tag =  sequenceOfTags[subIndx];
      std::vector<size_t> expandedIndices_sub;
      expandedIndices_sub.resize(Ms[subIndx].N());//local2Global[tag].size());
      std::iota(expandedIndices_sub.begin(),expandedIndices_sub.end(),0);
      expandedIndices_bddc[subIndx] = expandedIndices_sub;


      std::vector<size_t> compressedIndex_sub;
      compressedIndex_sub.resize(Ms[subIndx].N());//local2Global[tag].size());
      std::iota(compressedIndex_sub.begin(),compressedIndex_sub.end(),0);
      compressedIndex_bddc[subIndx] = compressedIndex_sub;

      std::vector<size_t> expandedIndices_pre_sub;
      expandedIndices_pre_sub.resize(Ms[subIndx].N());//local2Global[tag].size());
      std::iota(expandedIndices_pre_sub.begin(),expandedIndices_pre_sub.end(),0);
      expandedIndices_pre_bddc[subIndx] = expandedIndices_pre_sub;

      // std::cout << "subIndx: "<< subIndx << " local2Global[tag].size(): "<<local2Global[tag].size()<< std::endl;
    }

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
    std::vector<StateUe> collocationU(grid.points().N(),component<0>(x)); // the (possibly non-polynomial) initial guess    
    // --------------------------------------------------------------------------------------------
    // right hand sides at the collocation points
    // --------------------------------------------------------------------------------------------   
    std::vector<std::vector<Vector>> ru_bddc(n_subdomains);  // right hand sides at the collocation points for BDDC

    for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
    {
      std::vector<Vector> ru_sub(grid.points().N());
      ru_bddc[subIndx] = ru_sub;
    }

    // --------------------------------------------------------------------------------------------
    // do SDC sweeps until accuracy and the number of sweeps is reached 
    // --------------------------------------------------------------------------------------------   
    bool accurate = false; 
    double normU2, sdcContraction = options.sdc_contraction;
    std::vector<double> sweepNorm;
    std::vector<double> sweepNorm_bddc;
    bool debug = false;

    std::cerr <<"sweep\t"<<"||du||\t\t" << "||u||\t\t" <<"sdcContraction\t\t"<<"\n";   

    do
    {
      // --------------------------------------------------------------------------------------------
      // increase sweep
      // --------------------------------------------------------------------------------------------  
      ++statistics.totalSweeps;
      sweep++;
      if (sweep==0)   
      { 
        Cellfltr.set_cells(elementIdx);
      }
      
      // --------------------------------------------------------------------------------------------
      // set the time step
      // --------------------------------------------------------------------------------------------  
      eq.time(t);

      if(sweep ==0)      
        size_e_adaptivity = variableSet.degreesOfFreedom(0,1);

      // --------------------------------------------------------------------------------------------  
      // number of restricted degrees of freedom
      // --------------------------------------------------------------------------------------------  
      size_t const nrDofs = expandedIndices.size();
      // --------------------------------------------------------------------------------------------  
      // Loop over collocation nodes (excluding the initial point) and compute the right hand sides in ru
      // --------------------------------------------------------------------------------------------  
      assemblyRhsTimer.resume();
      
      computeRHS_BDDC_split(steps,x,
                      Cellfltr,
                      number_cells,
                      collocationU,
                      size_e,
                      grid,sweep,t,
                      expandedIndices,eq,eqSemi,
                      assembler,options,ru_bddc,
                      expandedIndices_bddc,
                      uAll, dt,
                      map_II, 
                      map_GammaGamma_noDuplicate, 
                      n_subdomains, 
                      sharedDofsKaskade, 
                      weights,
                      sequenceOfTags,
                      IG_seq,
                      nDofs,
                      ms_double_ass,
                      local2Global,
                      global2Local,
                      map_indices);

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
      // compute restricted rhs contributions
      // -------------------------------------------------------------------------------------------- 
      // compute restricted rhs contributions. Note that we need the WHOLE residual, hence working with the restricted
      // matMuu would not work. Hence we implement the matrix-vector multiplication (with all columns but a row subset)
      // on our own: M*(u_i-u_{i+1})
      std::vector<std::vector<Vector>> matMdiffu_bddc;
      std::vector<std::vector<Vector>> initial_bddc;
      std::vector<std::vector<Vector>> duVec_bddc;
      matMdiffu_bddc.resize(n_subdomains);
      initial_bddc.resize(n_subdomains);
      duVec_bddc.resize(n_subdomains);
      for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
      {
        int tag =  sequenceOfTags[subIndx];
        std::vector<Vector> matMdiffu_sub(grid.points().N(),Vector(expandedIndices_bddc[subIndx].size()));//local2Global[subIndx].size())); // 
        duVec_bddc[subIndx] = matMdiffu_sub;
        initial_bddc[subIndx] = matMdiffu_sub;
        for (int j=0; j<expandedIndices_bddc[subIndx].size(); ++j)
        {
          for (int i=0; i<collocationU.size()-1; ++i)
          {
            matMdiffu_sub[i][j] = 0;
          }

          auto row = Ms[subIndx][expandedIndices_bddc[subIndx][j]];
          for (auto ci=row.begin(); ci!=row.end(); ++ci)
          {
            for (int i=0; i<collocationU.size()-1; ++i)
            {
              matMdiffu_sub[i][j] += *ci * ((collocationU[i].coefficients()[local2Global[tag][ci.index()]])       
                                           -(collocationU[i+1].coefficients()[local2Global[tag][ci.index()]]));  
            }
          }
        }
        matMdiffu_bddc[subIndx] = matMdiffu_sub;
      }
      // -------------------------------------------------------------------------------------------- 
      // compute restricted matrices
      // -------------------------------------------------------------------------------------------- 
      std::vector<SparseMatrix> matMuu_bddc;
      matMuu_bddc.resize(n_subdomains);

      std::vector<SparseMatrix> matStiffuu_bddc;
      matStiffuu_bddc.resize(n_subdomains);

      for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
      {
        SparseMatrix matMuu_sub(expandedIndices_bddc[subIndx],compressedIndex_bddc[subIndx],Ms[subIndx]);
        SparseMatrix matStiffuu_sub(expandedIndices_bddc[subIndx],compressedIndex_bddc[subIndx],Ks[subIndx]);
        matMuu_bddc[subIndx] = matMuu_sub;
        matStiffuu_sub *=(-1.0/options.dt); // cancle out the dt coeffient -dt
        matStiffuu_bddc[subIndx] = matStiffuu_sub;
      }

      // -------------------------------------------------------------------------------------------- 
      // Reaction matrix f_u. Only the restricted subgrid dofs are considered
      // -------------------------------------------------------------------------------------------- 
      assemblyReactionTimer.resume();
      typedef Dune::BDMatrix<typename SparseMatrix::block_type> DiagonalMatrix;
      std::vector<std::vector<DiagonalMatrix>> allMatFu_bddc(n_subdomains);
      for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
      {
        std::vector<DiagonalMatrix> allMatFu_sub(grid.points().N(),DiagonalMatrix(expandedIndices_bddc[subIndx].size()));
        allMatFu_bddc[subIndx] = allMatFu_sub;
      }
      assemblyReactionTimer.stop();
      // -------------------------------------------------------------------------------------------- 
      // perform SDC sweep
      // -------------------------------------------------------------------------------------------- 
      // writeVTKFile(x,output + "/x-step-"+paddedString(steps));
      // printuAll(x,uAll,options, output + "/x-step-"+paddedString(steps),"u_pre");
      std::string name = "test_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep);
      sweepNorm_bddc.push_back( sdcIterationStepBDDCSPLIT(BDDC_SDC_with_initial, BDDC_verbose, BDDC_SDC_verbose, sweep, cg_solver,iter_cg_with_bddc,tol,matlab_dir,
                                                            ifa,
                                                            n_subdomains,
                                                            grid,
                                                            Shat,
                                                            map_II, 
                                                            map_GammaGamma_noDuplicate,  
                                                            sharedDofsKaskade, 
                                                            weights,
                                                            map_indices,
                                                            A_,
                                                            options.dt,
                                                            expandedIndices,
                                                            map_index_to_subdomain,
                                                            sequenceOfTags,
                                                            local2Global,
                                                            global2Local,
                                                            IG_seq,
                                                            matMuu_bddc,
                                                            matStiffuu_bddc,
                                                            ru_bddc,
                                                            allMatFu_bddc,
                                                            matMdiffu_bddc,
                                                            duVec_bddc,
                                                            initial_bddc));
  //     // std::cout << "--------------------------------------------------------------------------------"<<std::endl;
  //     // std::cout <<"sweepNorm_bddc.back() ===> "<<sweepNorm_bddc.back() <<std::endl;
  //     // -------------------------------------------------------------------------------------------- 
  //     // update solution: add Newton correction
  //     // --------------------------------------------------------------------------------------------
  //     // for (int subIdx=0; subIdx<n_subdomains; ++subIdx)
  //     // {
  //     //   auto ui = subs[subIdx].getSolution();
  //     //   auto dui = ui; subs[subIdx].getCorrection(dui);
  //     //   int subIdx_size = sub_length_var[subIdx];     

  //     //   for (int local = 0; local < subIdx_size; ++local)
  //     //   {
  //     //     double val = component<0>(u).coefficients()[local2Global[subIdx][local]] + ui[local];
  //     //     component<0>(u).coefficients()[local2Global[subIdx][local]] = val;
  //     //   }
  //     // }

  //     // for (int ii=1; ii<collocationUe.size(); ii++) // start loop at 1 as initial values is not changed
  //     // {
  //     //   State du(x);
  //     //   State du_mark(x);
  //     //   du=0;
  //     //   du_mark=0;
  //     //   for (int j=0; j<expandedIndices.size(); ++j)
  //     //   {
  //     //     size_t ej = expandedIndices[j];
  //     //     collocationUe[ii].coefficients()[ej] += duVec[ii][j];
  //     //     at_c<0>(du.data).coefficients()[ej] = duVec[ii][j];
  //     //     if(duVec[ii][j]!=0) at_c<0>(du_mark.data).coefficients()[ej] = 1.0;    
  //     //   }
        
  //     //   // if(steps==std::floor(maxSteps/2)) printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
  //     //   // if(steps==std::floor(maxSteps/2)) printuAll(du_mark,uAll,options.order, output + "/du_mark_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du_mark");
  //     //   if(options.plot)
  //     //     printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
  //     // }

  //     // // --------------------------------------------------------------------------------------------
  //     // //
  //     // // --------------------------------------------------------------------------------------------

      for (int ii=1; ii<collocationU.size(); ii++) // start loop at 1 as initial values is not changed
      {
        State du(x);
        State du_mark(x);
        du=0;
        du_mark=0;
        for (int j=0; j<expandedIndices.size(); ++j)
        {
          size_t ej = expandedIndices[j];
          int subIndx =  map_index_to_subdomain[ej];
          int tag =  sequenceOfTags[subIndx];
          collocationU[ii].coefficients()[ej] += duVec_bddc[subIndx][ii][global2Local[tag][ej]]; // FIX ME sinec the du is shrinked!!!

          // std::cout << "global: " << ej << ":" << "subIndx: "<< subIndx << " =>  " << global2Local[subIndx][ej] << " value " << duVec_bddc[subIndx][ii][global2Local[subIndx][ej]]<< std::endl;
          at_c<0>(du.data).coefficients()[ej] = duVec_bddc[subIndx][ii][j];
          if(duVec_bddc[subIndx][ii][j]!=0) at_c<0>(du_mark.data).coefficients()[ej] = 1.0;    
        }
        
        // if(steps==std::floor(maxSteps/2)) printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
        // if(steps==std::floor(maxSteps/2)) printuAll(du_mark,uAll,options.order, output + "/du_mark_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du_mark");
        if(options.plot)
          printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
      }

  //     // for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
  //     // {
        
  //     //   for (int ii=1; ii<collocationU.size(); ii++) // start loop at 1 as initial values is not changed
  //     //   {
  //     //     State du(x);
  //     //     State du_mark(x);
  //     //     du=0;
  //     //     du_mark=0;

  //     //     for (int j=0; j<expandedIndices_bddc[subIndx].size(); ++j)
  //     //     {
  //     //       size_t ej = expandedIndices_bddc[subIndx][j];
  //     //       int subIndx_size = sub_length_var[subIndx];
  //     //       size_t index_local = global2Local[subIndx][ej];

  //     //       if(index_local < subIndx_size){
  //     //         collocationU[ii].coefficients()[ej] += duVec_bddc[subIndx][ii][index_local];
  //     //         //collocationU[ii].coefficients()[local2Global[subIndx][ej]] += duVec_bddc[subIndx][ii][index_local];
  //     //       }
  //     //     }

  //     //     if(options.plot)
  //     //       printuAll(du,uAll,options.order, output + "/du_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep)+"_col_"+paddedString(ii),"du");
  //     //   }
  //     // }

      // --------------------------------------------------------------------------------------------
      // sdcContraction
      // --------------------------------------------------------------------------------------------      
      State sol_tmp(x); 
      component<0>(sol_tmp) = collocationU.back();
      if(options.plot) printuAll(sol_tmp,uAll,options.order, output + "/newton_update_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep),"u_update");
      // printuAll(sol_tmp,uAll,options.order, output + "/emiSDCBDDC"+paddedString(steps)+"_sweep_"+paddedString(sweep),"u");
      if(steps==std::floor(maxSteps/2) and options.plot) printuAll(sol_tmp,uAll,options.order, output + "/newton_update_steps_"+paddedString(steps)+"_sweep_"+paddedString(sweep),"u_update");
      // --------------------------------------------------------------------------------------------
      // Fix the norm
      // --------------------------------------------------------------------------------------------    
      BlockVectorX tmp(nDofs), sol_tmp_copy(nDofs);
      tmp *= 0; 
      sol_tmp_copy *= 0;
      sol_tmp.write(sol_tmp_copy.begin());
      M_.mv(sol_tmp_copy,tmp);
      normU2 = sol_tmp_copy * tmp;

      if (sweepNorm_bddc.size()>1)
      {
        double c = sweepNorm_bddc.back()/sweepNorm_bddc[sweepNorm_bddc.size()-2];
        sdcContraction = std::sqrt(c*sdcContraction);
      }
      // std::cerr << sweep <<"\t"<< expandedIndices.size()  <<"\t"<< sweepNorm_bddc.back() << "\t" <<std::sqrt(normU2) << "\t"<<sdcContraction << "\t\t" << Cellfltr.get_size()<<"\n";   
      std::cerr << sweep <<"\t"<< sweepNorm_bddc.back() << "\t" <<std::sqrt(normU2) << "\t"<<sdcContraction <<"\n";   

      // // --------------------------------------------------------------------------------------------  
      // // select degrees of freedom to take into account in the next sweep. This is a 
      // // further reduction of the up to now used dofs
      // // --------------------------------------------------------------------------------------------  
      int count = 0;
      int discount = 0;
      if (options.tolSelect > 0)
      {
        std::map<int, int> map_expanded_indices;

        State markSelectedDOF(x);
        markSelectedDOF*=0;

        std::vector<size_t> newExpandedIndices;
        std::vector<size_t> newCompressedIndex;
        std::set<size_t> set_ExpandedIndices;
      
        expandedIndices_pre.assign(expandedIndices.begin(),expandedIndices.end());

        newCompressedIndex.resize(compressedIndex.size(),compressedIndex.size());
        compressedIndex.assign(newCompressedIndex.begin(),newCompressedIndex.end());

        for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
        {
          for (int ii=0; ii<duVec_bddc[subIndx].back().size(); ++ii)
          {
            double duMax = 0;
            for (auto const& duj: duVec_bddc[subIndx])
              duMax = std::max(duMax,std::abs(duj[ii]));

            //if(sdcContraction>=1 || sdcContraction*duMax/(1-sdcContraction) > options.tolSelect)
            if(sdcContraction>1 || sdcContraction*duMax/(1-sdcContraction) > options.tolSelect)
            {
              size_t ej = expandedIndices[ii];
              std::set<int> s = index2IndexsSet[local2Global[subIndx][ej]];
              std::set<int>::iterator it;
              for (it = s.begin(); it != s.end(); ++it) {
                set_ExpandedIndices.insert(*it);
              }
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
          Cellfltr.set_cells(selected_cell_idx);
          at_c<0>(markSelectedDOF.data).coefficients()[ej_next] = 1.0;
          size_e_adaptivity++;
        }
        expandedIndices.assign(newExpandedIndices.begin(),newExpandedIndices.end());



        for (int i = 0; i < expandedIndices.size(); ++i)
        {
          // std::cout << expandedIndices[i] <<std::endl;
          compressedIndex[expandedIndices[i]] = i;
        }
        count = expandedIndices.size();
        discount = compressedIndex.size()-expandedIndices.size();
      
        if(options.plot and compressedIndex.size()!=expandedIndices.size()) printuAll(markSelectedDOF,uAll,options.order, output + "/Selected-steps_"+paddedString(steps)+"-sweep_"+paddedString(sweep),"SelectedDOF");

        expandedIndices_pre.assign(expandedIndices.begin(),expandedIndices.end());


        for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
        {

          std::vector<size_t> newCompressedIndex_bddc;
          newCompressedIndex_bddc.resize(compressedIndex_bddc[subIndx].size(),compressedIndex_bddc[subIndx].size());
          compressedIndex_bddc[subIndx].assign(newCompressedIndex_bddc.begin(),newCompressedIndex_bddc.end());

          std::set<size_t> newExpandedIndices_bddc;

          std::unordered_map<int, int>::iterator it; 
          for (int i = 0; i < expandedIndices.size(); ++i)
          {
            size_t ej = expandedIndices[i];
            it = global2Local[subIndx].find(ej);
            if (it != global2Local[subIndx].end()){
              newExpandedIndices_bddc.insert(global2Local[subIndx][ej]);
            }
          }

          expandedIndices_bddc[subIndx].assign(newExpandedIndices_bddc.begin(),newExpandedIndices_bddc.end());

          for (int i = 0; i < expandedIndices_bddc[subIndx].size(); ++i)
          {
            compressedIndex_bddc[subIndx][expandedIndices_bddc[subIndx][i]] = i;
          }

          // for (int i = 0; i < expandedIndices_bddc[subIndx].size(); ++i)
          // {
          //   std::cout<< i << ": "<< expandedIndices_bddc[subIndx][i]<<std::endl;
          // }
       
          // for (int i = 0; i < compressedIndex_bddc[subIndx].size(); ++i)
          // {
          //   std::cout<< "compressedIndex_bddc" << ": "<< compressedIndex_bddc[subIndx][i]<<std::endl;
          // }

         // std::cout << "==============================================\n"; 
        }
        std::cout << "size_e_adaptivity " << size_e_adaptivity <<std::endl;
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
      //   at_c<0>(spatialError.data) = collocationU.back();

      //   projectHierarchically(variableSet,spatialError);
      //   at_c<0>(spatialError.data) -= collocationU.back();
      //   at_c<0>(xnext.data) = collocationU.back();
        
      //   // perform mesh adaptation
      //   accurate = embeddedErrorEstimator(variableSet,spatialError,xnext,IdentityScaling(),tolX,gridManager,0);
      //   if (!accurate) {
      //     size = variableSet.degreesOfFreedom(0,nvars);
      //     //std::cout << "\t\t\t  accurate = " << accurate << ",   dofs after mesh refinement: " << size << std::endl;
      //     sweepNorm.clear(); // old SDC sweep norm on coarser grid cannot be compared to sweeps on finer grid.
      //     reassemble = true;
          
      //     if (options.rosenbrockRefinementStyle)
      //       for (int i=1; i<collocationU.size();   ++i)
      //       {
      //         collocationU[i] = collocationU[0];
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
        std::vector<StateUe> newUe(collocationU.size()+1,collocationU[0]);
        for (int i=0; i<newUe.size(); ++i)
        {
          newUe[i] = 0;
          for (int j=0; j<collocationU.size(); ++j)
          {
            newUe[i].axpy(p[i][j],collocationU[j]);
          }
        }
        collocationU.swap(newUe);
        for (int subIndx = 0; subIndx < n_subdomains; ++subIndx)
        {
          ru_bddc[subIndx].push_back(ru_bddc[subIndx].front());     
          allMatFu_bddc[subIndx].push_back(allMatFu_bddc[subIndx].front()); 
          matMdiffu_bddc[subIndx].push_back(matMdiffu_bddc[subIndx].front());
          duVec_bddc[subIndx].push_back(duVec_bddc[subIndx].front());
          initial_bddc[subIndx].push_back(initial_bddc[subIndx].front()); 
        }

        std::cerr << "time points in ladder method " << grid.points() << '\n';
      }

      if(sweepNorm_bddc.back() < options.SDC_TOL){
        std::cout << "options.SDC_TOL"<<std::endl;
        // std::cout << "break sweep loop: " << sweepNorm.back()   << " sweepNorm.back()*sdcContraction/(1-sdcContraction):  " << sweepNorm.back()*sdcContraction/(1-sdcContraction) << " sdcContraction/(1-sdcContraction):   " << sdcContraction/(1-sdcContraction) << std::endl;
        break;
      }    

      if(count==0 and discount!=0){
        std::cout << "count==0"<<std::endl;
        break;
      }

    }
    //while ( sweep+1<options.maxSweeps && (sweep+1<options.minSweeps ||  sdcContraction>1 || !accurate) );
    while ( sweep+1<options.maxSweeps && (sweep+1<options.minSweeps ||  sdcContraction>1 || sweepNorm_bddc.back()*sdcContraction/(1-sdcContraction)>options.aTol || !accurate) );
    // while ( sweep+1<options.maxSweeps && (sweep+1<options.minSweeps) );
    

    // --------------------------------------------------------------------------------------------
    // extract final-time value
    // --------------------------------------------------------------------------------------------   
    auto step_test(x);
    component<0>(x) = collocationU.back();
    auto updated_sol(x);
    component<0>(updated_sol) -= component<0>(step_test);
    sol_BDDC_SDC *= 0;
    updated_sol.write(sol_BDDC_SDC.begin());
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
    if(options.plot) printuAll(x,uAll,options.order, output + "/emiSDCBDDC"+paddedString(steps),"u");
    // printuAll(x,uAll,options.order, output + "/emiSDCBDDC"+paddedString(steps),"u");
  }
  auto recordTime2 = high_resolution_clock::now();
  MaxStepsTimer.stop();
  std::cout << " ****************************************************************************************** " <<std::endl;
  duration<double, std::milli> ms_double = recordTime2 - recordTime1;
  std::cout <<"ms_double:\t" << ms_double.count() << ", " << ms_double_ass.count() << " ms\n";
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
  printuAll(x,uAll,options.order, output+"/emiSDCBDDCLast","u");  
 return x;
}

#endif
