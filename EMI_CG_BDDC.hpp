#ifndef INTEGRATE_BDDC_CG_HH
#define INTEGRATE_BDDC_CG_HH

// #include "EMI_write_utility.hpp"

template <class VEntry>
void writeSolution_vec(Dune::BlockVector<VEntry> const& b, std::string const& basename, int precision=16)
{
  std::string fname = basename + ".m";
  std::ofstream f(fname.c_str());
  f.precision(precision);
  
  // Write vector.
  f << "function [b] = " << "sol" << '\n'
  << " b = [\n";
  for (size_t i=0; i<b.N(); ++i)
    f << b[i] << std::endl;
  f << "];\n";
}



template <class Grid, class Functional, class VariableSet, class Spaces, class elementType, class Vector, class Options, class Matrix>
typename VariableSet::VariableSet semiImplicit_CG_BDDC(	GridManager<Grid>& gridManager,
  			                                                Functional& F,
  			                                                VariableSet const& variableSet, 
  			                                                Spaces const& spaces,
  			                                                Grid const& grid, 
                                                        Options const& options,
  			                                                std::string out,
  			                                                bool direct,
  			                                                typename VariableSet::VariableSet u,
  			                                                elementType & uAll,
                                                        Vector & sol_bddc,
                                                        std::vector<std::vector<LocalDof>> sharedDofsKaskade,
                                                        int interfaceTypes,
                                                        std::vector<Matrix> As,
                                                        std::vector<int> sequenceOfTags,          
                                                        std::map<int,std::set<int>> map_II,
                                                        std::map<int,std::set<int>> map_GammaGamma_noDuplicate,  
                                                        std::vector<Vector> weights,
                                                        bool cg_solver,
                                                        int iter_cg_with_bddc,
                                                        std::map<int,std::unordered_map<int, int>> local2Global,
                                                        double tol,
                                                        std::map<int, int> map_t2l,
                                                        std::map<int, int> map_indices,
                                                        bool BDDC_verbose,
                                                        std::map<int,std::vector<int>> IG_seq,
                                                        std::string matlab_dir,
                                                        bool write_to_file
                                                        )
{

  double dt = options.dt;
  int maxSteps = std::floor(options.T/options.dt);
  int order  = options.order;
  std::cerr << "semiImplicitCGJacobi: dt = " << dt << " maxSteps " << maxSteps << "\n";

  using namespace boost::fusion;
  typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Functional> >  SemiLinearization;
	typedef VariationalFunctionalAssembler<SemiLinearization> Assembler;

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

  if(options.plot) writeVTK(uAll,out+"/emiBDDCInitial"+paddedString(0,2),
                    IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

  SemiImplicitEulerStep<Functional>  eq(&F,dt);
  double const tau = dt;
  eq.setTau(tau);
  
  size_t nDofs = variableSet.degreesOfFreedom(0,nvars);

	auto du(u);
  auto step_test(u);
	du *= 0;
  auto u_pre(u);

  auto temp_ui(u);
  auto temp_corr(u);
  auto temp_RawCorr(u);
  auto temp_res(u);
  auto temp_resRes(u);
	
  Vector u_semi(nDofs);

  assembler.assemble(SemiLinearization(eq,u,u,du),options.assemblyThreads);
  AssembledGalerkinOperator<Assembler> A(assembler);
  Matrix LHS = assembler.template get<Matrix>(false); 
  int n_subdomains = sequenceOfTags.size();
  std::vector<int> subdomSize(n_subdomains);
  for (int subIdx=0; subIdx<n_subdomains; ++subIdx){
    int tag = sequenceOfTags[subIdx];
    subdomSize[subIdx] = As[subIdx].N();  
  }

  InterfaceAverages<1,int> ifa(sharedDofsKaskade,subdomSize,interfaceTypes);


  using TransmissionScalar = double;
  using BddcSubdomain = Subdomain<1,double,double,SpaceTransfer<1,double,TransmissionScalar>>;
  //using BddcSubdomain = Subdomain<1>;
  std::vector<std::unique_ptr<BddcSubdomain>> subsptr(n_subdomains);
  std::vector<int> activeIds(n_subdomains);
  std::iota(activeIds.begin(), activeIds.end(), 0);
  // parallelFor(0,n_subdomains,[&](int subIdx)
  // {
  //   int tag = sequenceOfTags[subIdx];
  //   subsptr[subIdx] = std::make_unique<BddcSubdomain>(subIdx,As[subIdx],ifa);
  //   //if(subIdx!=3) activeIds.push_back(subIdx);
  //   //if(subIdx<38) 
  //   activeIds.push_back(subIdx);
  //   // if(subIdx!=3) activeIds.push_back(subIdx);
  // });
  for (int subIdx = 0; subIdx < n_subdomains; ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::cout << subIdx << " -> " << tag << " -> "<< activeIds[subIdx] << std::endl;
    subsptr[subIdx] = std::make_unique<BddcSubdomain>(subIdx,As[subIdx],ifa);
    //activeIds[subIdx] = subIdx;
  }
  
  // parallelFor(0,n_subdomains,[&](int subIdx)
  // {
  //   activeIds.push_back(subIdx);
  // });

  for (int time_step=0; time_step<maxSteps; ++time_step) 
  {
    std::cout << "\n ---------- itr "<< time_step << " ---------- "<< std::endl;
    // set nnz to the number of structural nonzero elements of the matrix to be assembled below
    size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
    size_t  size = variableSet.degreesOfFreedom(0,nvars);

    // ------------------------------------------------------------------------------------
    // update the rhs
    // ------------------------------------------------------------------------------------
    F.Mass_stiff(10);
    timer.start("updating rhs");
    SemiImplicitEulerStep<Functional>  eq(&F,options.dt);
    double const tau = options.dt;
    eq.setTau(tau);

    assembler.assemble(SemiLinearization(eq,u,u,du),Assembler::RHS,options.assemblyThreads);
    timer.stop("updating rhs");

    // ------------------------------------------------------------------------------------
    // update rhs
    // ------------------------------------------------------------------------------------
    timer.start("updating sub rhs");
    
    Vector rhs_vec_test(nDofs);
    auto rhs = assembler.rhs();
    rhs.write(rhs_vec_test.begin());

    Vector rhs_petsc_test(nDofs);
    rhs.write(rhs_petsc_test.begin());

    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, rhs_vec_test, rhs_petsc_test);
    writeSolution_sol(rhs_petsc_test,matlab_dir+"/rhs_bddc_"+std::to_string(time_step)); 
    std::vector<Vector> Fs(n_subdomains);
        // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // fill the sub matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
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
          // std::cout << "("<< index<< ","<< coef << "), ";          
          Fs_subIdx[i] = coef*rhs_petsc_test[index];
        }
        // std::cout << "\n";
      }
      Fs[subIdx] = Fs_subIdx;
    }

     if(write_to_file and false)
    {
      for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
      {
        std::string path = std::to_string(subIdx);
        writeToMatlabPath(As[subIdx],Fs[subIdx],"A_kaskade_shrinked"+subIdx,matlab_dir, true);      
      }  
    }
    
    timer.stop("updating sub rhs");

    timer.start("alg subdom creation");
    
    std::vector<BddcSubdomain> subs;
    for (auto& sp: subsptr){
      subs.push_back(*sp);
    }

    timer.stop("alg subdom creation");
  
    timer.start("bddc creation");
    BDDCSolver<BddcSubdomain> bddcSolver(subs,ifa.coarseConstraints(),activeIds,cg_solver,BDDC_verbose );
    bddcSolver.update_rhs(Fs);
    timer.stop("bddc creation");

    // ------------------------------------------------------------------------------------
    // call bddc solver
    // ------------------------------------------------------------------------------------
    std::cout << "\n\n";
    std::vector<double> resNorm;
    for (int k=0; k<iter_cg_with_bddc; ++k)
    {
      timer.start("BDDC solve");
      resNorm.push_back(bddcSolver.solve());      
      timer.stop("BDDC solve");
      if(false)
      {
        for (int subIdx=0; subIdx<n_subdomains; ++subIdx)
        {
          int tag = sequenceOfTags[subIdx];
          auto ui = subs[subIdx].getSolution();
          auto corr = ui; subs[subIdx].getCorrection(corr);
          auto rawCorr = ui; subs[subIdx].getRawCorrection(rawCorr);
          auto res = ui; subs[subIdx].getResidual(res);
          auto resRes = ui; resRes = subs[subIdx].getRestrictedResidual();


          int subIdx_size = map_t2l[tag];     
          for (int local = 0; local < subIdx_size; ++local)
          {
            component<0>(temp_ui).coefficients()[local2Global[tag][local]] = ui[local];
            component<0>(temp_corr).coefficients()[local2Global[tag][local]] = corr[local];
            component<0>(temp_RawCorr).coefficients()[local2Global[tag][local]] = rawCorr[local];
            component<0>(temp_res).coefficients()[local2Global[tag][local]] = res[local];
            component<0>(temp_resRes).coefficients()[local2Global[tag][local]] = resRes[local];
          }
        }
        uAll = component<0>(temp_ui);
        writeVTK(uAll,out+"/emiBDDC_Sol"+paddedString(time_step,2)+"_itr_"+paddedString(k,2),
         IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

        uAll = component<0>(temp_corr); 
        writeVTK(uAll,out+"/emiBDDC_Corr"+paddedString(time_step,2)+"_itr_"+paddedString(k,2),
         IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

        uAll = component<0>(temp_RawCorr); 
        writeVTK(uAll,out+"/emiBDDC_RawCorr"+paddedString(time_step,2)+"_itr_"+paddedString(k,2),
         IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

        uAll = component<0>(temp_res);
        writeVTK(uAll,out+"/emiBDDC_Res"+paddedString(time_step,2)+"_itr_"+paddedString(k,2),
         IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

        uAll = component<0>(temp_resRes);
        writeVTK(uAll,out+"/emiBDDC_resRes"+paddedString(time_step,2)+"_itr_"+paddedString(k,2),
         IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");        
      }


      if(resNorm.back()<tol)
        break;
    } 
    // ------------------------------------------------------------------------------------
    // update solution
    // ------------------------------------------------------------------------------------
    {
      for (int subIdx=0; subIdx<n_subdomains; ++subIdx)
      {
        int tag = sequenceOfTags[subIdx];
        auto ui = subs[subIdx].getSolution();
        auto dui = ui; subs[subIdx].getCorrection(dui);
        int subIdx_size = map_t2l[tag];     
        for (int local = 0; local < subIdx_size; ++local)
        {
          double val = component<0>(u).coefficients()[local2Global[tag][local]] + ui[local];
          component<0>(u).coefficients()[local2Global[tag][local]] = val;
          component<0>(step_test).coefficients()[local2Global[tag][local]] = ui[local];
        }
      }
    }

    sol_bddc *= 0;
    step_test.write(sol_bddc.begin());

    uAll = component<0>(u);

    Vector sol_petsc_test(nDofs);
    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_bddc, sol_petsc_test);
    writeSolution_sol(sol_petsc_test,matlab_dir+"/sol_bddc_"+std::to_string(time_step)); 

    if(options.plot) writeVTK(uAll,out+"/emiBDDC"+paddedString(time_step,2),
             IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

    writeVTK(uAll,out+"/emiBDDC"+paddedString(time_step,2),
             IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
    int lookback = std::min(10,iter_cg_with_bddc-1);
    double contraction = std::pow(resNorm.back()/resNorm[resNorm.size()-lookback],1.0/lookback);
    std::cout << "Estimated contraction factor: " << contraction << ". (kappa ~ " << (1+contraction)/(1-contraction) << ").\n";
  }
  writeVTK(uAll,out+"/emiBDDCLast",
             IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
	return u;
}
#endif