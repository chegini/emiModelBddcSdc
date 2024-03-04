#ifndef INTEGRATE_CG_HH
#define INTEGRATE_CG_HH

template <class Grid, class Functional, class VariableSet, class Spaces, class elementType, class Vector, class Options, class Matrix>
typename VariableSet::VariableSet semiImplicit_CG_BDDC(	GridManager<Grid>& gridManager,
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
                                                        Vector & sol_semi,
                                                        std::vector<std::vector<LocalDof>> sharedDofsKaskade,
                                                        int interfaceTypes,
                                                        int n_subdomains,
                                                        std::vector<Matrix> As,
                                                        std::map<int,std::set<int>> IGamma,
                                                        std::vector<int> sequenceOfTags, 
                                                        std::vector<int> startingIndexOfTag,           
                                                        std::map<int,std::set<int>> map_II,
                                                        std::map<int,std::set<int>> map_GammaGamma,  
                                                        std::map<int,std::map<int,std::set<int>>> map_GammaNbr, 
                                                        std::map<int,std::vector<int>> sequanceOfsubdomains,
                                                        std::vector<Vector> weights,
                                                        std::vector<Vector> Fs,
                                                        bool cg_solver,
                                                        int iter_cg_with_bddc,
                                                        std::vector<std::unordered_map<int, int>> local2Global,
                                                        std::vector<std::unordered_map<int, int>> global2Local,
                                                        double tol,
                                                        std::map<int, int> map_t2l,
                                                        std::map<int, int> map_indices,
                                                        bool BDDC_verbose
                                                        )
{

  // double dt = options.dt;
  // int maxSteps = std::floor(options.T/options.dt);
  // int order  = options.order;
  // std::cerr << "semiImplicitCGJacobi: dt = " << dt << " maxSteps " << maxSteps << "\n";

  // using namespace boost::fusion;
  // typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Functional> >  SemiLinearization;
	// typedef VariationalFunctionalAssembler<SemiLinearization> Assembler;

  // int number_cells = gridManager.grid().size(0);

	// auto& timer = Timings::instance();
	// Assembler assembler(spaces);

	// constexpr int nvars = Functional::AnsatzVars::noOfVariables;
	// constexpr int neq = Functional::TestVars::noOfVariables;

  // typedef typename Functional::OriginVars::template CoefficientVectorRepresentation<0,neq>::type LinearSpace;

  // size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
  // size_t  size = variableSet.degreesOfFreedom(0,nvars);
  // std::cout << " size  "<< size <<"   nnz  " << nnz << std::endl;

  // uAll = component<0>(u);

  // if(options.plot) writeVTK(uAll,out+"/emiBDDCInitial"+paddedString(0,2),
  //                   IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

  // SemiImplicitEulerStep<Functional>  eq(&F,dt);
  // double const tau = dt;
  // eq.setTau(tau);
  
  // size_t nDofs = variableSet.degreesOfFreedom(0,nvars);

	// auto du(u);
  // auto step_test(u);
	// du *= 0;
  // auto du_cg(u);
  // du_cg *= 0;

	// Vector u_semi(nDofs);

  // assembler.assemble(SemiLinearization(eq,u,u,du),options.assemblyThreads);
  // AssembledGalerkinOperator<Assembler> A(assembler);
  // Matrix LHS = assembler.template get<Matrix>(false); 
  // std::vector<int> subdomSize(n_subdomains);
  // std::vector<int> subdomLocalSize(n_subdomains);
  // for (int subIdx=0; subIdx<n_subdomains; ++subIdx){
  //   subdomSize[subIdx] = As[subIdx].N();  
  //   subdomLocalSize[subIdx] = IGamma[subIdx].size();
  // }

  // InterfaceAverages<1,int> ifa(sharedDofsKaskade,subdomSize,interfaceTypes);
  // InterfaceAverages<1,int> ifa(sharedDofsKaskade,subdomSize,interfaceTypes);

  // for (int time_step=0; time_step<maxSteps; ++time_step) 
  // {
  //   std::cout << "\n ---------- itr "<< time_step << " ---------- "<< std::endl;
  //   // set nnz to the number of structural nonzero elements of the matrix to be assembled below
  //   size_t  nnz = assembler.nnz(0,neq,0,nvars,false);
  //   size_t  size = variableSet.degreesOfFreedom(0,nvars);

  //   // ------------------------------------------------------------------------------------
  //   // update the rhs
  //   // ------------------------------------------------------------------------------------
  //   F.Mass_stiff(10);
  //   timer.start("updating rhs");
  //   SemiImplicitEulerStep<Functional>  eq(&F,options.dt);
  //   double const tau = options.dt;
  //   eq.setTau(tau);

  //   assembler.assemble(SemiLinearization(eq,u,u,du),Assembler::RHS,options.assemblyThreads);
  //   timer.stop("updating rhs");

  //   // ------------------------------------------------------------------------------------
  //   // update rhs
  //   // ------------------------------------------------------------------------------------
  //   timer.start("updating sub rhs");
    
  //   Vector rhs_vec_test(nDofs);
  //   auto rhs = assembler.rhs();
  //   rhs.write(rhs_vec_test.begin());

  //   std::vector<Vector> Fs;
  //   construct_Fs(sequenceOfTags, startingIndexOfTag, map_II, map_GammaGamma, map_GammaNbr, sequanceOfsubdomains, rhs_vec_test_new, sharedDofsKaskade, weights, map_indices, Fs);
  //   timer.stop("updating sub rhs");

  //   timer.start("alg subdom creation");
  //   using TransmissionScalar = double;
  //   using BddcSubdomain = Subdomain<1,double,double,SpaceTransfer<1,double,TransmissionScalar>>;
  //   //using BddcSubdomain = Subdomain<1>;
  //   std::vector<std::unique_ptr<BddcSubdomain>> subsptr(n_subdomains);
    
  //   // for (int i = 0; i < n_subdomains; ++i)
  //   // {
  //   //   std::cout << i << " " << As[i].N() << " " << Fs[i].size()<<std::endl;
  //   // }
 
  //   parallelFor(0,n_subdomains,[&](int i)
  //   {
  //     subsptr[i] = std::make_unique<BddcSubdomain>(i,As[i],Fs[i],ifa);
  //     // subsptr[i] = std::make_unique<Subdomain<1>>(i,As[i],Fs[i],ifa);
  //   });

  //   std::vector<BddcSubdomain> subs;
  //   for (auto& sp: subsptr)
  //     subs.push_back(std::move(*sp));

  //   timer.stop("alg subdom creation");
  
  //   timer.start("bddc creation");
  //   BDDCSolver<BddcSubdomain> bddcSolver(subs,ifa.coarseConstraints(),cg_solver,BDDC_verbose );
  //   // BDDCSolver<Subdomain<1>> bddcSolver(subs,ifa.coarseConstraints());
  //   timer.stop("bddc creation");

  //   // ------------------------------------------------------------------------------------
  //   // call bddc solver
  //   // ------------------------------------------------------------------------------------
  //   std::cout << "\n\n";
  //   std::vector<double> resNorm;
  //   for (int k=0; k<iter_cg_with_bddc; ++k)
  //   {
  //     timer.start("BDDC solve");
  //     resNorm.push_back(bddcSolver.solve());      
  //     timer.stop("BDDC solve");
  //     if(resNorm.back()<tol)
  //       break;
  //   } 
  //   // ------------------------------------------------------------------------------------
  //   // update solution
  //   // ------------------------------------------------------------------------------------
  //   {
  //     for (int subIdx=0; subIdx<n_subdomains; ++subIdx)
  //     {
  //       auto ui = subs[subIdx].getSolution();
  //       auto dui = ui; subs[subIdx].getCorrection(dui);
  //       int subIdx_size = sub_length_var[subIdx];     

  //       for (int local = 0; local < subIdx_size; ++local)
  //       {
  //         double val = component<0>(u).coefficients()[local2Global[subIdx][local]] + ui[local];
  //         component<0>(u).coefficients()[local2Global[subIdx][local]] = val;
  //       }
  //     }
  //   }

  //   uAll = component<0>(u);

  //   if(options.plot) writeVTK(uAll,out+"/emiBDDC"+paddedString(time_step,2),
  //            IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

  //   writeVTK(uAll,out+"/emiBDDC"+paddedString(time_step,2),
  //            IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
  //   int lookback = std::min(10,iter_cg_with_bddc-1);
  //   double contraction = std::pow(resNorm.back()/resNorm[resNorm.size()-lookback],1.0/lookback);
  //   std::cout << "Estimated contraction factor: " << contraction << ". (kappa ~ " << (1+contraction)/(1-contraction) << ").\n";
  // }
  // writeVTK(uAll,out+"/emiBDDCLast",
  //            IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
  // writeVectorTofile(u,"matlab_dir/emiBDDCLast");
	return u;
}
#endif