/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2021-2021 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* author: Fatemeh chegini
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "EMI.hh"
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
int main(int argc, char* argv[])
{
  using namespace Kaskade::BDDC;
  using namespace boost::fusion;

  std::cout << "Start subdomain tutorial program" << std::endl;

  constexpr int dim = SPACEDIM;
  int refinements, order, solver, refinements_sol, interfaceTypes, iter_cg_with_bddc;
  double penalty, sigma_i, sigma_e, C_m, R, R_extra, tol, dt;
  bool  direct, onlyLowerTriangle, vtk_, timing, test_newCof;
  bool run_implicit_CG, run_implicit_CG_SDC, run_implicit_CG_BDDC, run_implicit_CG_SDC_BDDC, run_implicit_CG_SDC_BDDC_first_Sweep;
  std::string inputfile, early_excited, extra_set, intra_set, dir_out, matlab_dir;
  bool cg_semi, plot, withSplitFace,cg_solver;
  bool test_mesh_data, write_to_file;
  bool BDDC_SDC_with_initial, BDDC_verbose;
  int verbose, assemblyThreads;
  CardiacIntegrationOptions options;
  if (getKaskadeOptions(argc,argv,Options
  ("input",                    inputfile,                           "./input/example4subc_mesh.vtu","subdomain definition")
  ("extra_set",                extra_set,                           "./input/example4subc_list_extracellular.txt","subdomain definition")
  ("intra_set",                intra_set,                           "./input/example4subc_list_intracellular.txt","subdomain definition")
  ("excited",                  early_excited,                       "./input/example4subc_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/example4subc_2extra_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/example4subc_2extra_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/example4subc_2extra_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/example4subc_2extra_early_excited.txt","subdomain definition")
  ("dir",                      dir_out,                             "./output","subdomain definition")
  ("matlab_dir",               matlab_dir,                          "./matlab_dir","subdomain definition")
  ("refine",                   refinements,                         0,"uniform mesh refinements")
  ("refine_sol",               refinements_sol,                     0,"uniform mesh refinements")
  ("order",                    order,                               1,"polynomial ansatz order")
  ("penalty",                  penalty,                             1e6,"penalty factor for Dirichlet/Nitsche boundary conditions")
  ("sigma_i",                  sigma_i,                             3.0,"intracellular conductivity")
  ("sigma_e",                  sigma_e,                             20.0,"extracellular conductivity")
  ("C_m",                      C_m,                                 1.0,"membrane capacitance")
  ("R",                        R,                                   0.1,"the conductances of the gap junction") 
  ("R_extra",                  R_extra,                             1e-7,"the conductances between two extracellulatr regions")
  ("direct",                   direct,                              false,"if true, use a direct solver")
  ("cg",                       cg_semi,                             true, "cg")
  ("solver",                   solver,                              2,"0=UMFPACK, 1=PARDISO 2=MUMPS 3=SUPERLU 4=UMFPACK32/64 5=UMFPACK64")
  ("onlyLowerTriangle",        onlyLowerTriangle,                   true, "onlyLowerTriangle")
  ("interfacetypes",           interfaceTypes,                      7,"bit flags for coarse interfaces to include: 1 corner 2 edge 3 face")
  ("vtk",                      vtk_,                                true,"write VTK output")
  ("iter_cg_bddc",             iter_cg_with_bddc,                   1000,"number of BDDC iterations")
  ("timing",                   timing,                              true,"whether to write timing info")
  ("test",                     test_mesh_data,                      false,"debug mode")
  ("run_implicit_CG",          run_implicit_CG,                     true, "run linearly semi-implicit method + CG")
  ("run_implicit_CG_SDC",      run_implicit_CG_SDC,                 false, "run linearly semi-implicit method + CG + Jacobi + SDC")
  ("run_implicit_CG_BDDC",     run_implicit_CG_BDDC,                true, "run linearly semi-implicit method + CG + Jacobi + SDC ")
  ("run_implicit_CG_SDC_BDDC", run_implicit_CG_SDC_BDDC,            false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("run_implicit_CG_SDC_BDDC_first_Sweep", run_implicit_CG_SDC_BDDC_first_Sweep,false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("test_newCof",              test_newCof,                         false,"to test the coefficients")
  ("withSplitFace",            withSplitFace,                       false,"split faces in BDDC")  
  ("cg_solver",                cg_solver,                           true,"split faces in BDDC")  
  ("write_to_file",            write_to_file,                       false,"write to matlab file")  
  ("maxSteps",                 options.maxSteps,                    10,  "max number of time steps")
  ("vtk",                      options.writeVTK,                    1,  "write VTK output files 0=none, 1=time steps 2=sweeps")
  ("T_",                       options.T,                           0.01,  "final time[ms]")
  ("dt",                       options.dt,                          0.01,  "time step size[ms]")
  ("orderU",                   options.order,                       1,  "FE ansatz order for transmembrane voltage & action potential")
  ("atol",                     options.aTol,                        1e-15,  "absolute L^2 tolerance")
  ("stol",                     options.tolSelect,                   0.0,  "L^inf tol for DoF selection")
  ("maxCGIter",                options.maxCGIter,                   10000,  "maximum number of IterateType::CG iterations in linear solver (0=direct solver)")
  ("cgTol",                    options.cgTol,                       1e-8,  "absolute IterateType::CG energy error tolerance")
  ("adapt",                    options.adapt,                       false,  "do adaptivity or not")
  ("sweeps",                   options.minSweeps,                   1000,  "minimal number of SDC sweeps")
  ("maxSweeps",                options.maxSweeps,                   1000,  "maximal number of SDC sweeps")
  ("nColloc",                  options.nCollocU,                    1,  "number of collocation points in time")
  ("nCollocStart",             options.nCollocUstart,               1,  "start sweeps with that many collocation points")
  ("verbose",                  options.verbosity,                   1,  "output density")
  ("sweepType",                options.sweepType,                   1,  "0: Euler, 1: LU")
  ("nReactionSweeps",          options.nReactionSweeps,             0,  "number of post-sweep Euler steps for reaction nonlinearity")
  ("minRefLevel",              options.minRefLevel,                 0,  "keep this level on refinement") 
  ("nThreads",                 options.assemblyThreads,             28,  "# of threads in assembler") 
  ("rosenbrockRefinement",     options.rosenbrockRefinementStyle,   false,  "do spatial refinement first, before SDC") 
  ("plot",                     options.plot,                        false,  "final time[ms]")
  ("CG_shift",                 options.CG_shift,                    true,  "shift the cg update")
  ("SDC_TOL",                  options.SDC_TOL,                     1e-5,  "SDC tolerance")
  ("sdc_contraction",          options.sdc_contraction,             0.2,  "SDC constraction")  
  ("BDDC_SDC_with_initial",    BDDC_SDC_with_initial,               true,  "BDDC with initial guess from previous collocation sol") 
  ("BDDC_verbose",             BDDC_verbose,                        true,  "SDC tolerance") 
  )) return 0;
  tol =  options.cgTol;
  if (mkdir("output", 0777) == -1)
    std::cerr << "Creating directory 'output':  " << strerror(errno) << std::endl;
  else
    std::cout << "Directory called 'output' created" << std::endl;

  if (mkdir("matlab_dir", 0777) == -1)
    std::cerr << "Creating directory 'matlab_dir':  " << strerror(errno) << std::endl;
  else
    std::cout << "Directory called 'output' created" << std::endl;

  std::string out = dir_out;

  std::cout << "dt: " << options.dt <<std::endl;
  std::cout << "order: " << order <<std::endl;
  std::cout << "refinement: "<< refinements << std::endl;
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;

  auto& timer = Timings::instance();

  // ------------------------------------------------------------------------------------------------------------
  // get the material from the mesh  
  // ------------------------------------------------------------------------------------------------------------
  timer.start("reading FE grid");
  VTKReader vtk(inputfile);
  GridManager<Grid> gridManager(vtk.createGrid<Grid>());
  using H1SpaceMaterial = FEFunctionSpace<DiscontinuousLagrangeMapper<double,LeafView>>;
  using cellMaterial = H1SpaceMaterial::Element<1>::type;
  H1SpaceMaterial materialSpace(gridManager,gridManager.grid().leafGridView(), 0);
  cellMaterial material(materialSpace);
  vtk.getCoefficients("domain",material);

  if(refinements) std::cout << "---------- refinements ---------- "<< std::endl;
  gridManager.globalRefine(refinements);
  writeVTK(material,out+"/Materials",IoOptions(),"conductivity");
  timer.stop("reading FE grid");

  // ------------------------------------------------------------------------------------------------------------
  // the membrane model & stress model
  // ------------------------------------------------------------------------------------------------------------
  // typedef TenTusscher Membrane;
  typedef AlievPanfilov Membrane;
  Membrane membrane;
  // std::cout << "Using membrane model " << membrane.name() << ".\n";

  // ------------------------------------------------------------------------------------------------------------
  // read extracellular materials
  // ------------------------------------------------------------------------------------------------------------
  std::ifstream file_extra_list(extra_set);
  int n_extra_set;
  file_extra_list >> n_extra_set;
  std::vector<int> arr_extra_set(n_extra_set);
  getSubdomain(arr_extra_set, file_extra_list);

  std::ifstream file_intra_list(intra_set);
  int n_intra_set;
  file_intra_list >> n_intra_set;
  std::vector<int> arr_intra_set(n_intra_set);
  getSubdomain(arr_intra_set, file_intra_list);

  std::ifstream file_excited_region(early_excited);
  int n_excited_region;
  file_excited_region >> n_excited_region;
  std::vector<int> arr_excited_region(n_excited_region);
  getSubdomain(arr_excited_region, file_excited_region);

  Dune::FieldVector<double,dim> zero(0);


  // assume that only extra cellular has the tag number zero
  // in case of having more than one extracellular subdomain, 
  // for the aszast space, we consider only one extra cellular subdoamjn, thats' why we 
  // map all the extracellular tags to zero,
  // however in order to assemble matrices for BDDC, we need to consider the origonal tags numbers.  
  std::map<int, Dune::FieldVector<double,1>> map_OriginalTag_anzastTag;
  for (int i = 0; i < arr_extra_set.size(); ++i)
  {
    Dune::FieldVector<double,1> extra_zero(0);
    extra_zero[0] = 0.0;
    map_OriginalTag_anzastTag.insert({ arr_extra_set[i], extra_zero }); 
  }
  for (int i = 0; i < arr_intra_set.size(); ++i)
  {
    Dune::FieldVector<double,1> intra_zero(0);
    intra_zero[0] = arr_intra_set[i];
    map_OriginalTag_anzastTag.insert({ arr_intra_set[i], intra_zero }); 
  }


  FEFunctionSpace uSpace(gridManager,PiecewiseContinuousLagrangeMapper( gridManager.grid().leafGridView(),
                                                                        order,
                                                                        [&](auto cell) { 
                                                                            if(arr_extra_set.size()>1){
                                                                              return map_OriginalTag_anzastTag[static_cast<int>(material.value(cell,zero))];
                                                                            }
                                                                            return material.value(cell,zero);
                                                                          }));
  

  L2Space<Grid> outSpace(gridManager,gridManager.grid().leafGridView(),order);

  auto spaces = makeSpaceList(&uSpace);
  auto variableSetDesc = makeVariableSetDescription(spaces,
                                  boost::fusion::make_vector(Variable<SpaceIndex<0>,Components<1>,VariableId<0>>("u")));

  L2Space<Grid>::Element_t<1> uAll(outSpace);

  using VariableSetDesc = decltype(variableSetDesc);
  using GridType = decltype(gridManager.grid());
  using SPACE = decltype(spaces);
  using Functional = EMI_model<double,VariableSetDesc,cellMaterial,GridType,SPACE,Membrane>;

  Functional F( material,
                gridManager.grid(),
                spaces,
                penalty,
                sigma_i,
                sigma_e,
                C_m,  
                R,
                R_extra);
  F.extracellular_materials(arr_extra_set);
  typedef typename Functional::OriginVars::VariableSet State;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  constexpr int neq = Functional::TestVars::noOfVariables;
  std::cout << " neq "<< neq <<std::endl;
  std::cout << " nvars "<< nvars <<std::endl;

  using LinearSpace = VariableSetDesc::CoefficientVectorRepresentation<0,neq>::type;

  //construct Galerkin representation
  using SemiLinearization = SemiLinearizationAtInner<SemiImplicitEulerStep<Functional>>;
  using Assembler = VariationalFunctionalAssembler<SemiLinearization>;
  using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
  typedef Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>> Matrix;
  size_t nDofs = variableSetDesc.degreesOfFreedom(0,nvars);

  Assembler assembler(spaces);

  auto up = variableSetDesc.variableSet();
  auto u = variableSetDesc.variableSet();
  auto uM = variableSetDesc.variableSet();
  std::cout << " The numeber of dofs: "<< nDofs <<std::endl;
  std::cout << " The numeber of cells :: "<<gridManager.grid().size(0)<<std::endl;
  size_t dof_size = variableSetDesc.degreesOfFreedom(0, 1);
   // ------------------------------------------------------------------------------------
  // initila the data
  // ------------------------------------------------------------------------------------
  F.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);

  dt = options.dt;
  SemiImplicitEulerStep<Functional>  eq(&F,dt);
  double const tau = dt;
  eq.setTau(tau);
  auto du(u);
  du *= 0;
  assembler.assemble(SemiLinearization(eq,u,u,du),Assembler::RHS,options.assemblyThreads);
  auto rhs_oiginal = assembler.rhs();
  Vector rhs_vec_original(nDofs);
  rhs_oiginal.write(rhs_vec_original.begin());


  // ------------------------------------------------------------------------------------
  // semi implicit + CG methods
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG){
      Vector sol_semi(nDofs);
      Functional F_semi( material,
                  gridManager.grid(),
                  spaces,
                  penalty,
                  sigma_i,
                  sigma_e,
                  C_m,  
                  R,
                  R_extra);
      F_semi.extracellular_materials(arr_extra_set);
      F_semi.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      timer.start("linearly semi implicit method");
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implict approach" << std::endl;
      std::cout << "---------------------------------------------" << std::endl;


      uAll = component<0>(u);
      u = semiImplicit_CG_Jacobi( gridManager,
                                  F_semi,
                                  variableSetDesc,
                                  spaces,
                                  gridManager.grid(),
                                  options,
                                  out,
                                  cg_semi, 
                                  direct,
                                  u,
                                  uAll,
                                  sol_semi
                                  );  
      timer.stop("linearly semi implicit method");

      {
        Vector sol_semi_to_petsc(sol_semi);
        sol_semi_to_petsc = 0;
        // petsc_structure_rhs(dof_set,II, GammaGamma, sol_semi, n_subdomains,sol_semi_to_petsc);
        if(write_to_file) writeSolution(sol_semi_to_petsc,matlab_dir+"/sol");
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // Extract the mesh data
  // - II, GammaGamma, IGamma, GammaGamma_W_Nbr, gamma_nbrs, sequanceOfsubdomains 
  // - e2i
  // - i2e
  // - i2i
  // - icoord
  // - itT
  // - map_t2l
  // - map_sT2l
  // - map_II
  // - map_IGamma
  // - map_GammaGamma
  // ------------------------------------------------------------------------------------

  std::vector<std::vector<int>> e2i(gridManager.grid().size(0)); //element to indices
  std::vector<std::set<int>> e2e(gridManager.grid().size(0));    //element to element
  std::vector<std::set<int>> i2e(dof_size);                      //index to element, for the cell Filter
  std::vector<std::set<int>> i2i(dof_size);                      //index to index

  std::vector<std::vector<double>> icoord(dof_size);             //coordinates of each dofs
  std::vector<int> i2T(dof_size);                                //index to tags
  std::map<int, int> map_t2l;                                    //map: tag to lenth
  std::map<int, int> map_sT2l;                                   //map: sequance of each tag to length
  std::map<int,int> map_nT2oT;                                   //map: new Tag to original Tag
                              
  std::map<int,std::set<int>> map_II;                            // II
  std::map<int,std::set<int>> map_IGamma;                        // IGamma
  std::map<int,std::set<int>> map_GammaGamma;                    // GammaGamma
  std::map<int,std::set<int>> map_GammaGamma_W_Nbr;              // GammaGamma with nbr
  std::map<int,std::map<int,std::set<int>>> map_GammaNbr;        // GammaNbr
  std::map<int,std::vector<int>> sequanceOfsubdomains;           // sequence of neighboring tags for each subdomain

  mesh_data_structure(boost::fusion::at_c<0>(u.data),  
                      material, 
                      e2i, i2e, e2e, i2i, icoord, i2T, map_t2l, map_sT2l, map_II, map_IGamma, map_GammaGamma, map_GammaGamma_W_Nbr, map_GammaNbr);
  int n_subdomains = map_t2l.size();

  std::vector<int> sequenceOfTags(n_subdomains);
  std::vector<int> startingIndexOfTag(n_subdomains);

  computed_sequenceOfTags(map_t2l,map_GammaNbr, sequenceOfTags, startingIndexOfTag, map_nT2oT, sequanceOfsubdomains);

  // ------------------------------------------------------------------------------------
  // compute the data petsc from the mesh data
  // - local2Global
  // - global2Local
  // - globalIndices
  // - map_indices
  // - map_i2sub
  // - i2iSet
  // ------------------------------------------------------------------------------------
  std::map<int,std::unordered_map<int, int>> local2Global;
  std::map<int,std::unordered_map<int, int>> global2Local;
  std::map<int,std::vector<int>> globalIndices;
  subdomain_indices(sequenceOfTags, map_II, map_GammaGamma, map_GammaNbr, local2Global, global2Local, globalIndices);

  std::map<int, int> map_indices;
  std::map<int, int> map_i2sub;
  map_kaskade2petcs(sequenceOfTags, map_II, map_GammaGamma, map_indices, map_i2sub);
  
  removeInnerIndices_i2i(i2i);
  std::set<std::set<int>> i2iSet(i2i.begin(),i2i.end());  //index to index only those has more than one neighours on the interfaces
  if(false)
  {
    std::set<std::set<int>>::iterator it;
    for (it = i2iSet.begin(); it != i2iSet.end(); ++it) {
      std::set<int> s = *it;
      std::set<int>::iterator itr;
      for (itr = s.begin(); itr != s.end(); ++itr) {
        std::cout << *itr << " ";
      }
      std::cout <<std::endl;
    }
  }

  std::vector<std::vector<LocalDof>> sharedDofsKaskade;
  compute_sharedDofsKaskade(sequenceOfTags, map_indices, map_II, map_GammaGamma, map_GammaNbr, write_to_file, matlab_dir, sharedDofsKaskade);
  std::cout << "generated sub matrices of cell by cell for BDDC in petsc(data for Kaskade)~!!!!!\n\n\n\n" << std::endl;

  int mesh_dim = SPACEDIM==2? 2:3;
  write_Dirichlet_and_coordinates(boost::fusion::at_c<0>(u.data), e2i, map_indices, icoord, dof_size, mesh_dim, write_to_file, matlab_dir);

  // ------------------------------------------------------------------------------------
  // - i2iSet
  // ------------------------------------------------------------------------------------
  assembler.assemble(SemiLinearization(eq,u,u,du),options.assemblyThreads);
  AssembledGalerkinOperator<Assembler> Ass(assembler); 
  // ------------------------------------------------------------------------------------
  // construct mass and stiffness matrix from semi-implicit structure
  // ------------------------------------------------------------------------------------
  Matrix A_;
  Matrix M_;
  Matrix K_;

  A_ = assembler.template get<Matrix>(false);
  assembler.assemble(SemiLinearization(eq,u,u,du),Assembler::RHS,options.assemblyThreads);
  auto rhs = assembler.rhs();
  // writeToMatlab(assembler,matlab_dir+"/matrixA_", "A");  

  F.Mass_stiff(1);
  SemiImplicitEulerStep<Functional>  eqM(&F,options.dt);
  eqM.setTau(0);
  assembler.assemble(SemiLinearization(eqM,u,u,du), Assembler::MATRIX, options.assemblyThreads);  
  M_ = assembler.template get<Matrix>(false);
  // writeToMatlab(assembler,matlab_dir+"/matrixM_", "M"); 

  // get stiffness 
  F.Mass_stiff(0);
  SemiImplicitEulerStep<Functional>  eqK(&F,options.dt);
  eqK.setTau(1);
  assembler.assemble(SemiLinearization(eqK,u,u,du), Assembler::MATRIX, options.assemblyThreads); 
  K_ = assembler.template get<Matrix>(false);
  K_*=(-options.dt);
  // writeToMatlab(assembler,matlab_dir+"/matrixK_", "K"); 
  
  // ------------------------------------------------------------------------------------ 
  // compute rhs based on petsc structure
  // ------------------------------------------------------------------------------------ 
  Vector rhs_vec_test(nDofs);
  rhs.write(rhs_vec_test.begin());
  Vector rhs_petsc_test(nDofs);
  rhs.write(rhs_petsc_test.begin());
  petsc_structure_rhs(sequenceOfTags, startingIndexOfTag, map_II, map_GammaGamma, rhs_vec_original,rhs_petsc_test);
  // ------------------------------------------------------------------------------------ 
  // compute rhs of each based on petsc structure
  // ------------------------------------------------------------------------------------ 
  std::vector<Vector> Fs_petcs;
  petsc_structure_rhs_petsc(sequenceOfTags, startingIndexOfTag, map_II, map_GammaGamma, map_GammaNbr, rhs_vec_original, map_indices, sharedDofsKaskade, Fs_petcs);
  std::vector<Matrix> subMatrices;
  std::vector<Matrix> subMatrices_M;
  std::vector<Matrix> subMatrices_K;
  // ------------------------------------------------------------------------------------ 
  // compute rhs based on petsc structure
  // ------------------------------------------------------------------------------------
  std::vector<Vector> weights; 
  construct_submatrices_petsc(map_nT2oT,
                              gridManager,
                              F,
                              variableSetDesc, 
                              spaces,
                              gridManager.grid(), 
                              u,
                              dt,
                              sequenceOfTags, 
                              startingIndexOfTag,
                              map_II,
                              map_GammaGamma,
                              map_GammaNbr,
                              sequanceOfsubdomains,
                              map_indices, 
                              A_,K_,M_,
                              rhs_petsc_test,
                              nDofs,
                              options.assemblyThreads,
                              write_to_file,
                              matlab_dir,
                              Fs_petcs,
                              weights,
                              subMatrices,
                              subMatrices_M,
                              subMatrices_K);


  if(write_to_file) generate_Interror_and_Interfaces_indices(sequenceOfTags, map_II, map_GammaGamma, map_GammaGamma_W_Nbr, map_indices, matlab_dir);

  return 0;
}
