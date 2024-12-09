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
#include <thread>

// A sample function representing work on a subdomain
void processSubdomain(int subdomainId) {
    std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simulate work
    std::cout << "Processed subdomain: " << subdomainId << " on thread: "
              << std::this_thread::get_id() << "\n";
}

template<class i3tuple>
void processSubdomain(int subIdx, const std::vector<int>& sequenceOfTags, 
                      const std::map<int, std::vector<int>>& sequanceOfsubdomainsKaskade,std::map<int, std::vector<i3tuple>> & MapSharedDofsKaskadeTuple, std::mutex& mapMutex) {
    int tag = sequenceOfTags[subIdx];  // Get the tag for the current subdomain
    const std::vector<int>& tmp = sequanceOfsubdomainsKaskade.at(tag);  // Get the list of subdomains for the tag

    for (int i = 0; i < tmp.size(); ++i) {
        auto it = MapSharedDofsKaskadeTuple.find(tmp[i]);  // Find the subdomain in the map
        
        std::lock_guard<std::mutex> lock(mapMutex); // Lock the map for thread safety

        if (it != MapSharedDofsKaskadeTuple.end()) {
            // If the subdomain already exists, add the new tuple
            std::vector<i3tuple>& values = it->second;
            values.push_back(i3tuple(subIdx, i, tmp[i]));
        } else {
            // Otherwise, create a new entry for the subdomain
            MapSharedDofsKaskadeTuple[tmp[i]] = {i3tuple(subIdx, i, tmp[i])};
        }
    }
}


int main(int argc, char* argv[])
{
  using namespace Kaskade::BDDC;
  using namespace boost::fusion;

  std::cout << "Start subdomain tutorial program" << std::endl;
  
  constexpr int dim = SPACEDIM;
  int refinements, order, solver, refinements_sol, interfaceTypes, iter_cg_with_bddc;
  double penalty, sigma_i, sigma_e, C_m, R, R_extra, tol, dt;
  bool  direct, onlyLowerTriangle, vtk_, timing, test_newCof;
  bool run_implicit_CG, run_implicit_CG_SDC, run_implicit_CG_BDDC, run_implicit_CG_SDC_BDDC, run_implicit_CG_BDDC_Fused, run_implicit_CG_SDC_BDDC_first_Sweep;
  bool run_implicit_CG_SDC_BDDC_all_collocation_once;
  bool run_implicit_CG_SDC_BDDC_smallest_collocation;
  bool run_implicit_CG_SDC_BDDC_all_collocation_once_update;
  std::string inputfile, early_excited, extra_set, intra_set, dir_out, matlab_dir;
  bool cg_semi, plot, withSplitFace,cg_solver;
  bool test_mesh_data, write_to_file;
  bool BDDC_SDC_with_initial, BDDC_verbose, BDDC_SDC_verbose;
  int verbose, assemblyThreads;
  CardiacIntegrationOptions options;
  if (getKaskadeOptions(argc,argv,Options
  // ("input",                    inputfile,                           "./input/coarse_4elem.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/coarse_4elem_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/coarse_4elem_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/coarse_4elem_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/coarse_emi.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/coarse_emi_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/coarse_emi_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/coarse_emi_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/coarse_kaskade_e.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/coarse_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/coarse_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/coarse_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/example4subc_3extra_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/example4subc_3extra_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/example4subc_3extra_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/example4subc_3extra_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/example3subc_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/example3subc_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/example3subc_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/example3subc_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/example4subc_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/example4subc_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/example4subc_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/example4subc_early_excited.txt","subdomain definition")
  ("input",                    inputfile,                           "./input/example4subc_2extra_mesh.vtu","subdomain definition")
  ("extra_set",                extra_set,                           "./input/example4subc_2extra_list_extracellular.txt","subdomain definition")
  ("intra_set",                intra_set,                           "./input/example4subc_2extra_list_intracellular.txt","subdomain definition")
  ("excited",                  early_excited,                       "./input/example4subc_2extra_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/example4subc_join_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/example4subc_join_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/example4subc_join_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/example4subc_join_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/example4subc_2extra_mesh_old.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/example4subc_2extra_list_extracellular_old.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/example4subc_2extra_list_intracellular_old.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/example4subc_2extra_early_excited_old.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/twoCells3d_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/twoCells3d_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/twoCells3d_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/twoCells3d_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/twoCells3d_mesh_new.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/twoCells3d_list_extracellular_new.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/twoCells3d_list_intracellular_new.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/twoCells3d_early_excited_new.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/tenCells3d_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/tenCells3d_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/tenCells3d_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/tenCells3d_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/2Cells3d_2extra_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/2Cells3d_2extra_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/2Cells3d_2extra_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/2Cells3d_2extra_early_excited.txt","subdomain definition") 
  // ("input",                    inputfile,                           "./input/cube3D.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/cube3D_extra.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/cube3D_intra.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/cube3D_early_eacited.txt","subdomain definition") 
  // ("input",                    inputfile,                           "./input/emiGrid.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/emiGrid_extraCellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/emiGrid_intraCellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid3nx1ny1nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular3nx1ny1nz.txt","subdomain definition")
    // ("intra_set",                intra_set,                           "./input/intracellular3nx1ny1nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid3nx3ny1nz.vtu","subdomain definition")
    // ("intra_set",                intra_set,                           "./input/extracellular3nx3ny1nz.txt","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/intracellular3nx3ny1nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid3nx3ny2nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular3nx3ny2nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular3nx3ny2nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid3nx3ny3nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular3nx3ny3nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular3nx3ny3nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid4nx4ny4nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular4nx4ny4nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular4nx4ny4nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid5nx5ny5nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular5nx5ny5nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular5nx5ny5nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid6nx6ny6nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular6nx6ny6nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular6nx6ny6nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid7nx7ny7nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular7nx7ny7nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular7nx7ny7nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
    // ("input",                    inputfile,                           "./input/emiGrid8nx8ny8nz.vtu","subdomain definition")
    // ("extra_set",                extra_set,                           "./input/extracellular8nx8ny8nz.txt","subdomain definition")
    // ("intra_set",                extra_set,                           "./input/intracellular8nx8ny8nz.txt","subdomain definition")
    // ("excited",                  early_excited,                       "./input/emiGrid_earlyExited.txt","subdomain definition") 
  // ("input",                    inputfile,                           "./input/emiGridnx5ny1nz1.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/emiGridnx5ny1nz1_extraCellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/emiGridnx5ny1nz1_intraCellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/emiGridnx5ny1nz1_earlyExited.txt","subdomain definition")       
  // ("input",                    inputfile,                           "./input/Materials.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_inner_new.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_gap_new.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_gap_membrane_new.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_gap_membrane_inner_new.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_membrane_only.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_inner_only.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_gap_only.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_gap_membrane.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_gap_membrane_inner.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_x.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/perturbed_all.vtu","subdomain definition")
  //("input",                    inputfile,                           "./input/10Cells3d_10extra_mesh_rescaled.vtu","subdomain definition")
  //("input",                    inputfile,                           "./input/10Cells3d_10extra_mesh_unconstructed.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/10Cells3d_10extra_mesh_refine2.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/10Cells3d_10extra_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/10Cells3d_10extra_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/10Cells3d_10extra_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/10Cells3d_10extra_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/20Cells3d_20extra_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/20Cells3d_20extra_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/20Cells3d_20extra_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/20Cells3d_20extra_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/40Cells3d_40extra_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/40Cells3d_40extra_list_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/40Cells3d_40extra_list_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/40Cells3d_40extra_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/robin_mesh.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/robin_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/robin_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/robin_early_excited.txt","subdomain definition")
  // ("input",                    inputfile,                           "./input/40cells3D.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/40cells3D_early_excitedtxt.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/40cells3D_list_extracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/40cells3D_list_intracellular.txt","subdomain definition")
 // ("input",                    inputfile,                           "./input/pepe_combi_domi.vtu","subdomain definition")
// //("input",                    inputfile,                           "./input/pepe_combi_domi_smaller.vtu","subdomain definition")
//   // ("input",                    inputfile,                           "./input/pepe_combi_domi_smaller_more.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/pepe_combi_domi_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/pepe_combi_domi_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/pepe_combi_domi_excited.txt","subdomain definition")
  //  ("input",                    inputfile,                           "./input/pepe_sep_domi.vtu","subdomain definition")
  // // ("input",                    inputfile,                           "./input/pepe_combi_domi_smaller.vtu","subdomain definition")
  // // ("input",                    inputfile,                           "./input/pepe_combi_domi_smaller_more.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/pepe_sep_domi_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/pepe_sep_domi_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/pepe_sep_domi_excited.txt","subdomain definition")
//  ("input",                    inputfile,                           "./input/robin_combi_domi.vtu","subdomain definition")
//  ("input",                    inputfile,                           "./input/robin_combi_domi_smaller.vtu","subdomain definition")
  // ("input",                    inputfile,                           "./input/robin_combi_domi_smaller_more.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/robin_combi_domi_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/robin_combi_domi_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/robin_combi_domi_excited.txt","subdomain definition")
  //   ("input",                    inputfile,                           "./input/robin_sep_domi.vtu","subdomain definition")
  // // ("input",                    inputfile,                           "./input/robin_sep_domi_smaller.vtu","subdomain definition")
  // // ("input",                    inputfile,                           "./input/robin_sep_domi_smaller_more.vtu","subdomain definition")
  // ("extra_set",                extra_set,                           "./input/robin_sep_domi_extracellular.txt","subdomain definition")
  // ("intra_set",                intra_set,                           "./input/robin_sep_domi_intracellular.txt","subdomain definition")
  // ("excited",                  early_excited,                       "./input/robin_sep_domi_excited.txt","subdomain definition")
  // ("dir",                      dir_out,                             "./output_rescaled0","subdomain definition")
  // ("matlab_dir",               matlab_dir,                          "./matlab_dir_rescaled0","subdomain definition")
  //   ("dir",                      dir_out,                             "./output_rescaled1","subdomain definition")
  // ("matlab_dir",               matlab_dir,                          "./matlab_dir_rescaled1","subdomain definition")
  //   ("dir",                      dir_out,                             "./output_rescaled2","subdomain definition")
  // ("matlab_dir",               matlab_dir,                          "./matlab_dir_rescaled2","subdomain definition")
  ("dir",                      dir_out,                             "./output","subdomain definition")
  ("matlab_dir",               matlab_dir,                          "./matlab_dir","subdomain definition")
  // ("dir",                      dir_out,                             "/scratch/htc/fchegini/output","subdomain definition")
  // ("matlab_dir",               matlab_dir,                          "/scratch/htc/fchegini/matlab_dir","subdomain definition")
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
  ("iter_cg_bddc",             iter_cg_with_bddc,                   100,"number of BDDC iterations")
  ("timing",                   timing,                              true,"whether to write timing info")
  ("test",                     test_mesh_data,                      false,"debug mode")
  ("run_implicit_CG",          run_implicit_CG,                     true, "run linearly semi-implicit method + CG")
  ("run_implicit_CG_SDC",      run_implicit_CG_SDC,                 false, "run linearly semi-implicit method + CG + Jacobi + SDC")
  ("run_implicit_CG_BDDC",     run_implicit_CG_BDDC,                true, "run linearly semi-implicit method + CG + Jacobi + SDC ")
  ("run_implicit_CG_BDDC_Fused",run_implicit_CG_BDDC_Fused,         false, "run linearly semi-implicit method + CG + Jacobi + SDC ")
  ("run_implicit_CG_SDC_BDDC", run_implicit_CG_SDC_BDDC,            false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("run_implicit_CG_SDC_BDDC_all_collocation_once", run_implicit_CG_SDC_BDDC_all_collocation_once_update,            false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("run_implicit_CG_SDC_BDDC_all_collocation_once", run_implicit_CG_SDC_BDDC_all_collocation_once,            false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("run_implicit_CG_SDC_BDDC_smallest_collocation", run_implicit_CG_SDC_BDDC_smallest_collocation,            false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("run_implicit_CG_SDC_BDDC_first_Sweep", run_implicit_CG_SDC_BDDC_first_Sweep,false, "run linearly semi-implicit method + CG + BDDC + SDC " )
  ("test_newCof",              test_newCof,                         false,"to test the coefficients")
  ("withSplitFace",            withSplitFace,                       false,"split faces in BDDC")  
  ("cg_solver",                cg_solver,                           true,"split faces in BDDC")  
  ("write_to_file",            write_to_file,                       true,"write to matlab file")  
  ("maxSteps",                 options.maxSteps,                    5,  "max number of time steps")
  ("vtk",                      options.writeVTK,                    1,  "write VTK output files 0=none, 1=time steps 2=sweeps")
  ("T_",                       options.T,                           0.01,  "final time[ms]")
  ("dt",                       options.dt,                          0.01,  "time step size[ms]")
  ("orderU",                   options.order,                       1,  "FE ansatz order for transmembrane voltage & action potential")
  ("atol",                     options.aTol,                        1e-15,  "absolute L^2 tolerance")
  ("stol",                     options.tolSelect,                   1e-5,  "L^inf tol for DoF selection")
  ("maxCGIter",                options.maxCGIter,                   10000,  "maximum number of IterateType::CG iterations in linear solver (0=direct solver)")
  ("cgTol",                    options.cgTol,                       1e-8,  "absolute IterateType::CG energy error tolerance")
  ("adapt",                    options.adapt,                       false,  "do adaptivity or not")
  ("sweeps",                   options.minSweeps,                   5,  "minimal number of SDC sweeps")
  ("maxSweeps",                options.maxSweeps,                   5,  "maximal number of SDC sweeps")
  ("nColloc",                  options.nCollocU,                    3,  "number of collocation points in time")
  ("nCollocStart",             options.nCollocUstart,               1,  "start sweeps with that many collocation points")
  ("verbose",                  options.verbosity,                   1,  "output density")
  ("sweepType",                options.sweepType,                   1,  "0: Euler, 1: LU")
  ("nReactionSweeps",          options.nReactionSweeps,             0,  "number of post-sweep Euler steps for reaction nonlinearity")
  ("minRefLevel",              options.minRefLevel,                 0,  "keep this level on refinement") 
  ("nThreads",                 options.assemblyThreads,             28,  "# of threads in assembler") 
  ("rosenbrockRefinement",     options.rosenbrockRefinementStyle,   false,  "do spatial refinement first, before SDC") 
  ("plot",                     options.plot,                        false,  "final time[ms]")
  ("CG_shift",                 options.CG_shift,                    true,  "shift the cg update")
  ("SDC_TOL",                  options.SDC_TOL,                     1e-6,  "SDC tolerance")
  ("sdc_contraction",          options.sdc_contraction,             0.2,  "SDC constraction")  
  ("BDDC_SDC_with_initial",    BDDC_SDC_with_initial,               true,  "BDDC with initial guess from previous collocation sol") 
  ("BDDC_verbose",             BDDC_verbose,                        true,  "SDC tolerance") 
  ("BDDC_SDC_verbose",         BDDC_SDC_verbose,                    true,  "SDC tolerance") 
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
  std::cout << "mesh: " << inputfile <<std::endl;
  std::cout << "dt: " << options.dt <<std::endl;
  std::cout << "order: " << order <<std::endl;
  std::cout << "refinement: "<< refinements << std::endl;
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;

  auto& timer = Timings::instance();

  // ---------------------------------------------------------------------------------------
  // REMOVE!!!
  // ---------------------------------------------------------------------------------------
  // unsigned int numThreads = std::thread::hardware_concurrency();
  unsigned int numThreads = std::max(1u, std::thread::hardware_concurrency());

  if (numThreads == 0) {
      std::cout << "Unable to detect the number of hardware threads.\n";
  } else {
      std::cout << "Number of hardware threads available: " << numThreads << "\n";
  }

  // int n = 10; // Total number of subdomains
  // unsigned int m = numThreads;  // Number of threads

  // {
  //   std::mutex mapMutex; // Mutex to protect shared resource MapSharedDofsKaskadeTuple
  //   typedef std::tuple<int, int, int> i3tuple; // Defining the i3tuple type

  //   // Map to store the shared dofs kaskade tuples
  //   std::map<int, std::vector<i3tuple>> MapSharedDofsKaskadeTuple;

  //   std::vector<int> sequenceOfTags = {0, 1, 2, 3, 4};  // Example sequence of tags
  //   std::map<int, std::vector<int>> sequanceOfsubdomainsKaskade = {
  //       {0, {0, 1, 2}}, {1, {3, 4, 5}}, {2, {6, 7, 8}},
  //       {3, {9, 10, 11}}, {4, {12, 13, 14}}
  //   }; // Example map of subdomains for each tag

  //   int n = sequenceOfTags.size();
  //   int m = 4; // Number of threads
  //   std::vector<std::thread> threads;

  //   // Launch threads to process subdomains concurrently
  //   for (int subIdx = 0; subIdx < n; ++subIdx) {
  //       threads.emplace_back([subIdx, &sequenceOfTags, &sequanceOfsubdomainsKaskade, &MapSharedDofsKaskadeTuple, &mapMutex]() {
  //           processSubdomain(subIdx, sequenceOfTags, sequanceOfsubdomainsKaskade, MapSharedDofsKaskadeTuple, mapMutex);
  //       });
  //   }

  //   // Wait for all threads to finish
  //   for (auto& t : threads) {
  //       if (t.joinable()) {
  //           t.join();
  //       }
  //   }

  //   // Optionally, print the contents of MapSharedDofsKaskadeTuple for debugging
  //   for (const auto& entry : MapSharedDofsKaskadeTuple) {
  //       std::cout << "Key: " << entry.first << "\n";
  //       for (const auto& val : entry.second) {
  //           std::cout << "  (" << std::get<0>(val) << ", " 
  //                     << std::get<1>(val) << ", " 
  //                     << std::get<2>(val) << ")\n";
  //       }
  //   }
    
  // }
  // // Sequential execution
  // auto sequentialStart = std::chrono::high_resolution_clock::now();

  // for (int i = 0; i < n; ++i) {
  //     processSubdomain(i); // Process each subdomain one by one
  // }

  // auto sequentialEnd = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> sequentialDuration = sequentialEnd - sequentialStart;
  // std::cout << "Sequential execution took: " << sequentialDuration.count() << " seconds.\n";

  // // Multithreaded execution
  // auto parallelStart = std::chrono::high_resolution_clock::now();

  // std::vector<std::thread> threads; // To manage threads
  // int chunkSize = (n + m - 1) / m; // Calculate chunk size (ceil(n/m))

  // for (int i = 0; i < m; ++i) {
  //     int start = i * chunkSize;
  //     int end = std::min(start + chunkSize, n);

  //     // Launch a thread to process the range of subdomains
  //     threads.emplace_back([start, end]() {
  //         for (int j = start; j < end; ++j) {
  //             processSubdomain(j);
  //         }
  //     });
  // }

  // // Wait for all threads to finish
  // for (auto& t : threads) {
  //     if (t.joinable()) {
  //         t.join();
  //     }
  // }

  // auto parallelEnd = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> parallelDuration = parallelEnd - parallelStart;
  // std::cout << "Parallel execution took: " << parallelDuration.count() << " seconds.\n";

  // // Comparing results
  // double speedup = sequentialDuration.count() / parallelDuration.count();
  // std::cout << "Speedup: " << speedup << "x\n";

  // return 0;
  // ------------------------------------------------------------------------------------------------------------
  // get the material from the mesh  
  // ------------------------------------------------------------------------------------------------------------
  timer.start("reading FE grid");
  VTKReader vtk(inputfile);
  GridManager<Grid> gridManager(vtk.createGrid<Grid>());
  gridManager.enforceConcurrentReads(true);
  using H1SpaceMaterial = FEFunctionSpace<DiscontinuousLagrangeMapper<double,LeafView>>;
  using cellMaterial = H1SpaceMaterial::Element<1>::type;
  H1SpaceMaterial materialSpace(gridManager,gridManager.grid().leafGridView(), 0);
  cellMaterial material(materialSpace);
  vtk.getCoefficients("domain",material);
  std::cout << "---------- refinements ---------- "<< std::endl;

  std::cout << "sigma_i = " << sigma_i << " sigma_e = "<< sigma_e << std::endl;
  gridManager.globalRefine(refinements);
  if(options.plot) writeVTK(material,out+"/Materials",IoOptions(),"domain");
  writeVTK(material,out+"/Materials",IoOptions(),"domain");
  timer.stop("reading FE grid");
  // return 0;
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
  std::cout<< "n_extra_set: " << n_extra_set <<std::endl;
  std::vector<int> arr_extra(n_extra_set);
  getSubdomain(arr_extra, file_extra_list);

  std::ifstream file_intra_list(intra_set);
  int n_intra_set;
  file_intra_list >> n_intra_set;
  std::cout<< "n_intra_set: " << n_intra_set <<std::endl;
  std::vector<int> arr_intra_set(n_intra_set);
  getSubdomain(arr_intra_set, file_intra_list);

  std::ifstream file_excited_region(early_excited);
  int n_excited_region;
  file_excited_region >> n_excited_region;
  std::cout<< "n_excited_region: " << n_excited_region <<std::endl;
  std::vector<int> arr_excited_region(n_excited_region);
  getSubdomain(arr_excited_region, file_excited_region);

  Dune::FieldVector<double,dim> zero(0);


  // assume that only extra cellular has the tag number zero
  // in case of having more than one extracellular subdomain, 
  // for the aszast space, we consider only one extra cellular subdoamjn, thats' why we 
  // map all the extracellular tags to zero,
  // however in order to assemble matrices for BDDC, we need to consider the origonal tags numbers.  
  std::map<int, Dune::FieldVector<double,1>> map_OriginalTag_anzastTag;
  for (int i = 0; i < arr_extra.size(); ++i)
  {
    Dune::FieldVector<double,1> extra_zero(0);
    extra_zero[0] = 0.0;
    map_OriginalTag_anzastTag.insert({ arr_extra[i], extra_zero }); 
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
                                                                            if(arr_extra.size()>1){
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
  F.extracellular_materials(arr_extra);
  typedef typename Functional::OriginVars::VariableSet State;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  constexpr int neq = Functional::TestVars::noOfVariables;

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
  std::cout << "The numeber of dofs: "<< nDofs <<std::endl;
  std::cout << "The numeber of cells :: "<<gridManager.grid().size(0)<<std::endl;
  size_t dof_size = variableSetDesc.degreesOfFreedom(0, 1);
  // ------------------------------------------------------------------------------------
  // initila the data
  // ------------------------------------------------------------------------------------
  F.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
  uAll = component<0>(u);
  //if(options.plot) 
  writeVTK(uAll,out+"/initialTest",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
  dt = options.dt;
  SemiImplicitEulerStep<Functional>  eq(&F,dt);
  double const tau = dt;
  eq.setTau(tau);
  auto du(u);
  du *= 0;
  std::cout << "Start: assembler.assemble(SemiLinearization "<<std::endl;
  assembler.assemble(SemiLinearization(eq,u,u,du),Assembler::RHS,options.assemblyThreads);
  std::cout << "END: assembler.assemble(SemiLinearization "<<std::endl;
  auto rhs_oiginal = assembler.rhs();




  // ------------------------------------------------------------------------------------
  // semi implicit + CG methods
  // ------------------------------------------------------------------------------------
  // {
  //   if(run_implicit_CG){

  //     std::vector<int> sequenceOfTags(10);//n_subdomains);
  //     std::map<int, int> map_indices;
  //     std::map<int,std::set<int>> map_II;                            // II
  //     std::map<int,std::set<int>> map_GammaGamma_noDuplicate;        // GammaGamma_nodup

  //     Vector sol_semi(nDofs);
  //     Functional F_semi( material,
  //                 gridManager.grid(),
  //                 spaces,
  //                 penalty,
  //                 sigma_i,
  //                 sigma_e,
  //                 C_m,  
  //                 R,
  //                 R_extra);
  //     F_semi.extracellular_materials(arr_extra);
  //     F_semi.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
  //     uAll = component<0>(u);
  //     if(options.plot) writeVTK(uAll,out+"/initialSemiF",
  //              IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
  //     timer.start("linearly semi implicit method");
  //     std::cout << "---------------------------------------------" << std::endl;
  //     std::cout << "semi implict approach" << std::endl;
  //     std::cout << "---------------------------------------------" << std::endl;


  //     uAll = component<0>(u);
  //     u = semiImplicit_CG_Jacobi( gridManager,
  //                                 F_semi,
  //                                 variableSetDesc,
  //                                 spaces,
  //                                 gridManager.grid(),
  //                                 options,
  //                                 out,
  //                                 cg_semi, 
  //                                 direct,
  //                                 u,
  //                                 uAll,
  //                                 sol_semi,
  //                                 matlab_dir,
  //                                 sequenceOfTags, 
  //                                 map_indices, 
  //                                 map_II, 
  //                                 map_GammaGamma_noDuplicate);  
  //     timer.stop("linearly semi implicit method");

  //     // {
  //     //   Vector sol_semi_to_petsc(sol_semi);
  //     //   sol_semi_to_petsc = 0;
  //     //    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_semi,sol_semi_to_petsc);
  //     //   if(write_to_file) writeSolution(sol_semi_to_petsc,matlab_dir+"/sol");
  //     // }
  //   }
  // }
  // return 0;

  // ------------------------------------------------------------------------------------
  // Extract the mesh data
  // - II, GammaGamma, IGamma, GammaGamma_W_Nbr, gamma_nbrs, sequenceOfsubdomains 
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
  std::vector<std::set<int>> i2t(dof_size);                      //index to tags
  std::vector<std::set<int>> i2i(dof_size);                      //index to index
  std::set<int> interface_extra_dofs;                            // set of dofs on the extracellular interfaces 

  std::map<std::pair<int, int>, std::vector<double>> coord;      //coordinates of each dofs: key index and tag
  std::map<int, std::vector<double>> coord_globalIndex;         //coordinates of each dofs: key index
  std::vector<int> i2Tag(dof_size);                              //index to tags
  std::set<int> tags;                                            // list of tags
  std::map<int, int> map_t2l;                                    //map: tag to lenth
  std::map<int, int> map_sT2l;                                   //map: sequance of each tag to length
  std::map<int,int> map_nT2oT;                                   //map: new Tag to original Tag
                              
  std::map<int,std::set<int>> map_II;                            // II
  std::map<int,std::set<int>> map_IGamma;                        // IGamma
  std::map<int,std::set<int>> map_IGamma_noDuplicate;            // IGamma
  std::map<int,std::set<int>> map_GammaGamma;                    // GammaGamma
  std::map<int,std::set<int>> map_GammaGamma_noDuplicate;        // GammaGamma_nodup
  std::map<int,std::set<int>> map_GammaGamma_W_Nbr;              // GammaGamma with nbr
  std::map<int,std::set<int>> map_GammaNbr_Nbr;                  // GammaGamma only nbr without out the extra neighors...
  std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate;      // GammaGamma only nbr without out the extra neighors...
  std::map<int,std::map<int,std::set<int>>> map_GammaNbr;        // GammaNbr
  std::map<int,std::vector<int>> sequenceOfsubdomains;           // sequence of neighboring tags for each subdomain

  std::map<int,bool> map_markCorners;                              // if the indices on the corners when we have more than one extracellular
  std::set<int> corners;
  mesh_data_structure(boost::fusion::at_c<0>(u.data),  
                      material, arr_extra,
                      e2i, i2e, i2t, e2e, i2i, coord, coord_globalIndex, i2Tag, tags,
                      map_t2l, map_sT2l, map_II, map_IGamma, map_GammaGamma, map_IGamma_noDuplicate, 
                      map_GammaGamma_noDuplicate, map_GammaGamma_W_Nbr, map_GammaNbr_Nbr, 
                      map_GammaNbr_Nbr_noDuplicate, map_GammaNbr, interface_extra_dofs,sequenceOfsubdomains);
  int n_subdomains = map_II.size();

  std::vector<int> sequenceOfTags(n_subdomains);
  std::vector<int> sequenceOfTags_extra(n_extra_set); // only extra cellular
  // extracellular is only even number
  std::map<int,int> startingIndexOfTag;
  std::map<int,int> Tag2IndexSub;
  computed_sequenceOfTags(map_t2l,map_IGamma_noDuplicate, map_GammaNbr, sequenceOfTags, sequenceOfTags_extra, startingIndexOfTag, Tag2IndexSub, map_nT2oT);

  if(false)
  {
  std::cout <<"==========================\n";
  std::cout <<"sequenceOfsubdomains\n";
  std::cout <<"==========================\n";
  for (int i = 0; i < sequenceOfTags.size(); ++i)
  {
    int tag = sequenceOfTags[i];
    std::vector<int> vec_nbr = sequenceOfsubdomains[tag];

    std::cout << tag<<": ";
    for (int nbr = 0; nbr < vec_nbr.size(); ++nbr)
    {
      std::cout << vec_nbr[nbr]<<" ";
    }
    std::cout <<"\n";
  }
  std::cout <<"==========================\n";
  }


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
  subdomain_indices(sequenceOfTags, map_II, map_GammaGamma_noDuplicate, map_GammaNbr_Nbr_noDuplicate, local2Global, global2Local, globalIndices);

  std::map<int, int> map_indices;
  std::map<int, int> map_index_to_subdomain;
  map_kaskade2petcs(sequenceOfTags, map_II, map_GammaGamma_noDuplicate, map_indices, map_index_to_subdomain);
  std::cout <<"DONE!!!" <<std::endl;
  // return 0;

  std::set<std::set<int>> i2iSet_;
  removeInnerIndices_i2i(i2i);
  std::set<std::set<int>> i2iSet(i2i.begin(),i2i.end());  //index to index only those has more than one neighours on the interfaces ?

  {
    double precision = 16;
    std::string fname = matlab_dir+"/corners.m";
    //std::ofstream f(fname.c_str());

    std::set<std::set<int>>::iterator it;
    for (it = i2iSet.begin(); it != i2iSet.end(); ++it) 
    {
      std::set<int> s = *it;
      std::set<int>::iterator itr;

      for (itr = s.begin(); itr != s.end(); ++itr) 
      {
        if(s.size()>2) corners.insert(map_indices[*itr]);
        //if(s.size()>2) f << map_indices[*itr] << " "; 
        //if(false) std::cout << map_indices[*itr]<< " ";
      }
      //if(s.size()>2) f << "\n";
      //if(false) std::cout <<std::endl;
    }
  }

  marked_corners(arr_extra, sequenceOfTags, map_indices, i2t, i2iSet, map_GammaNbr_Nbr, matlab_dir, map_markCorners);

  std::cout << "write Dirichlet and coordinates!" << std::endl;
  std::set<int> dofsDirichlet;
  int mesh_dim = SPACEDIM==2? 2:3;
  write_Dirichlet_and_coordinates(boost::fusion::at_c<0>(u.data), material, arr_extra, e2i, map_indices, coord,coord_globalIndex, dof_size, mesh_dim, write_to_file, matlab_dir, dofsDirichlet);
  

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
  {
    assembler.assemble(SemiLinearization(eq,u,u,du),options.assemblyThreads);
    auto rhs = assembler.rhs();
    // writeToMatlab(assembler,matlab_dir+"/matrixA_RHS_", "A");  
  }

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
  Vector rhs_vec_original(nDofs);
  rhs_oiginal.write(rhs_vec_original.begin());
  Vector rhs_vec_test(nDofs);
  rhs.write(rhs_vec_test.begin());
  Vector rhs_petsc_test(nDofs);
  rhs.write(rhs_petsc_test.begin());
  petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, rhs_vec_original,rhs_petsc_test);
  // ------------------------------------------------------------------------------------ 
  // compute rhs of each based on petsc structure
  // ------------------------------------------------------------------------------------ 
  std::vector<std::vector<LocalDof>> sharedDofsKaskade;
  compute_sharedDofsKaskade_moreExtraCells(sequenceOfTags, map_indices, map_II, map_GammaGamma, map_GammaGamma_W_Nbr, map_GammaNbr_Nbr, write_to_file, matlab_dir, sharedDofsKaskade);

  // ------------------------------------------------------------------------------------ 
  // compute rhs based on petsc structure
  // ------------------------------------------------------------------------------------
  std::cout << "generated sub matrices of EMI model for BDDC in petsc!" << std::endl;
  std::vector<Vector> Fs_petcs(n_subdomains);
  petsc_structure_rhs_subdomain_petsc(sequenceOfTags, map_II, map_GammaGamma_noDuplicate, rhs_vec_original, map_indices, sharedDofsKaskade, Fs_petcs);

  std::vector<int> cells(gridManager.grid().size(0)); // vector with size ints.
  std::iota (std::begin(cells), std::end(cells), 0);

  std::set<int> cells_set(cells.begin(),cells.begin());
  CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), cells_set, tags,material); 

  std::vector<Matrix> subMatrices(n_subdomains);
  std::vector<Matrix> subMatrices_M(n_subdomains);
  std::vector<Matrix> subMatrices_K(n_subdomains);
  std::vector<Vector> weights(n_subdomains);
  // construct_submatrices_petsc
  construct_submatrices_petsc_parallel(arr_extra,
                              map_nT2oT,
                              gridManager,
                              F,
                              Cellfltr,
                              variableSetDesc, 
                              spaces,
                              gridManager.grid(), 
                              u,
                              dt,
                              sequenceOfTags, 
                              startingIndexOfTag,
                              map_II,
                              map_GammaGamma,
                              map_GammaGamma_noDuplicate,
                              map_GammaNbr,
                              map_GammaNbr_Nbr_noDuplicate,
                              sequenceOfsubdomains,
                              map_indices,
                              map_markCorners,
                              cells_set,
                              tags,
                              i2Tag,
                              i2t,
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
  Cellfltr.select_based_on_tag(false);
  if(write_to_file) generate_Interror_and_Interfaces_indices(sequenceOfTags, map_II, map_GammaGamma_noDuplicate, map_GammaGamma_W_Nbr, map_indices, matlab_dir);

  // ------------------------------------------------------------------------------------ 
  // compute submatrices and rhs based on Kaskade structure
  // ------------------------------------------------------------------------------------
  std::cout << "generated sub matrices of EMI model for BDDC in kaskade!" << std::endl;
  std::vector<Matrix> As(n_subdomains);
  std::vector<Matrix> Ms(n_subdomains);
  std::vector<Matrix> Ks(n_subdomains);

  std::vector<Vector> Fs(n_subdomains);
  std::map<int,std::vector<int>> IG_seq;
  std::vector<std::vector<LocalDof>> sharedDofsKaskade_new;
  std::map<int,int> T2Index;

  // construct_As
  // Sequential execution
  auto sequentialStart_As = std::chrono::high_resolution_clock::now();

  std::vector<int> dofsDirichlet_vec(dofsDirichlet.begin(), dofsDirichlet.end());
  std::cout << "=========================================================="<<std::endl;
  //construct_As_parallel
  construct_As_parallel( arr_extra, 
                sequenceOfTags, 
                map_II, 
                map_GammaGamma, 
                map_GammaGamma_noDuplicate, 
                map_GammaNbr_Nbr_noDuplicate, 
                rhs_petsc_test, 
                weights, 
                map_indices, 
                matlab_dir,
                write_to_file,
                subMatrices,
                subMatrices_M,
                subMatrices_K, 
                options.dt,
                As, 
                Ms,
                Ks,
                Fs, 
                IG_seq, 
                dofsDirichlet_vec,
                sharedDofsKaskade_new,
                T2Index);

  auto sequentialEnd_AS = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> sequentialDuration_AS = sequentialEnd_AS - sequentialStart_As;
  std::cout << "Sequential execution took: " << sequentialDuration_AS.count() << " seconds.\n";

  //if(write_to_file)
  {
    
    int max_subdomain = 0;
    int min_subdomain = 1e10;
    for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
    {

      //std::cout << "As[subIdx].N()-> "<< As[subIdx].N() << " Ms[subIdx].N()-> " << Ms[subIdx].N()  << " Ks[subIdx].N()-> "  << Ks[subIdx].N()  << std::endl;
      //std::cout << subIdx << " As[subIdx].N()-> "<< As[subIdx].N() << std::endl;
      if(As[subIdx].N()>max_subdomain)
        max_subdomain = As[subIdx].N();

      if(As[subIdx].N()<min_subdomain)
        min_subdomain = As[subIdx].N();

      int tag = sequenceOfTags[subIdx];
      std::string path = std::to_string(subIdx);
      // writeToMatlabPath(As[subIdx],Fs[subIdx],"A_kaskade_shrinked"+path,matlab_dir, true);      
    }  
    std::cout << "A_.N() = " <<A_.N() << " max_subdomain "<< max_subdomain<< " max_subdomain "<< min_subdomain<< std::endl;
  }

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
      F_semi.extracellular_materials(arr_extra);
      F_semi.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);
      if(options.plot) writeVTK(uAll,out+"/initialSemiF",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");
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
                                  sol_semi,
                                  matlab_dir,
                                  sequenceOfTags, 
                                  map_indices, 
                                  map_II, 
                                  map_GammaGamma_noDuplicate);  
      timer.stop("linearly semi implicit method");

      {
        Vector sol_semi_to_petsc(sol_semi);
        sol_semi_to_petsc = 0;
         petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_semi,sol_semi_to_petsc);
        if(write_to_file) writeSolution(sol_semi_to_petsc,matlab_dir+"/sol");
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // semi implicit + CG + SDC methods
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG_SDC)
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implicit with SDC + CG                   " << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      Vector sol_SDC(nDofs);
      Functional F_SDC( material,
                  gridManager.grid(),
                  spaces,
                  penalty,
                  sigma_i,
                  sigma_e,
                  C_m,  
                  R,
                  R_extra);
      F_SDC.extracellular_materials(arr_extra);
      CardiacIntegrationStatistics statistics;
      F_SDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);

      if(options.plot) writeVTK(uAll,out+"/SDCInitial",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

      std::cout <<" test CellFilter!!!!\n";
      std::set<int> s_temp;
      for (int i = 0; i < gridManager.grid().size(0); ++i) 
        s_temp.insert(i);

      // CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), s_temp); 
      CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), cells_set, tags,material); 
      u = semiImplicit_CG_Jacobi_SDC( gridManager,
                                      F_SDC,
                                      Cellfltr,
                                      variableSetDesc,
                                      spaces,
                                      gridManager.grid(),
                                      sol_SDC,
                                      u,
                                      i2e,
                                      options,
                                      statistics,
                                      out,
                                      uAll,
                                      i2i
                                      );
      Vector sol_sdc_to_petsc(sol_SDC);
      sol_sdc_to_petsc = 0; 
      petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_SDC,sol_sdc_to_petsc);
      if(write_to_file) writeSolution(sol_sdc_to_petsc,matlab_dir+"/sol_sdc"); 
    }
  }

  // ------------------------------------------------------------------------------------
  // semi implicit + CG + BDDC methods
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG_BDDC)
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implicit with CG + BDDC                  " << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      Vector sol_BDDC(nDofs);
      Functional F_BDDC(  material,
                          gridManager.grid(),
                          spaces,
                          penalty,
                          sigma_i,
                          sigma_e,
                          C_m,  
                          R,
                          R_extra);
      F_BDDC.extracellular_materials(arr_extra);
      F_BDDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);

      u = semiImplicit_CG_BDDC( gridManager,
                                F_BDDC,
                                variableSetDesc,
                                spaces,
                                gridManager.grid(),
                                options,
                                out,
                                direct,
                                u,
                                uAll,
                                sol_BDDC,
                                sharedDofsKaskade_new,
                                interfaceTypes,
                                As,
                                sequenceOfTags, 
                                map_II,
                                map_GammaGamma_noDuplicate,
                                weights,
                                cg_solver,
                                iter_cg_with_bddc,
                                local2Global,
                                tol,
                                map_t2l,
                                map_indices,
                                BDDC_verbose,
                                IG_seq,
                                matlab_dir,
                                write_to_file
                                );  
        Vector sol_bddc_to_petsc(sol_BDDC);
        sol_bddc_to_petsc = 0; 
        petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_BDDC,sol_bddc_to_petsc);
        if(write_to_file) writeSolution(sol_bddc_to_petsc,matlab_dir+"/sol_bddc"); 
    }  
  }

  // // ------------------------------------------------------------------------------------
  // // semi implicit + CG + BDDC methods Fused
  // // ------------------------------------------------------------------------------------
  // {
  //   // if(run_implicit_CG_BDDC_Fused)
  //   // {
  //   //   std::map<int,std::set<int>> map_II_fused;                        // II
  //   //   std::map<int,std::set<int>> map_GammaGamma_fused;                // GammaGamma
  //   //   std::map<int,std::set<int>> map_GammaGamma_noDuplicate_fused;    // GammaGamma_nodup
  //   //   std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate_fused;  // GammaGamma only nbr without out the extra neighors...
  //   //   std::map<int,Matrix> subMatrices_fused;
  //   //   std::map<int, int> map_t2l_fused;                                    //map: tag to lenth
  //   //   merge_inner_interface_bddc_fused( sequenceOfTags_extra,
  //   //                                     T2Index,
  //   //                                     map_II,
  //   //                                     map_GammaGamma,
  //   //                                     map_GammaGamma_noDuplicate,
  //   //                                     map_GammaNbr_Nbr_noDuplicate,
  //   //                                     subMatrices,
  //   //                                     map_II_fused,
  //   //                                     map_GammaGamma_fused,
  //   //                                     map_GammaGamma_noDuplicate_fused,
  //   //                                     map_GammaNbr_Nbr_noDuplicate_fused,
  //   //                                     map_t2l_fused,
  //   //                                     subMatrices_fused);





  //   //   std::vector<Matrix> As_fused(n_subdomains);
  //   //   std::vector<Vector> Fs_fused(n_subdomains);
  //   //   std::vector<std::vector<int>> IG_seq_fused(n_subdomains);
  //   //   std::vector<std::vector<LocalDof>> sharedDofsKaskade_new_fused;
  //   //   construct_As_fused( arr_extra, 
  //   //                       sequenceOfTags_extra, 
  //   //                       map_II_fused, 
  //   //                       map_GammaGamma_fused, 
  //   //                       map_GammaGamma_noDuplicate_fused, 
  //   //                       map_GammaNbr_Nbr_noDuplicate_fused, 
  //   //                       rhs_petsc_test, 
  //   //                       weights, 
  //   //                       map_indices, 
  //   //                       matlab_dir,
  //   //                       write_to_file,
  //   //                       subMatrices_fused, 
  //   //                       As_fused, 
  //   //                       Fs_fused, 
  //   //                       IG_seq_fused, 
  //   //                       sharedDofsKaskade_new_fused);
  //   //   std::cout << "---------------------------------------------" << std::endl;
  //   //   std::cout << "semi implicit with CG + BDDC + Fused subdomain" << std::endl;
  //   //   std::cout << "---------------------------------------------" << std::endl;
  //   //   Vector sol_BDDC(nDofs);
  //   //   Functional F_BDDC(  material,
  //   //                       gridManager.grid(),
  //   //                       spaces,
  //   //                       penalty,
  //   //                       sigma_i,
  //   //                       sigma_e,
  //   //                       C_m,  
  //   //                       R,
  //   //                       R_extra);
  //   //   F_BDDC.extracellular_materials(arr_extra);
  //   //   F_BDDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
  //   //   uAll = component<0>(u);

  //   //   u = semiImplicit_CG_BDDC_fused( gridManager,
  //   //                             F_BDDC,
  //   //                             variableSetDesc,
  //   //                             spaces,
  //   //                             gridManager.grid(),
  //   //                             options,
  //   //                             out,
  //   //                             direct,
  //   //                             u,
  //   //                             uAll,
  //   //                             sol_BDDC,
  //   //                             sharedDofsKaskade_new_fused,
  //   //                             interfaceTypes,
  //   //                             As,               // pass the right one 
  //   //                             sequenceOfTags,   // sequenceOfTags for fused!
  //   //                             sequenceOfTags_extra, 
  //   //                             map_II_fused,                          // update 
  //   //                             map_GammaGamma_noDuplicate_fused,      // update
  //   //                             map_GammaNbr_Nbr_noDuplicate_fused,    // update 
  //   //                             weights,
  //   //                             cg_solver,
  //   //                             iter_cg_with_bddc,
  //   //                             local2Global,  // update
  //   //                             tol,
  //   //                             map_t2l_fused,       // update 
  //   //                             map_indices,
  //   //                             BDDC_verbose,
  //   //                             IG_seq_fused,        // update
  //   //                             matlab_dir,
  //   //                             write_to_file
  //   //                             );  
  //   //     Vector sol_bddc_to_petsc(sol_BDDC);
  //   //     sol_bddc_to_petsc = 0; 
  //   //     petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_BDDC,sol_bddc_to_petsc);
  //   //     if(write_to_file) writeSolution(sol_bddc_to_petsc,matlab_dir+"/sol_bddc"); 
  //   // }  
  // }

  // ------------------------------------------------------------------------------------
  // semi implicit + CG + SDC + BDDC methods
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG_SDC_BDDC)
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implicit with SDC + BDDC + CG           " << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      Vector sol_BDDC_SDC(nDofs);
      Functional F_BDDC_SDC(material,
                            gridManager.grid(),
                            spaces,
                            penalty,
                            sigma_i,
                            sigma_e,
                            C_m,  
                            R,
                            R_extra);
      F_BDDC_SDC.extracellular_materials(arr_extra);
      CardiacIntegrationStatistics statistics;
      F_BDDC_SDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);

      if(options.plot) writeVTK(uAll,out+"/emiSDCBDDCInitial",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

      std::cout <<" test CellFilter!!!!\n";
      std::set<int> s_temp;
      for (int i = 0; i < gridManager.grid().size(0); ++i) 
        s_temp.insert(i);

      CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), cells_set, tags,material); 
      u = semiImplicit_CG_BDDC_SDC( gridManager,
                                    F_BDDC_SDC,
                                    Cellfltr,
                                    variableSetDesc,
                                    spaces,
                                    gridManager.grid(),
                                    u,
                                    i2e,
                                    options,
                                    statistics,
                                    out,
                                    uAll,
                                    i2i,
                                    cg_semi,
                                    direct,
                                    matlab_dir,
                                    sol_BDDC_SDC,
                                    sharedDofsKaskade_new,
                                    interfaceTypes,
                                    n_subdomains,
                                    A_,
                                    M_,
                                    K_,
                                    As,
                                    Ms,
                                    Ks,
                                    IG_seq,
                                    map_IGamma,
                                    sequenceOfTags,
                                    map_II, 
                                    map_GammaGamma_noDuplicate, 
                                    weights,
                                    Fs_petcs,
                                    cg_solver,
                                    iter_cg_with_bddc,
                                    local2Global,
                                    global2Local,
                                    map_index_to_subdomain,
                                    tol,
                                    map_t2l, 
                                    map_indices,
                                    BDDC_SDC_with_initial,
                                    BDDC_verbose,
                                    BDDC_SDC_verbose);

    Vector sol_bddc_sdc_step(sol_BDDC_SDC);
    sol_bddc_sdc_step = 0; 
    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_BDDC_SDC,sol_bddc_sdc_step);
    if(write_to_file) writeSolution(sol_bddc_sdc_step,matlab_dir+"/sol_bddc_sdc"); 


    }  
  }

  // ------------------------------------------------------------------------------------
  // semi implicit + CG + SDC + BDDC methods all collocation_once
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG_SDC_BDDC_all_collocation_once_update)
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implicit with SDC + BDDC + CG  all collocations  update       " << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      Vector sol_BDDC_SDC_all_coll(nDofs);
      Functional F_BDDC_SDC(material,
                            gridManager.grid(),
                            spaces,
                            penalty,
                            sigma_i,
                            sigma_e,
                            C_m,  
                            R,
                            R_extra);
      F_BDDC_SDC.extracellular_materials(arr_extra);
      CardiacIntegrationStatistics statistics;
      F_BDDC_SDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);

      if(options.plot) writeVTK(uAll,out+"/emiSDCBDDCInitial",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

      std::cout <<" test CellFilter!!!!\n";
      std::set<int> s_temp;
      for (int i = 0; i < gridManager.grid().size(0); ++i) 
        s_temp.insert(i);

      CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), cells_set, tags,material); 
      u = semiImplicit_CG_BDDC_SDC_allCollocations_once_update( gridManager,
                                    F_BDDC_SDC,
                                    Cellfltr,
                                    variableSetDesc,
                                    spaces,
                                    gridManager.grid(),
                                    u,
                                    i2e,
                                    options,
                                    statistics,
                                    out,
                                    uAll,
                                    i2i,
                                    cg_semi,
                                    direct,
                                    matlab_dir,
                                    sol_BDDC_SDC_all_coll,
                                    sharedDofsKaskade_new,
                                    interfaceTypes,
                                    n_subdomains,
                                    A_,
                                    M_,
                                    K_,
                                    As,
                                    Ms,
                                    Ks,
                                    IG_seq,
                                    map_IGamma,
                                    sequenceOfTags,
                                    Tag2IndexSub,
                                    i2Tag,
                                    map_II, 
                                    map_GammaGamma_noDuplicate, 
                                    weights,
                                    Fs_petcs,
                                    cg_solver,
                                    iter_cg_with_bddc,
                                    local2Global,
                                    global2Local,
                                    map_index_to_subdomain,
                                    tol,
                                    map_t2l, 
                                    map_indices,
                                    BDDC_SDC_with_initial,
                                    BDDC_verbose,
                                    BDDC_SDC_verbose);

    Vector sol_bddc_sdc_step(sol_BDDC_SDC_all_coll);
    sol_bddc_sdc_step = 0; 
    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_BDDC_SDC_all_coll,sol_bddc_sdc_step);
    if(write_to_file) writeSolution(sol_bddc_sdc_step,matlab_dir+"/sol_bddc_sdc_small_all_coll"); 

    }  
  }

  // ------------------------------------------------------------------------------------
  // semi implicit + CG + SDC + BDDC methods all collocation_once
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG_SDC_BDDC_all_collocation_once)
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implicit with SDC + BDDC + CG  all collocations         " << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      Vector sol_BDDC_SDC_all_coll(nDofs);
      Functional F_BDDC_SDC(material,
                            gridManager.grid(),
                            spaces,
                            penalty,
                            sigma_i,
                            sigma_e,
                            C_m,  
                            R,
                            R_extra);
      F_BDDC_SDC.extracellular_materials(arr_extra);
      CardiacIntegrationStatistics statistics;
      F_BDDC_SDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);

      if(options.plot) writeVTK(uAll,out+"/emiSDCBDDCInitial",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

      std::cout <<" test CellFilter!!!!\n";
      std::set<int> s_temp;
      for (int i = 0; i < gridManager.grid().size(0); ++i) 
        s_temp.insert(i);

      CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), cells_set, tags,material); 
      u = semiImplicit_CG_BDDC_SDC_allCollocations_once( gridManager,
                                    F_BDDC_SDC,
                                    Cellfltr,
                                    variableSetDesc,
                                    spaces,
                                    gridManager.grid(),
                                    u,
                                    i2e,
                                    options,
                                    statistics,
                                    out,
                                    uAll,
                                    i2i,
                                    cg_semi,
                                    direct,
                                    matlab_dir,
                                    sol_BDDC_SDC_all_coll,
                                    sharedDofsKaskade_new,
                                    interfaceTypes,
                                    n_subdomains,
                                    A_,
                                    M_,
                                    K_,
                                    As,
                                    Ms,
                                    Ks,
                                    IG_seq,
                                    map_IGamma,
                                    sequenceOfTags,
                                    map_II, 
                                    map_GammaGamma_noDuplicate, 
                                    weights,
                                    Fs_petcs,
                                    cg_solver,
                                    iter_cg_with_bddc,
                                    local2Global,
                                    global2Local,
                                    map_index_to_subdomain,
                                    tol,
                                    map_t2l, 
                                    map_indices,
                                    BDDC_SDC_with_initial,
                                    BDDC_verbose,
                                    BDDC_SDC_verbose);

    Vector sol_bddc_sdc_step(sol_BDDC_SDC_all_coll);
    sol_bddc_sdc_step = 0; 
    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_BDDC_SDC_all_coll,sol_bddc_sdc_step);
    if(write_to_file) writeSolution(sol_bddc_sdc_step,matlab_dir+"/sol_bddc_sdc_small_all_coll"); 

    }  
  }

    // ------------------------------------------------------------------------------------
  // semi implicit + CG + SDC + BDDC methods all collocation_once
  // ------------------------------------------------------------------------------------
  {
    if(run_implicit_CG_SDC_BDDC_smallest_collocation)
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "semi implicit with SDC + BDDC + CG  smallest collocation         " << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
      Vector sol_BDDC_SDC_small_coll(nDofs);
      Functional F_BDDC_SDC(material,
                            gridManager.grid(),
                            spaces,
                            penalty,
                            sigma_i,
                            sigma_e,
                            C_m,  
                            R,
                            R_extra);
      F_BDDC_SDC.extracellular_materials(arr_extra);
      CardiacIntegrationStatistics statistics;
      F_BDDC_SDC.scaleInitialValue<0>(InitialValue(0,material,arr_excited_region),u);
      uAll = component<0>(u);

      if(options.plot) writeVTK(uAll,out+"/emiSDCBDDCInitial",
               IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),"u");

      std::cout <<" test CellFilter!!!!\n";
      std::set<int> s_temp;
      for (int i = 0; i < gridManager.grid().size(0); ++i) 
        s_temp.insert(i);

      CellFilter Cellfltr(boost::fusion::at_c<0>(u.data), cells_set, tags,material); 
      u = semiImplicit_CG_BDDC_SDC_smallest_collocation( gridManager,
                                    F_BDDC_SDC,
                                    Cellfltr,
                                    variableSetDesc,
                                    spaces,
                                    gridManager.grid(),
                                    u,
                                    i2e,
                                    options,
                                    statistics,
                                    out,
                                    uAll,
                                    i2i,
                                    cg_semi,
                                    direct,
                                    matlab_dir,
                                    sol_BDDC_SDC_small_coll,
                                    sharedDofsKaskade_new,
                                    interfaceTypes,
                                    n_subdomains,
                                    A_,
                                    M_,
                                    K_,
                                    As,
                                    Ms,
                                    Ks,
                                    IG_seq,
                                    map_IGamma,
                                    sequenceOfTags,
                                    map_II, 
                                    map_GammaGamma_noDuplicate, 
                                    weights,
                                    Fs_petcs,
                                    cg_solver,
                                    iter_cg_with_bddc,
                                    local2Global,
                                    global2Local,
                                    map_index_to_subdomain,
                                    tol,
                                    map_t2l, 
                                    map_indices,
                                    BDDC_SDC_with_initial,
                                    BDDC_verbose,
                                    BDDC_SDC_verbose);

    Vector sol_bddc_sdc_step(sol_BDDC_SDC_small_coll);
    sol_bddc_sdc_step = 0; 
    petsc_structure_rhs(sequenceOfTags, map_indices, map_II, map_GammaGamma_noDuplicate, sol_BDDC_SDC_small_coll,sol_bddc_sdc_step);
    if(write_to_file) writeSolution(sol_bddc_sdc_step,matlab_dir+"/sol_bddc_sdc_small_coll"); 

    }  
  }

  return 0;
}


