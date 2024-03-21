
#include <iostream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <bits/stdc++.h>
#include <filesystem>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <boost/timer/timer.hpp>


#include "utilities/threading.hh"
#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/gridBasics.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "fem/variables.hh"
#include "fem/norms.hh"
#include "io/vtk.hh"
#include "io/vtkreader.hh"
// #include "io/matlab.hh"
#include "linalg/direct.hh"
#include "linalg/matrixBlocks.hh"
#include "fem/spaces.hh"
#include "utilities/gridGeneration.hh"  
#include "utilities/kaskopt.hh"         
#include "timestepping/semieuler.hh"
#include "timestepping/sdc.hh"
#include "linalg/dynamicMatrix.hh"
#include "mg/bddc.hpp"
#include "utilities/timing.hh"
#include "fem/diffops/membraneModels.hh"
#include "fem/partitionedspace.hh" 
#include "linalg/jacobiPreconditioner.hh"
#include "mg/bddc.hpp"
#include "linalg/apcg.hh"
using namespace Kaskade;

#include "EMI_model.hh"
#include "EMI_mesh_utility.hpp"
#include "EMI_mesh_data.hpp"

#include "EMI_CG_Jacobi.hpp"
#include "EMI_CG_BDDC.hpp"

using namespace std;


#ifndef SPACEDIM
#define SPACEDIM 2
#endif

template <typename Material>
struct InitialValue 
{
  using Scalar = double;
  static constexpr int components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  InitialValue(int c, Material const& material_,const std::vector<int> & arr_excited_region_): component(c), material(material_) {
    arr_excited_region.assign(arr_excited_region_.begin(), arr_excited_region_.end()); 
  }
  
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> const& localCoordinate) const
  {
    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);

    constexpr int dim = SPACEDIM;
    Dune::FieldVector<double,dim> zero(0);
    int cell_mat = material.value(cell,zero);

    if(std::find(arr_excited_region.begin(), arr_excited_region.end(), cell_mat) != arr_excited_region.end()){
      return 0.5;
    } 
    return 0.0;
  }

private:
  int component;
  Material const& material;
  std::vector<int> arr_excited_region;
};

struct CardiacIntegrationOptions
{
  double T;             // final time 
  double dt;            // time step size (target)
  int maxSteps;         // upper limit on time step number
  int nCollocU;         // number of collocation points for transmembrane voltage & co
  int nCollocUstart;    // start with that many collocation points
  int verbosity;        // whether to report more or less
  int minSweeps;        // perform at least this many sweeps
  int maxSweeps;        // perform at most this many sweeps
  int sweepType;        // 0: Euler, 1: LU
  int writeVTK;         // output of solutions
  Kaskade::IoOptions::OutputType outputType; // ascii or binary
  bool adapt;           // whether to perform adaptive grid refinement
  int minRefLevel;      // don't coarsen below this level
  int assemblyThreads;  // number of assembler threads to use
  double aTol;          // absolute tolerance (L^2 norm)
  double tolSelect;     // tolerance for node selection (L^\infty norm)
  int order;           // FE order of transmembrane voltage
  bool rosenbrockRefinementStyle; // 
  std::string prefix;   // target directory for vtk output
  int maxCGIter;        // max number of IterateType::CG iterations
  double cgTol;         // IterateType::CG energy tolerance
  double tol_BDDC;      // BDDC tol
  int nReactionSweeps;  // number of extra pointwise sweep steps
  bool plot;            // flag to plot the output 
  bool CG_shift; 
  double SDC_TOL; // SDC tolerance
  double sdc_contraction; // sdc contraction
};

struct CardiacIntegrationStatistics
{
  double avgDofs;
  int    maxLevel;
  int    totalSweeps;
  double sdcTime;
};

template <class FSElement, class Material>
class CellFilter {

public:
    CellFilter(FSElement& fse_, std::set<int> marked_cells_, std::set<int> tags_, Material const& material_): fse(fse_), material(material_) {
      marked_cells.insert(marked_cells_.begin(), marked_cells_.end());
      tags.insert(tags_.begin(), tags_.end());
    }    

    typedef typename FSElement::Space ImageSpace;
    typedef typename ImageSpace::Grid Grid;

    template <class Cell>   
    bool operator()(Cell const& cell) const {
      

      auto cellIndex = fse.space().indexSet().index(cell);
      auto cell_itr = marked_cells.find(cellIndex);
      
      Dune::FieldVector<double,ImageSpace::dim> zero(0.0);
      int material_var = material.value(cell,zero);
      auto tag_itr = tags.find(material_var);

      if(based_on_tag and tag_itr != tags.end())
      {
        std::cout << "CellFilter based on tag: "<< material_var << " tag " << " tags.size() "<<tags.size()<< std::endl;
        return tag_itr != tags.end();
      }else if(based_on_tag){
          return false;
      }
    
      return cell_itr != marked_cells.end();
    }

    void set_cells(std::set<int> marked_cells_update){
      marked_cells.clear();  
      marked_cells.insert(marked_cells_update.begin(), marked_cells_update.end());
    }

    void set_tags(std::set<int> tags_){
      tags.clear();  
      tags.insert(tags_.begin(), tags_.end());
      // std::set<int>::iterator tag_itr;
      // for (tag_itr = tags.begin(); tag_itr != tags.end(); ++tag_itr) {
      //   std::cout << "set_tags: "<< *tag_itr <<std::endl;
      // }
    }

    void select_based_on_tag(bool based_on_tag_)
    {
      based_on_tag = based_on_tag_;
    }

    void get_cells(){
      std::set<int>::iterator it_cell;
      for (it_cell = marked_cells.begin(); it_cell != marked_cells.end(); ++it_cell) {
          std::cout << " cells: "<<*it_cell <<std::endl;
      }
    }
    int get_size(){
      return marked_cells.size();
    }
private:
    FSElement fse;
    std::set<int> marked_cells;
    Material const& material;
    std::set<int> tags;
    bool based_on_tag = false;
};

template <class State, class elementType>
void printuAll(State const& x, elementType & uAll, int order, std::string filename, std::string varname)
{

  uAll*=0;
  uAll = component<0>(x);
  writeVTK(uAll,filename,
                    IoOptions().setOrder(order).setPrecision(7).setDataMode(IoOptions::nonconforming),varname);
}

void writeDoubleTofile(double b, std::string name){

  std::string fname = name+".m";
  std::ofstream f(fname.c_str());
  f.precision(16);

  f << "function [b] = " << "solution" << '\n' << " b = \n";
  f << b;

  f << "\n;";
}

void writeVectorTofile(auto b, std::string name){

  std::string fname = name+".m";
  std::ofstream f(fname.c_str());
  f.precision(16);

  int size_rhs = component<0>(b).size();

  f << "function [b] = " << "solution" << '\n'
    << " b = [\n";

  for (int k = 0; k < size_rhs; ++k)
  {
   f << std::setprecision(16) << component<0>(b)[k] <<'\n';
  }
  
  f << "\n];";
}

template <class VEntry>
void writeSolution(Dune::BlockVector<VEntry> const& b, std::string const& basename, int precision=16)
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

template <class Matrix, class Vector>
void writeMatrixVectorTofile(Matrix & A, Vector & b, std::string name){

  std::string fname = name+".m";
  std::ofstream f(fname.c_str());
  f.precision(16);

  MatrixAsTriplet mmat(A);

  f << "function [A,b] = " << "sdcMatrix" << '\n'
  << " b = [\n";

  for (int k = 0; k < b.size(); ++k)
  {
    f << b[k] <<'\n';
  }
  f << "\n];\ndata = [\n";

  for (size_t i=0; i<mmat.ridx.size(); ++i)
  {
    f << mmat.ridx[i]+1 << ' ' << mmat.cidx[i]+1 << ' ' << mmat.data[i] << '\n';
  }
  f << "];\nA = sparse(data(:,1),data(:,2),data(:,3)," << A.N() << "," << A.M() << ");\n";
}
