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

#ifndef EMIMODEL_HH
#define EMIMODEL_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "fem/fixdune.hh"
#include "fem/functional_aux.hh"

#include <algorithm>
#include "fem/variables.hh"
#include "utilities/linalg/scalarproducts.hh"


// A hash function used to hash a pair of any kind
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
      auto hash1 = std::hash<T1>{}(p.first);
      auto hash2 = std::hash<T2>{}(p.second);

      if (hash1 != hash2) {
          return hash1 ^ hash2;             
      }
       
      // If hash1 == hash2, their XOR is zero.
        return hash1;
    }
};


template <class Scalar_, class VarSet, class Material, class Grid_, class SPACE, class MembraneModel>
class EMI_model : public FunctionalBase<WeakFormulation>//Kaskade::VariationalFunctional>//WeakFormulation>
{
  using Self=EMI_model<Scalar_,VarSet,Material,Grid_,SPACE,MembraneModel>;
public:
  using Scalar = Scalar_;
  using AnsatzVars = VarSet;
  using OriginVars = VarSet;
  using TestVars = VarSet;

  using GridView = typename AnsatzVars::GridView;

  static constexpr int dim = AnsatzVars::Grid::dimension;

  static constexpr int ueIdx = 100;
  static constexpr int uIdx = 0;
  static constexpr int uSpaceIdx = spaceIndex<AnsatzVars,uIdx>;

  // convenient template aliases
  template <int row> using TestComponents = std::integral_constant<size_t,TestVars::template Components<row>::m>;
  template <int row> using AnsatzComponents = std::integral_constant<size_t,AnsatzVars::template Components<row>::m>;

 typedef MembraneModel  Membrane;

  class DomainCache
  {
  public:
    DomainCache(EMI_model const& F_, typename AnsatzVars::VariableSet const& vars_,int flags=7):
      F(F_), 
      vars(vars_)
    {
    }

    template <class Cell>
    void moveTo(Cell const& cell)
    {
      Dune::FieldVector<double,dim> zero(0);
      cell_mat = F.material.value(cell,zero);
      extra_cell = false;
      if(std::find(F.arr_extra.begin(), F.arr_extra.end(), cell_mat) != F.arr_extra.end()) 
        extra_cell = true;
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;
      u  = component<uIdx>(vars).value(at_c<uSpaceIdx>(evaluators));
      du  = component<uIdx>(vars).derivative(at_c<uSpaceIdx>(evaluators));
      f = 0.0;
    }
    Scalar d0() const
    {
      return 0.0;
    }

    template<int row>
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>  
    d1(VariationalArg<Scalar,dim,TestComponents<row>::value> const& argT) const
    {

      if(extra_cell)
      {
        return -F.sigma_e*sp(du,argT.derivative) + f*argT.value;
      }
      else
      {
        return -F.sigma_i*sp(du,argT.derivative) + f*argT.value;
      }

      // never get here
      return 0.0;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m>
    d2(VariationalArg<Scalar,dim,TestComponents<row>::value> const& argT,
       VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& argA) const
    {
      if (row != col) // subdomains do not couple in the domain
        return 0;

      // Remember that right hand side is Laplace operator div(sigma*nabla u), and this is,
      // after multiplication with test function phi and integration by parts, represented by
      // - nabla phi * sigma * nabla u. Note the minus sign!

      Scalar p = -1;
      
      if( F.mass==1){
        p = 0;    
      } 
      if(F.mass>1)
        p = 1;

      if(extra_cell)  
      {
        return -p*F.sigma_e*sp(argT.derivative,argA.derivative);
      }
      else
      {
        return -p*F.sigma_i*sp(argT.derivative,argA.derivative);
      }

      // never get here
      return 0.0;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m> 
    b2(Kaskade::VariationalArg<Scalar,dim> const &argT,
       Kaskade::VariationalArg<Scalar,dim> const &argA) const
    {
      return 0.0;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> u, f;
    Dune::FieldMatrix<Scalar,AnsatzComponents<uIdx>::value,dim> du;
    Kaskade::LinAlg::EuclideanScalarProduct sp;
    int cell_mat;
    bool extra_cell;
  };

  class BoundaryCache
  {
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;
  public:
    BoundaryCache(EMI_model const& F_, typename AnsatzVars::VariableSet const& vars_,int flags=7)
    : F(F_), vars(vars_), e(nullptr), gamma(1e-5)
    {
    }

    void moveTo(FaceIterator const& cell)
    {
      e = &cell;

      Dune::FieldVector<double,dim> zero(0);
      cell_mat = F.material.value((*cell).inside(),zero);
      extra_cell = false;
      if(std::find(F.arr_extra.begin(), F.arr_extra.end(), cell_mat) != F.arr_extra.end())
        extra_cell = true;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators) 
    {
      using namespace boost::fusion;      
      ue = component<uIdx>(vars).value(at_c<uSpaceIdx>(evaluators));
      ue0 = 0;
    }
    
    template<int row>
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m> 
    d1(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      if(extra_cell) 
      {
        return -gamma*(ue-ue0)* arg.value;
      }

      return 0.0;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m>
    d2(VariationalArg<Scalar,dim,TestComponents<row>::value> const &argT, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &argA) const
    {
      double p = -1;

      if(F.mass>1)
        p = 1;

      if(F.mass==1)
        return 0;

      if(extra_cell) 
      {
        return -p*gamma*argT.value*argA.value;
      }

      return 0.0;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m> 
    b2(Kaskade::VariationalArg<Scalar,dim> const &argT,
       Kaskade::VariationalArg<Scalar,dim> const &argA) const
    {
      return 0.0;
    }

  private:
    EMI_model const& F;
    typename AnsatzVars::VariableSet const& vars;
    FaceIterator const* e;
    Scalar gamma;
    Dune::FieldVector<Scalar,AnsatzComponents<uIdx>::value> ue, ue0;
    int cell_mat;
    bool extra_cell;
  };

  class InnerBoundaryCache 
  {
  public:
    InnerBoundaryCache(EMI_model const& F_,
                       typename AnsatzVars::VariableSet const& vars_,
                       int flags=7):
      F(F_),
      vars(vars_), penalty(1), order(1),
      C_m(F_.C_m), 
      R(F_.R)
    {
    }

    template <class FaceIterator>
    void moveTo(FaceIterator const& f)
    {
      face = *f;      
      int in = (F.grid).leafIndexSet().index(face.inside());
      int out = (F.grid).leafIndexSet().index(face.outside());
      // w = F.set_gate_variable()[std::make_pair(in, out)];

      Dune::FieldVector<double,dim> zero(0);
      cellDomain = F.material.value(face.inside(),zero);
      neighbourDomain = F.material.value(face.outside(),zero);
    }
    
    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators, Evaluators const& neighbourEvaluators)
    {

      using namespace boost::fusion;

      // Get the subdomain of the current cell and the neighbour cell. Ideally, this would be done in
      // moveTo(), but we query the support (i.e. which potential has dofs on the current cell)
      // simply by asking the evaluators, which are only available in evaluateAt(). Other implementations
      // not relying on evaluators are of course possible.
 
      if (cellDomain==neighbourDomain)    // no membrane - skip
        return;

      // std::cout << "on membrane position " << face.geometry().global(x)[0]
      //           << ", outward normal " << face.centerUnitOuterNormal() << "\n";

      // Extract the potential in our cell and the neighbour cell.
      Scalar cellU = vars.template value<uIdx>(evaluators);
      Scalar neighbourU = vars.template value<uIdx>(neighbourEvaluators);

      // Initialize transmembrane voltage as our potential minus neighbour potential.
      Scalar v = cellU - neighbourU;

      two_extra = false;
      if (std::find(F.arr_extra.begin(), F.arr_extra.end(), cellDomain) != F.arr_extra.end() and 
          std::find(F.arr_extra.begin(), F.arr_extra.end(), neighbourDomain) != F.arr_extra.end()) // two extracellular regions
      {
        two_extra = true;
      }

      memInterface = false;
      if ((std::find(F.arr_extra.begin(), F.arr_extra.end(), cellDomain) != F.arr_extra.end() and 
           std::find(F.arr_extra.begin(), F.arr_extra.end(), neighbourDomain) == F.arr_extra.end()) || 
          (std::find(F.arr_extra.begin(), F.arr_extra.end(), cellDomain) == F.arr_extra.end() and 
           std::find(F.arr_extra.begin(), F.arr_extra.end(), neighbourDomain) != F.arr_extra.end())) // membrane
      {
        memInterface = true;
      }


      if (cellDomain > neighbourDomain)
        std::tie(current,dcurrentCell,dcurrentNeighbour) = ionCurrent(v);
      else
      {
        std::tie(current,dcurrentCell,dcurrentNeighbour) = ionCurrent(-v);
        current = -current;
      }
    }

    std::tuple<Scalar,Scalar,Scalar> ionCurrent(Scalar v) const
    {
      // two extra celllular neighbors
      if (two_extra)
      {
        v = 0;
        // compute fake gap junction between two extra cellluar regions
        Scalar gap_junc = F.GapJunctionLinearExtraCell(v);
        Scalar dgap_junc = F.dGapJunctionLinearExtraCell();
        return std::make_tuple(gap_junc,0,0);
      }

      // nonlinear ion current
      if (memInterface) // membrane
      {
        // w = 0;
        // membrane current properly
        // ionic model in kaskade implemented by -1
        Scalar Ionic = -F.membrane().current(v,0);//w
        Scalar dIonic = -F.membrane().current_du(v,0);//w

        // ignore the derivative of ionic from the lhs
        // return std::make_tuple(Ionic, dIonic, -dIonic);
        return std::make_tuple(Ionic, 0, 0);
      }

      // compute gap junction ->  linear gap junction  
      Scalar gap_junc = F.GapJunctionLinear(v);
      Scalar dgap_junc = F.dGapJunctionLinear();

      // ignore the derivative of ionic from the lhs
      return std::make_tuple(gap_junc,0,0);
    }

    template<int row>
    Dune::FieldVector<Scalar,1>
    d1(Kaskade::VariationalArg<Scalar,dim> const& argT) const
    {
      // Test function support is limited to the own subdomain, i.e. in the extracellular domain,
      // there are no intracellular test function degrees of freedom - check this.
      assert(row==cellDomain);

      // Integration only over membrane - if cells on both sides of the face lie in the same
      // subdomain, we have nothing to do.
      if (cellDomain==neighbourDomain)
        return 0;

      // This is the flux integral on the membrane.
      auto result = -argT.value * current;

      return result;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,1,1>
    d2(Kaskade::VariationalArg<Scalar,dim> const &argT,
       Kaskade::VariationalArg<Scalar,dim> const &argA, bool centerCell) const
    {
      assert(row==cellDomain);

      // As in d1, nothing to do if the face is not part of the membrane
      if (cellDomain==neighbourDomain)
        return 0;

      if(F.mass==1)
        return 0; // only mass

      Scalar dCurrent = centerCell? dcurrentCell: dcurrentNeighbour;

      
      double p = -1;
      
      if(F.mass>1)
        p = 1;

      // Compute derivative of the flux integrand.
      Scalar result = -p * argT.value * dCurrent * argA.value;

      return result;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m> 
    b2(Kaskade::VariationalArg<Scalar,dim> const &argT,
       Kaskade::VariationalArg<Scalar,dim> const &argA, bool centerCell) const
    {
      assert(row==cellDomain);
  
      if (cellDomain==neighbourDomain)
        return 0;

      if(two_extra)
        return 0; // only mass

      if( F.mass==0){
        return 0; // only stiffness
      } 
      
      if(F.mass_submatrix){
        if((row == col) and (cellDomain==F.row_submatrix) and (neighbourDomain==F.col_submatrix)  )
        {
          // std::cout << " F.row "<< F.row_submatrix << " F.col" << F.col_submatrix <<std::endl;
          Scalar sign = centerCell? 1: -1;

          // Compute mass matrix
          Scalar result = C_m * sign * argT.value * argA.value;    
          return result;
        }
        return 0;
      }
        
      Scalar sign = centerCell? 1: -1;

      // Compute mass matrix
      Scalar result = C_m * sign * argT.value * argA.value;    
      return result;
    }

  private:
    EMI_model const& F;
    typename AnsatzVars::VariableSet const& vars;
    Scalar penalty, C_m, eta_0, eta_1, eta_2, eta_3, v_th, v_pk, R;
    Scalar v1;
    Scalar v2;
    Scalar v12;
    Scalar w = 0;
    int order = 1;
    int cellDomain, neighbourDomain;
    Scalar current, dcurrentCell, dcurrentNeighbour;
    Face<GridView> face;
    bool memInterface;
    bool two_extra;
  };

public:
  template <int row, int col>
  struct D2: public FunctionalBase<VariationalFunctional>::D2<row,col>
  {
    static constexpr bool present = true;
    static constexpr bool symmetric = true;
  };

  template <int row, int col>
  struct B2: public FunctionalBase<VariationalFunctional>::B2<row,col>
  {
    static bool const present   = true;
    static bool const symmetric = true;
    static bool const constant  = false;
  };


  /**
   * Given an initial value, this is transfered to a properly scaled
   * finite element iterate.
   */
  template <int row, class WeakFunctionView>
  void scaleInitialValue(WeakFunctionView const& us, typename AnsatzVars::VariableSet& u) const 
  {
    interpolateGloballyWeak<Volume>(boost::fusion::at_c<row>(u.data),
                                    ScaledFunction<WeakFunctionView>(row==uIdx,us,*this));
  }

  template <class WeakFunctionView>
  struct ScaledFunction 
  {
    using Scalar = typename WeakFunctionView::Scalar;
    static int const components = WeakFunctionView::components;
    using ValueType = Dune::FieldVector<Scalar,components>;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& us_,
                   EMI_model<Scalar_,VarSet,Material,Grid_,SPACE,MembraneModel> const& f_): doScaling(doScaling_), us(us_), f(f_) {}

    template <class Cell>
    int order(Cell const&) const { return std::numeric_limits<int>::max(); }

    template <class Cell>
    ValueType value(Cell const& cell,
                    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
    {
      ValueType u = 0;
      if(doScaling)
        u = us.value(cell,localCoordinate);
      return u;
    }

  private:
    bool doScaling;
    WeakFunctionView const& us;
    EMI_model<Scalar_,VarSet,Material,Grid_,SPACE,MembraneModel> const& f;
  };

  bool considerFace(Face<GridView> const& face) const
  {
    Dune::FieldVector<Scalar,dim> zero(0.0);
    return material.value(face.inside(),zero) != material.value(face.outside(),zero);
  }

  Scalar GapJunctionLinear(Scalar v) const
  { 
    return v/R; 
  }

  Scalar dGapJunctionLinear() const
  { 
    return 1/R; 
  }

  Scalar GapJunctionLinearExtraCell(Scalar v) const
  { 
    return v/R_extra; 
  }

  Scalar dGapJunctionLinearExtraCell() const
  { 
    return 1/R_extra; 
  }

  Scalar time() const { return t; }
  void time(Scalar tnew)  { t=tnew; }

  void Mass_stiff(int mass_)
  {
    mass = mass_;
  }

  void set_mass_submatrix(bool mass_submatrix_)
  {
    mass_submatrix = mass_submatrix_;
  }

  void set_row_col_subdomain(int row, int col)
  {
    row_submatrix = row;
    col_submatrix = col;
  }

  Membrane const& membrane() const { return memb; }

  void set_gate_variable( std::unordered_map<std::pair<int, int>, double, hash_pair> & map_faces_gate_gate_)
  {
    map_faces_gate_gate = map_faces_gate_gate_;
  }

  std::unordered_map<std::pair<int, int>, double, hash_pair> set_gate_variable() const { return map_faces_gate_gate; }

  void extracellular_materials( std::vector<int> &arr_extra_ )
  {
    arr_extra.assign(arr_extra_.begin(), arr_extra_.end());  
  }

  EMI_model( Material const& material_,
                  Grid_ const& grid_,
                  SPACE const& space_,
                  double penalty_=1e6,
                  double sigma_i_=1.7,
                  double sigma_e_=3.0,
                  double C_m_=1e-2, 
                  double R_=0.015,
                  double R_extra_=1e-9)
  : penalty(penalty_), 
    C_m(C_m_), 
    sigma_i(sigma_i_),
    sigma_e(sigma_e_),
    R(R_),
    R_extra(R_extra_),
    material(material_),
    grid(grid_),
    space(space_)
  {}

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder-2;
  }
private:
  double penalty, sigma_i, sigma_e, C_m, t, R, R_extra;
  bool linear_ionic;
  Material const& material;
  int mass = 10; // by default set it to any nomber bigger than 1, to get mass matrix set it to 1, to get stiffness set it to 0 
  bool mass_submatrix = false;
  int row_submatrix;
  int col_submatrix;
  Grid_ const& grid;
  SPACE const& space;

  MembraneModel memb;
  std::unordered_map<std::pair<int, int>, double, hash_pair> map_faces_gate_gate;
  std::vector<int> arr_extra;
};
#endif
