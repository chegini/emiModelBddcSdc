#ifndef EMI_MESH_UTILITY_HH
#define EMI_MESH_UTILITY_HH

void getSubdomain(std::vector<int> & arr, std::ifstream & f){
  int n;
  int count = 0;
  while (count<arr.size())
  {
    f >> n;
    if(f.eof()) break;
    arr[count] = n;
    count+=1;
  }
}


struct GetGlobalCoordinate
{
   using Scalar = double;
   static int const components = 1;
   using ValueType = Dune::FieldVector<Scalar,components>;

   GetGlobalCoordinate(){}
   template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }

   template <class Cell>
   Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> value(Cell const& cell,Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
   {
    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
    return x;
   }
};


/**
 * \brief extract data of EMI material.
 *
 * \param fu global coordinate function
 * \param material data for each material sets with tag number
 * \param e2i element to indices
 * \param i2e index to element
 * \param i2i index to indices of other region with different matrial intersection
 * \param icoord index to global coordinates
 * \param i2T index to Tag
 * \param map_t2l tag to length
 */

template< class FSElement, class Function, class Material>
void getInnerInterfaceDofsForeachSubdomain(FSElement& fse,  
                                           Function const& fu, 
                                           Material const & material, 
                                           std::vector<std::vector<int>> & e2i, 
                                           std::vector<std::set<int>> & i2e,
                                           std::vector<std::set<int>> & i2i,
                                           std::vector<std::vector<double>> & icoord,
                                           std::vector<int> & i2T,
                                           std::map<int, int> & map_t2l  )
{

  typedef typename FSElement::Space ImageSpace;
  typedef typename ImageSpace::Grid Grid;

  DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, ImageSpace::sfComponents, 1>> globalValues;

  fse.coefficients() = typename ImageSpace::Scalar(0.0);

  typename ImageSpace::Evaluator isfs(fse.space()); // evalautor of finite space element

  auto const cend = fse.space().gridView().template end<0>(); //  cell end

  using ValueType = decltype(fu.value(*cend,Dune::FieldVector<typename Grid::ctype, ImageSpace::dim>()));
  
  std::map<int,int>::iterator it;

  // iterate over cells
  for (auto ci=fse.space().gridView().template begin<0>(); ci!=cend; ++ci)
  {
    auto eIndex = fse.space().indexSet().index(*ci); // get cell index
    isfs.moveTo(*ci); 

    auto const& localCoordinate(isfs.shapeFunctions().interpolationNodes());
    globalValues.setSize(localCoordinate.size(),1); // not used!


    using Cell = decltype(ci);
    auto dof_u = fse.space().mapper().globalIndices(*ci);
    int nrNodes = dof_u.size();

    Dune::FieldVector<double,ImageSpace::dim> zero(0.0);
    int material_var = material.value(*ci,zero);

    // iterate over nodes of each cell
    for (int i = 0; i < isfs.globalIndices().size(); ++i) 
    {
      int nIndex = isfs.globalIndices()[i];
      e2i[eIndex].push_back(nIndex); // e2n
      std::set<int> s_index = i2i[nIndex];
      s_index.insert(nIndex);
      i2i[nIndex] = s_index; 

      std::set<int> cell_indices = i2e[nIndex];
      cell_indices.insert(eIndex);
      i2e[nIndex] = cell_indices;

      auto x = fu.value(*ci,localCoordinate[i]);

      if(icoord[nIndex].size()==0){
        i2T[nIndex] = material_var;
        // count the number of dof for each subdomain
        it = map_t2l.find(material_var);
        if (map_t2l[material_var]!=0){
          int old = it->second;
          it->second = old+1;
        }else{
          map_t2l[material_var] = 1;
        }

        for (int j = 0; j < x.size(); ++j){
          icoord[nIndex].push_back(x[j]);
        }        
      }
    }
  }
}

/**
 * \brief extract data of EMI material.
 *
 * \param fu global coordinate function
 * \param material data for each material sets with tag number
 * \param icoord index to global coordinates
 * \param e2i element to indices
 * \param e2e element to element
 * \param i2i update index to indices of other region with different matrial intersection
 */

template< class FSElement, class Function, class Material>
void markedIndicesOnInterfacesForeachSubdomain(FSElement& fse,  
                                 Function const& fu, 
                                 Material const & material, 
                                 std::vector<std::vector<double>> & icoord,
                                 std::vector<std::vector<int>> & e2i, 
                                 std::vector<std::set<int>> & e2e,
                                 std::vector<std::set<int>> & i2i) 
{
  typedef typename FSElement::Space ImageSpace;
  typedef typename ImageSpace::Grid Grid;

  DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, ImageSpace::sfComponents, 1>> globalValues;

  fse.coefficients() = typename ImageSpace::Scalar(0.0);

  typename ImageSpace::Evaluator isfs(fse.space());

  auto const cend = fse.space().gridView().template end<0>();

  using ValueType = decltype(fu.value(*cend,Dune::FieldVector<typename Grid::ctype, ImageSpace::dim>()));
  std::vector<ValueType> fuvalue; // declare here to prevent reallocations
  for (auto ci=fse.space().gridView().template begin<0>(); ci!=cend; ++ci)
  {
    auto eIndex = fse.space().indexSet().index(*ci);
    isfs.moveTo(*ci);

    auto const& localCoordinate(isfs.shapeFunctions().interpolationNodes());
    globalValues.setSize(localCoordinate.size(),1);

    using Cell = decltype(ci);
    auto dof_u = fse.space().mapper().globalIndices(*ci);
    int nrNodes = dof_u.size();

    Dune::FieldVector<double,ImageSpace::dim> zero(0.0);
    int material_var = material.value(*ci,zero);

    std::set<int> s = e2e[eIndex];
    s.insert(eIndex);
    e2e[eIndex] = s; 

    for(auto const& intersection : intersections(fse.space().gridView(),*ci)){
      if(intersection.neighbor()){
        int eNbrIndex = fse.space().gridView().indexSet().index(intersection.outside());  
        if(material.value(*ci,zero)!=material.value(intersection.outside(),zero)) {

          for (int i = 0; i < e2i[eIndex].size(); ++i)
          {
            int nIndex_c1 = e2i[eIndex][i];
            std::set<int> s_index = i2i[nIndex_c1];
            for (int j = 0; j < e2i[eNbrIndex].size(); ++j)
            {
              int nIndex_c2 = e2i[eNbrIndex][j];

              if(ImageSpace::dim==2){
                if(icoord[nIndex_c1][0]==icoord[nIndex_c2][0] and
                   icoord[nIndex_c1][1]==icoord[nIndex_c2][1] ){
                  s_index.insert(nIndex_c2);
                  break;
                }                
              }
              else if(ImageSpace::dim==3){
                if(icoord[nIndex_c1][0]==icoord[nIndex_c2][0] and
                   icoord[nIndex_c1][1]==icoord[nIndex_c2][1] and
                   icoord[nIndex_c1][2]==icoord[nIndex_c2][2] ){
                  s_index.insert(nIndex_c2);
                  break;
                }                
              }
            }
            i2i[nIndex_c1] = s_index;
          }
          std::set<int> s = e2e[eIndex];
          s.insert(eNbrIndex);
          e2e[eIndex] = s; 
        }
      }
    }
  }
}

#endif
