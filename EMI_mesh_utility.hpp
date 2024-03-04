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
 * \param map_IGamma for each material, gives the all the global indices of the subdomain
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
                                           std::map<int, int> & map_t2l,
                                           std::map<int,std::set<int>> & map_IGamma)
{

  typedef typename FSElement::Space ImageSpace;
  typedef typename ImageSpace::Grid Grid;

  DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, ImageSpace::sfComponents, 1>> globalValues;

  fse.coefficients() = typename ImageSpace::Scalar(0.0);

  typename ImageSpace::Evaluator isfs(fse.space()); // evalautor of finite space element

  auto const cend = fse.space().gridView().template end<0>(); //  cell end

  using ValueType = decltype(fu.value(*cend,Dune::FieldVector<typename Grid::ctype, ImageSpace::dim>()));
  
  std::map<int,int>::iterator it;
  std::map<int,std::set<int>>::iterator it_igamma;

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

        // adding all the igamma for matreial 
        it_igamma = map_IGamma.find(material_var);
        if(it_igamma!= map_IGamma.end()){
          std::set<int> igamma = it_igamma->second;
          igamma.insert(nIndex);
          it_igamma->second = igamma;
        }else{
          std::set<int> igamma;
          igamma.insert(nIndex);
          map_IGamma[material_var] = igamma;
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
                                 std::vector<std::set<int>> & i2i,
                                 std::map<int,std::set<int>> & map_GammaGamma,
                                 std::map<int,std::set<int>> & map_GammaGamma_W_Nbr,
                                 std::map<int,std::map<int,std::set<int>>> & map_GammaNbr) 
{
  typedef typename FSElement::Space ImageSpace;
  typedef typename ImageSpace::Grid Grid;

  DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, ImageSpace::sfComponents, 1>> globalValues;

  fse.coefficients() = typename ImageSpace::Scalar(0.0);

  typename ImageSpace::Evaluator isfs(fse.space());

  auto const cend = fse.space().gridView().template end<0>();
  std::map<int,std::set<int>>::iterator it_gamma;
  std::map<int,std::set<int>>::iterator it_gammaNbr;

  std::map<int,std::map<int,std::set<int>>>::iterator it_GammaNbr;
  std::map<int,std::set<int>>::iterator it_GammaNbr_i;

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
          
          int material_nbr = material.value(intersection.outside(),zero);

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

                  // adding all the gammagamma for matreial 
                  it_gamma = map_GammaGamma.find(material_var);
                  if(it_gamma!= map_GammaGamma.end()){
                    std::set<int> gamma = it_gamma->second;
                    gamma.insert(nIndex_c1);
                    it_gamma->second = gamma;

                    // gamma + nbr
                    it_gammaNbr = map_GammaGamma_W_Nbr.find(material_var);
                    std::set<int> gammaNbr = it_gammaNbr->second;
                    gammaNbr.insert(nIndex_c1);
                    gammaNbr.insert(nIndex_c2);
                    it_gammaNbr->second = gammaNbr;

                    // nbr
                    it_GammaNbr = map_GammaNbr.find(material_var);
                    std::map<int,std::set<int>> gamma_nbr = it_GammaNbr->second;
                    it_GammaNbr_i = gamma_nbr.find(material_nbr);
                    if(it_GammaNbr_i!= gamma_nbr.end()){
                      std::set<int> gamma_nbr_i = it_GammaNbr_i->second;
                      gamma_nbr_i.insert(nIndex_c2);
                      it_GammaNbr_i->second = gamma_nbr_i;
                      it_GammaNbr->second = gamma_nbr;
                    }else{
                      std::set<int> gamma_nbr_i;
                      gamma_nbr_i.insert(nIndex_c2);
                      gamma_nbr[material_nbr] = gamma_nbr_i;
                      it_GammaNbr->second = gamma_nbr; 
                    }

                  }else{
                    std::set<int> gamma;
                    gamma.insert(nIndex_c1);
                    map_GammaGamma[material_var] = gamma;

                    // gamma + nbr
                    std::set<int> gammaNbr;
                    gammaNbr.insert(nIndex_c1);
                    gammaNbr.insert(nIndex_c2);
                    map_GammaGamma_W_Nbr[material_var] = gammaNbr;

                    // nbr
                    std::map<int,std::set<int>> gamma_nbr;
                    std::set<int> gamma_nbr_value;
                    gamma_nbr_value.insert(nIndex_c2);
                    gamma_nbr[material_nbr] = gamma_nbr_value;
                    map_GammaNbr[material_var] = gamma_nbr;

                  }
                  break;
                }                
              }
              else if(ImageSpace::dim==3){
                if(icoord[nIndex_c1][0]==icoord[nIndex_c2][0] and
                   icoord[nIndex_c1][1]==icoord[nIndex_c2][1] and
                   icoord[nIndex_c1][2]==icoord[nIndex_c2][2] ){
                  s_index.insert(nIndex_c2);

                  // adding all the gammagamma for matreial 
                  it_gamma = map_GammaGamma.find(material_var);
                  if(it_gamma!= map_GammaGamma.end()){
                    std::set<int> gamma = it_gamma->second;
                    gamma.insert(nIndex_c1);
                    it_gamma->second = gamma;

                    // gamma + nbr
                    it_gammaNbr = map_GammaGamma_W_Nbr.find(material_var);
                    std::set<int> gammaNbr = it_gammaNbr->second;
                    gammaNbr.insert(nIndex_c1);
                    gammaNbr.insert(nIndex_c2);
                    it_gammaNbr->second = gammaNbr;

                    // nbr
                    it_GammaNbr = map_GammaNbr.find(material_var);
                    std::map<int,std::set<int>> gamma_nbr = it_GammaNbr->second;
                    it_GammaNbr_i = gamma_nbr.find(material_nbr);
                    if(it_GammaNbr_i!= gamma_nbr.end()){
                      std::set<int> gamma_nbr_i = it_GammaNbr_i->second;
                      gamma_nbr_i.insert(nIndex_c2);
                      it_GammaNbr_i->second = gamma_nbr_i;
                      it_GammaNbr->second = gamma_nbr;
                    }else{
                      std::set<int> gamma_nbr_i;
                      gamma_nbr_i.insert(nIndex_c2);
                      gamma_nbr[material_nbr] = gamma_nbr_i;
                      it_GammaNbr->second = gamma_nbr;
                    }

                  }else{
                    std::set<int> gamma;
                    gamma.insert(nIndex_c1);
                    map_GammaGamma[material_var] = gamma;

                    // gamma +nbr
                    std::set<int> gammaNbr;
                    gammaNbr.insert(nIndex_c1);
                    gammaNbr.insert(nIndex_c2);
                    map_GammaGamma_W_Nbr[material_var] = gammaNbr;

                    // nbr
                    std::map<int,std::set<int>> gamma_nbr;
                    std::set<int> gamma_nbr_value;
                    gamma_nbr_value.insert(nIndex_c2);
                    gamma_nbr[material_nbr] = gamma_nbr_value;
                    map_GammaNbr[material_var] = gamma_nbr;  

                  }
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

template< class FSElement, class Function>
void markedIndicesForDirichlet(FSElement& fse,  
                                 Function const& fu, 
                                 std::vector<std::vector<int>> & cell2Indice, 
                                 std::set<int> & dofsDiriichlet )
{
  typedef typename FSElement::Space ImageSpace;
  typedef typename ImageSpace::Grid Grid;

  DynamicMatrix< Dune::FieldMatrix<typename ImageSpace::Scalar, ImageSpace::sfComponents, 1>> globalValues;

  fse.coefficients() = typename ImageSpace::Scalar(0.0);

  typename ImageSpace::Evaluator isfs(fse.space());

  auto const cend = fse.space().gridView().template end<0>();
  //std::cout <<  " cend localCoordinate! *ci " <<std::endl;

  using ValueType = decltype(fu.value(*cend,Dune::FieldVector<typename Grid::ctype, ImageSpace::dim>()));
  std::vector<ValueType> fuvalue; // declare here to prevent reallocations
  for (auto ci=fse.space().gridView().template begin<0>(); ci!=cend; ++ci)
  {
    auto cellIndex = fse.space().indexSet().index(*ci);
    std::cout << cellIndex <<std::endl;
    isfs.moveTo(*ci);

    auto const& localCoordinate(isfs.shapeFunctions().interpolationNodes());
    globalValues.setSize(localCoordinate.size(),1);

    using Cell = decltype(ci);
    auto dof_u = fse.space().mapper().globalIndices(*ci);
    int nrNodes = dof_u.size();

    for(auto const& intersection : intersections(fse.space().gridView(),*ci)){
      if(!intersection.neighbor()){
        for (int i = 0; i < cell2Indice[cellIndex].size(); ++i)
        {

          int index_c1 = cell2Indice[cellIndex][i];
          dofsDiriichlet.insert(index_c1);
        }
      }
    }
  }
}


#endif
