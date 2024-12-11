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


// Function template
template< class FSElement, class Function, class Material>
void getInnerInterfaceDofsForeachSubdomain(FSElement& fse,  
                                           Function const& fu, 
                                           Material const & material, 
                                           std::vector<std::vector<int>> & e2i, 
                                           std::vector<std::set<int>> & i2e,
                                           std::vector<std::set<int>> & i2t,
                                           std::vector<std::set<int>> & i2i,
                                           std::map<std::pair<int, int>, std::vector<double>> & coord,
                                           std::map<int, std::vector<double>> & coord_globalIndex, 
                                           std::vector<int> & i2T,
                                           std::map<int,std::set<int>> & map_IGamma, 
                                           int number_elem)
{
    typedef typename FSElement::Space ImageSpace;
    typedef typename ImageSpace::Grid Grid;

    // Thread-safe variables
    std::mutex mtx_e2i, mtx_i2e, mtx_i2t, mtx_i2i, mtx_coord, mtx_coord_globalIndex, mtx_i2T, mtx_map_IGamma;

    // Number of threads
    const int numThreads = std::thread::hardware_concurrency();
    const int cellsPerThread = (number_elem + numThreads - 1) / numThreads;

    // Thread function
    auto processCells = [&](int startIdx, int endIdx) {
        // Thread-local storage
        std::unordered_map<int, std::set<int>> local_i2e, local_i2t, local_i2i;
        std::unordered_map<std::pair<int, int>, std::vector<double>, 
                           boost::hash<std::pair<int, int>>> local_coord;
        std::unordered_map<int, std::vector<double>> local_coord_globalIndex;
        std::unordered_map<int, int> local_i2T;
        std::unordered_map<int, std::set<int>> local_map_IGamma;

        typename ImageSpace::Evaluator isfs(fse.space());
        auto const cend = fse.space().gridView().template end<0>();

        // Iterate over assigned cells
        for (int cellIdx = startIdx; cellIdx < endIdx; ++cellIdx) {
            auto ci = fse.space().gridView().template begin<0>();
            std::advance(ci, cellIdx);

            auto eIndex = fse.space().indexSet().index(*ci);
            isfs.moveTo(*ci);
            auto const& localCoordinate = isfs.shapeFunctions().interpolationNodes();

            using Cell = decltype(ci);
            auto dof_u = fse.space().mapper().globalIndices(*ci);
            int nrNodes = dof_u.size();

            Dune::FieldVector<double, ImageSpace::dim> zero(0.0);
            int material_var = material.value(*ci, zero);

            for (int i = 0; i < isfs.globalIndices().size(); ++i) {
                int nIndex = isfs.globalIndices()[i];
                std::pair<int, int> pairs = {nIndex, material_var};

                // Update local e2i
                {
                    std::lock_guard<std::mutex> lock(mtx_e2i);
                    e2i[eIndex].push_back(nIndex);
                }

                // Update i2i
                local_i2i[nIndex].insert(nIndex);

                // Update i2e
                local_i2e[nIndex].insert(eIndex);

                // Update i2t
                local_i2t[nIndex].insert(material_var);

                auto x = fu.value(*ci, localCoordinate[i]);

                // Update coord_globalIndex
                if (local_coord_globalIndex.find(nIndex) == local_coord_globalIndex.end()) {
                    for (double val : x) {
                        local_coord_globalIndex[nIndex].push_back(val);
                    }
                }

                // Update coord
                if (local_coord.find(pairs) == local_coord.end()) {
                    local_i2T[nIndex] = material_var;
                    for (double val : x) {
                        local_coord[pairs].push_back(val);
                    }

                    // Update map_IGamma
                    local_map_IGamma[material_var].insert(nIndex);
                }
            }
        }

        // Merge results into global structures
        {
            std::lock_guard<std::mutex> lock(mtx_i2e);
            for (const auto& [key, value] : local_i2e) {
                i2e[key].insert(value.begin(), value.end());
            }
        }
        {
            std::lock_guard<std::mutex> lock(mtx_i2t);
            for (const auto& [key, value] : local_i2t) {
                i2t[key].insert(value.begin(), value.end());
            }
        }
        {
            std::lock_guard<std::mutex> lock(mtx_i2i);
            for (const auto& [key, value] : local_i2i) {
                i2i[key].insert(value.begin(), value.end());
            }
        }
        {
            std::lock_guard<std::mutex> lock(mtx_coord);
            for (const auto& [key, value] : local_coord) {
                coord[key] = value;
            }
        }
        {
            std::lock_guard<std::mutex> lock(mtx_coord_globalIndex);
            for (const auto& [key, value] : local_coord_globalIndex) {
                coord_globalIndex[key] = value;
            }
        }
        {
            std::lock_guard<std::mutex> lock(mtx_i2T);
            for (const auto& [key, value] : local_i2T) {
                i2T[key] = value;
            }
        }
        {
            std::lock_guard<std::mutex> lock(mtx_map_IGamma);
            for (const auto& [key, value] : local_map_IGamma) {
                map_IGamma[key].insert(value.begin(), value.end());
            }
        }
    };

    // Launch threads
    std::vector<std::thread> threads;
    for (int t = 0; t < numThreads; ++t) {
        int startIdx = t * cellsPerThread;
        int endIdx = std::min(startIdx + cellsPerThread, number_elem);
        threads.emplace_back(processCells, startIdx, endIdx);
    }

    // Join threads
    for (auto& thread : threads) {
        thread.join();
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

// Multithreaded version of the function
template<class FSElement, class Function, class Material>
void markedIndicesOnInterfacesForeachSubdomain(
    FSElement& fse,
    Function const& fu,
    Material const& material,
    std::map<std::pair<int, int>, std::vector<double>>& coord,
    std::vector<std::vector<int>>& e2i,
    std::vector<std::set<int>>& e2e,
    std::vector<std::set<int>>& i2i,
    std::map<int, std::set<int>>& map_GammaGamma,
    std::map<int, std::set<int>>& map_GammaGamma_W_Nbr,
    std::map<int, std::map<int, std::set<int>>>& map_GammaNbr)
{
    typedef typename FSElement::Space ImageSpace;
    typedef typename ImageSpace::Grid Grid;

    auto const& gridView = fse.space().gridView();
    auto cbegin = fse.space().gridView().template begin<0>();
    auto const cend = fse.space().gridView().template end<0>(); //  cell end

    std::mutex mutex_e2e, mutex_i2i, mutex_map_GammaGamma, mutex_map_GammaGamma_W_Nbr, mutex_map_GammaNbr;

    size_t numThreads = std::thread::hardware_concurrency();
    size_t totalCells = std::distance(cbegin, cend);
    size_t chunkSize = (totalCells + numThreads - 1) / numThreads;

    auto worker = [&](size_t start, size_t end) {
      auto ci = cbegin;
      std::advance(ci, start);
        for (size_t idx = start; idx < end && ci != cend; ++idx, ++ci) 
        {
          auto eIndex = fse.space().indexSet().index(*ci);

          typename ImageSpace::Evaluator isfs(fse.space());
          isfs.moveTo(*ci);

          auto const& localCoordinate(isfs.shapeFunctions().interpolationNodes());
          auto dof_u = fse.space().mapper().globalIndices(*ci);
          int nrNodes = dof_u.size();

          Dune::FieldVector<double, ImageSpace::dim> zero(0.0);
          int material_var = material.value(*ci, zero);

          // Update e2e
          {
              std::lock_guard<std::mutex> lock(mutex_e2e);
              e2e[eIndex].insert(eIndex);
          }

          for (auto const& intersection : intersections(fse.space().gridView(), *ci)) {
              if (intersection.neighbor()) {
                  int eNbrIndex = fse.space().gridView().indexSet().index(intersection.outside());

                  if (material.value(*ci, zero) != material.value(intersection.outside(), zero)) {
                      int material_nbr = material.value(intersection.outside(), zero);

                      for (int i = 0; i < e2i[eIndex].size(); ++i) {
                          int nIndex_c1 = e2i[eIndex][i];
                          std::pair<int, int> pairs{nIndex_c1, material_var};

                          std::set<int> s_index;
                          {
                              std::lock_guard<std::mutex> lock(mutex_i2i);
                              s_index = i2i[nIndex_c1];
                          }

                          for (int j = 0; j < e2i[eNbrIndex].size(); ++j) {
                              int nIndex_c2 = e2i[eNbrIndex][j];
                              std::pair<int, int> pairs_nbr{nIndex_c2, material_nbr};

                              bool matched = false;
                              if (ImageSpace::dim == 2) {
                                  matched = coord[pairs][0] == coord[pairs_nbr][0] &&
                                            coord[pairs][1] == coord[pairs_nbr][1];
                              } else if (ImageSpace::dim == 3) {
                                  matched = coord[pairs][0] == coord[pairs_nbr][0] &&
                                            coord[pairs][1] == coord[pairs_nbr][1] &&
                                            coord[pairs][2] == coord[pairs_nbr][2];
                              }

                              if (matched) {
                                  s_index.insert(nIndex_c2);

                                  {
                                      std::lock_guard<std::mutex> lock(mutex_map_GammaGamma);
                                      map_GammaGamma[material_var].insert(nIndex_c1);
                                  }

                                  {
                                      std::lock_guard<std::mutex> lock(mutex_map_GammaGamma_W_Nbr);
                                      auto& gammaNbr = map_GammaGamma_W_Nbr[material_var];
                                      gammaNbr.insert(nIndex_c1);
                                      gammaNbr.insert(nIndex_c2);
                                  }

                                  {
                                      std::lock_guard<std::mutex> lock(mutex_map_GammaNbr);
                                      auto& gammaNbr = map_GammaNbr[material_var][material_nbr];
                                      gammaNbr.insert(nIndex_c2);
                                  }

                                  break;
                              }
                          }

                          {
                              std::lock_guard<std::mutex> lock(mutex_i2i);
                              i2i[nIndex_c1] = s_index;
                          }
                      }

                      {
                          std::lock_guard<std::mutex> lock(mutex_e2e);
                          e2e[eIndex].insert(eNbrIndex);
                      }
                  }
              }
          }
        }
    };

    std::vector<std::thread> threads;
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunkSize;
        size_t end = std::min(start + chunkSize, totalCells);
        threads.emplace_back(worker, start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }
}

template<class FSElement, class Function, class Material>
void markedIndicesForDirichlet(FSElement& fse,
                               Function const& fu,
                               Material const& material,
                               const std::vector<int>& arr_extra,
                               const std::vector<std::vector<int>>& cell2Indice,
                               std::set<int>& dofsDirichlet) {
    // Convert arr_extra to a set for fast lookups
    std::unordered_set<int> arr_extra_set(arr_extra.begin(), arr_extra.end());

    using ImageSpace = typename FSElement::Space;
    using Grid = typename ImageSpace::Grid;

    auto gridView = fse.space().gridView();
    auto cbegin = gridView.template begin<0>();
    auto cend = gridView.template end<0>();

    int totalCells = std::distance(cbegin, cend);
    unsigned int numThreads = std::max(1u, std::thread::hardware_concurrency());
    int cellsPerThread = (totalCells + numThreads - 1) / numThreads;

    // Mutex to synchronize updates to the shared dofsDirichlet set
    std::mutex mutex_dofsDirichlet;

    // Thread function to process a subset of cells
    auto processCells = [&](int startIdx, int endIdx) {
        std::unordered_set<int> localDofsDirichlet; // Thread-local set

        auto ci = cbegin;
        std::advance(ci, startIdx);
        for (int cellIdx = startIdx; cellIdx < endIdx && ci != cend; ++cellIdx, ++ci) {
            auto cellIndex = fse.space().indexSet().index(*ci);

            // Ensure cellIndex is within bounds of cell2Indice
            assert(cellIndex >= 0 && cellIndex < static_cast<int>(cell2Indice.size()));

            Dune::FieldVector<double, ImageSpace::dim> zero(0.0);
            int material_var = material.value(*ci, zero);

            if (arr_extra_set.find(material_var) != arr_extra_set.end()) {
                for (auto const& intersection : intersections(gridView, *ci)) {
                    if (!intersection.neighbor()) {
                        for (int index_c1 : cell2Indice[cellIndex]) {
                            localDofsDirichlet.insert(index_c1);
                        }
                    }
                }
            }
        }

        // Merge thread-local set into the global set
        std::lock_guard<std::mutex> lock(mutex_dofsDirichlet);
        dofsDirichlet.insert(localDofsDirichlet.begin(), localDofsDirichlet.end());
    };

    // Launch threads
    std::vector<std::thread> threads;
    for (unsigned int t = 0; t < numThreads; ++t) {
        int startIdx = t * cellsPerThread;
        int endIdx = std::min(startIdx + cellsPerThread, totalCells);
        threads.emplace_back(processCells, startIdx, endIdx);
    }

    // Join threads
    for (auto& thread : threads) {
        thread.join();
    }
}


#endif
