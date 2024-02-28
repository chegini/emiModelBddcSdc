#ifndef EMI_MESH_DATA_HH
#define EMI_MESH_DATA_HH

template< class FSElement, class Material>
void mesh_data_structure( FSElement& fse,  
                          Material const & material, 
                          std::vector<std::vector<int>> & e2i,
                          std::vector<std::set<int>> & i2e,
                          std::vector<std::set<int>> & e2e,
                          std::vector<std::set<int>> & i2i,
                          std::vector<std::vector<double>> & icoord,
                          std::vector<int>& i2T,
                          std::map<int, int> & map_t2l,
                          std::map<int, int> & map_sT2l,
                          std::map<int,std::set<int>> & map_II,
                          std::map<int,std::set<int>> & map_IGamma,
                          std::map<int,std::set<int>> & map_GammaGamma)
{
 	getInnerInterfaceDofsForeachSubdomain(fse,  
                                        GetGlobalCoordinate(),
                                        material, 
                                        e2i, 
                                        i2e,
                                        i2i,
                                        icoord,
                                        i2T,
                                        map_t2l,
                                        map_IGamma);

  markedIndicesOnInterfacesForeachSubdomain(fse,
                                            GetGlobalCoordinate(), 
                                            material,
                                            icoord,
                                            e2i,
                                            e2e,
                                            i2i,
                                            map_GammaGamma);



 int n_count = 0;
  {
    std::map<int, int>::iterator it;

    for (it = map_t2l.begin(); it != map_t2l.end(); it++)
    {
      map_sT2l[n_count] = it->second; 
      if(false){
	    std::cout 	<< "pre( "<<it->first    // string (key)
					<< ':'
					<< it->second << " ) -> "   // string's value 
					<< "next( "<< n_count    // string (key)
					<< ':'
					<< it->second << " ) "   // string's value 
					<< std::endl;      	
      } 
      n_count++;
    }
  }

  for ( const auto &IGamma : map_IGamma ) {
    std::vector<int> igamma(IGamma.second.begin(), IGamma.second.end());
    std::vector<int> gamma(map_GammaGamma[IGamma.first].begin(), map_GammaGamma[IGamma.first].end());

    std::set<int> I;
    std::set_difference(IGamma.second.begin(), IGamma.second.end(), 
                        map_GammaGamma[IGamma.first].begin(), map_GammaGamma[IGamma.first].end(),
                        std::inserter(I, I.end()));

    map_II[IGamma.first] = I;
    if(false) std::cout << IGamma.first << ": " << IGamma.second.size() << ", "<<map_GammaGamma[IGamma.first].size() <<", " << I.size() << std::endl;
    
  }

  if(false)
  {
		for (int i = 0; i < i2i.size(); ++i)
		{
		std::set<int> s = i2i[i];
		std::set<int>::iterator itr;

		if(s.size()>1){
		  std::cout << i << " : ";
		  for (itr = s.begin(); itr != s.end(); itr++) 
		  {
		    std::cout << *itr << " ";
		  }
		  std::cout << "\n";  
		}
		}

		std::cout << "igamma"<< ":\n";
		for ( const auto &IGamma : map_IGamma ) {
		std::cout << IGamma.first << ": ";

		std::set<int>::iterator it;
		for (it = IGamma.second.begin(); it != IGamma.second.end(); ++it) {
		   std::cout << *it << " ";
		}
		std::cout << "\n";
		}

		std::cout << "gamma"<< ":\n";
		for ( const auto &GammaGamma : map_GammaGamma ) {
		std::cout << GammaGamma.first << ": ";

		std::set<int>::iterator it;
		for (it = GammaGamma.second.begin(); it != GammaGamma.second.end(); ++it) {
		   std::cout << *it << " ";
		}
		std::cout << "\n";
		}

		std::cout << "I"<< ":\n";
		for ( const auto &II : map_II ) {
		std::cout << II.first << ": ";

		std::set<int>::iterator it;
		for (it = II.second.begin(); it != II.second.end(); ++it) {
		   std::cout << *it << " ";
		}
		std::cout << "\n";
		}


		std::cout << "t2l"<< ":\n";
		for ( const auto &t2l : map_t2l ) {
		std::cout << t2l.first << ": " << t2l.second << "\n";
		}
  }
}

#endif