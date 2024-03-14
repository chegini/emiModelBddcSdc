#ifndef EMI_MESH_DATA_HH
#define EMI_MESH_DATA_HH

#include <iostream>
#include <map>
#include "io/matlab.hh"
#include <unordered_map>
#include "mg/bddc.hh"
#include "io/matlab.hh"

#include "timestepping/semieuler.hh"

using namespace Kaskade;
using namespace Kaskade::BDDC;
using namespace std;

typedef pair<int, int> row_col;

#include <fstream>
#include <vector>

// struct vec {
//     size_t N;
//     std::vector<double> data;

//     vec(size_t N) : N(N), data(N) {};

//     vec(std::ifstream &in) {
//         in.read(reinterpret_cast<char*>(&N), sizeof(size_t));
//         data.resize(N);
//         in.read(reinterpret_cast<char*>(data.data()), N * sizeof(double));
//     }

//     void write(std::ofstream &out) {
//         out.write(reinterpret_cast<char*>(&N), sizeof(size_t));
//         out.write(reinterpret_cast<char*>(data.data()), N * sizeof(double));
//     }
// };

// struct csr_matrix {
//     size_t N, nnz;
//     std::vector<size_t> row_ptrs;
//     std::vector<size_t> col_idxs;
//     std::vector<double> values;

//     csr_matrix(size_t N, size_t nnz) : N(N), nnz(nnz), row_ptrs(N + 1, 0), col_idxs(nnz), values(nnz)
//     {}

//     csr_matrix(std::ifstream &in) {
//         in.read(reinterpret_cast<char*>(&N), sizeof(size_t));
//         in.read(reinterpret_cast<char*>(&nnz), sizeof(size_t));
//         row_ptrs.resize(N + 1);
//         col_idxs.resize(nnz);
//         values.resize(nnz);
//         in.read(reinterpret_cast<char*>(row_ptrs.data()), (N + 1) * sizeof(size_t));
//         in.read(reinterpret_cast<char*>(col_idxs.data()), nnz * sizeof(size_t));
//         in.read(reinterpret_cast<char*>(values.data()), nnz * sizeof(double));
//     }

//     void write(std::ofstream &out) {
//         out.write(reinterpret_cast<char*>(&N), sizeof(size_t));
//         out.write(reinterpret_cast<char*>(&nnz), sizeof(size_t));
//         out.write(reinterpret_cast<char*>(row_ptrs.data()), (N + 1) * sizeof(size_t));
//         out.write(reinterpret_cast<char*>(col_idxs.data()), nnz * sizeof(size_t));
//         out.write(reinterpret_cast<char*>(values.data()), nnz * sizeof(double));
//     }
// };

template <class Entry, class Index, class VEntry>
void writeToMatlabPath(NumaBCRSMatrix<Entry,Index> const& A, Dune::BlockVector<VEntry> const& b, std::string const& basename, std::string const& path, bool gen_rhs, int precision=16)
{
  std::string fname = path+"/"+basename + ".m";
  std::ofstream f(fname.c_str());
  f.precision(precision);

  // Write vector.
  f << "function [A,b] = " << basename << '\n'
  << " b = [\n";
  for (size_t i=0; i<b.N(); ++i)
    f << b[i] << std::endl;
  f << "];\n";
  
  // ------------------------------------
  // create random vecotr of size 100
  vec b_binary(b.N());
  for (size_t i = 0; i < b.N(); i++) {
      b_binary.data[i] = b[i];
  }
  // ------------------------------------

  // ------------------------------------


  // write matrix in triplet format
  f << "data = [\n";
  int n = Entry::rows;
  int m = Entry::cols;


  for (auto row = A.begin(); row!=A.end(); ++row)
  {
    for (int i=0; i<n; ++i)
    {
      auto ridx = row.index()*n + i;
      for (auto col=row.begin(); col!=row.end(); ++col)
      {
        for (int j=0; j<m; ++j)
        {  
          auto cidx = col.index()*m + j;
          f << ridx+1 << ' ' << cidx+1 << ' ' << (*col)[i][j] << std::endl;
        }
      }
    }
  }

  int nnz = 0;
  nnz = A.nonzeroes();
  // std::cout << "number of nonzero: " << nnz <<std::endl;

  csr_matrix A_binary(b.N(), 0);
  int count = 0;
  for (auto row = A.begin(); row!=A.end(); ++row)
  {
    for (int i=0; i<n; ++i)
    {
      auto ridx = row.index()*n + i;
      for (auto col=row.begin(); col!=row.end(); ++col)
      {
        for (int j=0; j<m; ++j)
        {
          auto val = (*col)[i][j];
          if (val != 0.0) {
            auto cidx = col.index()*m + j;
            A_binary.col_idxs.emplace_back(cidx);
            A_binary.values.emplace_back(val);
            count++;
          }
        }
      }
      A_binary.row_ptrs[ridx+1] = count;
    }
  }
  A_binary.nnz = count;
  f << "];\nA = sparse(data(:,1),data(:,2),data(:,3)," << n*A.N() << "," << m*A.M() << ");\n";
  // write the vector and matrix to disk
  std::string fname_vec = path+"/"+basename+"_vec.bin";
  std::string fname_A = path+"/"+basename+"_A.bin";
  std::ofstream out(fname_vec.c_str(), std::ios::binary);
  b_binary.write(out);
  out.close();

  out.open(fname_A.c_str(), std::ios::binary);
  A_binary.write(out);
  out.close();

  // read the vector and matrix from disk
  std::ifstream in(fname_vec.c_str(), std::ios::binary);
  vec b2(in);
  in.close();

  in.open(fname_A.c_str(), std::ios::binary);
  csr_matrix A2(in);
  in.close();

  // // compare the original and read vector and matrix
  // for (size_t i = 0; i < A.N(); i++) {
  //   if (b_binary.data[i] != b2.data[i]) {
  //       std::cerr << "Error: b.data[" << i << "] = " << b2.data[i] << " != " << b_binary.data[i] << std::endl;
  //   }
  // }

  // if (A_binary.N != A2.N || A_binary.nnz != A2.nnz) {
  //   std::cerr << "Error: A.N = " << A2.N << " != " << A_binary.N << " or A.nnz = " << A2.nnz << " != " << A_binary.nnz << std::endl;
  // }

  // for (size_t i = 0; i < A_binary.N + 1; i++) {
  //   if (A_binary.row_ptrs[i] != A2.row_ptrs[i]) {
  //     std::cerr << "Error: A.row_ptrs[" << i << "] = " << A2.row_ptrs[i] << " != " << A_binary.row_ptrs[i] << std::endl;
  //   }
  // }

  // for (size_t i = 0; i < A_binary.nnz; i++) {
  //   if (A_binary.col_idxs[i] != A2.col_idxs[i] || A_binary.values[i] != A2.values[i]) {
  //     std::cerr << "Error: A.col_idxs[" << i << "] = " << A2.col_idxs[i] << " != " << A_binary.col_idxs[i] << " or A.values[" << i << "] = " << A2.values[i] << " != " << A_binary.values[i] << std::endl;
  //   }
  // }
}


template <class Entry, class Index, class VEntry>
void writeToMatlabPath_final(NumaBCRSMatrix<Entry,Index> const& A, Dune::BlockVector<VEntry> const& b, std::string const& basename, std::string const& path, bool gen_rhs, int precision=16)
{
  vec b_binary(b.N());
  if(gen_rhs)
  {
    for (size_t i = 0; i < b.N(); i++) {
        b_binary.data[i] = b[i];
    }
  }

  // write matrix in triplet format
  int n = Entry::rows;
  int m = Entry::cols;

  int nnz = 0;
  nnz = A.nonzeroes();

  csr_matrix A_binary(b.N(), 0);
  int count = 0;
  for (auto row = A.begin(); row!=A.end(); ++row)
  {
    for (int i=0; i<n; ++i)
    {
      auto ridx = row.index()*n + i;
      for (auto col=row.begin(); col!=row.end(); ++col)
      {
        for (int j=0; j<m; ++j)
        {
          auto val = (*col)[i][j];
          if (val != 0.0) {
            auto cidx = col.index()*m + j;
            A_binary.col_idxs.emplace_back(cidx);
            A_binary.values.emplace_back(val);
            count++;
          }
        }
      }
      A_binary.row_ptrs[ridx+1] = count;
    }
  }
  A_binary.nnz = count;

  // write the vector and matrix to disk
  std::string fname_vec = path+"/"+basename+"_vec.bin";
  std::string fname_A = path+"/"+basename+"_A.bin";

  if(gen_rhs){
    std::ofstream out(fname_vec.c_str(), std::ios::binary);
    b_binary.write(out);
    out.close();
  }
  std::ofstream out(fname_A.c_str(), std::ios::binary);
  A_binary.write(out);
  out.close();
}

template< class FSElement, class Material>
void mesh_data_structure( FSElement& fse,  
                          Material const & material, 
                          std::vector<int> arr_extra,  
                          std::vector<std::vector<int>> & e2i,                         
                          std::vector<std::set<int>> & i2e,                            
                          std::vector<std::set<int>> & i2t,                            
                          std::vector<std::set<int>> & e2e,                             
                          std::vector<std::set<int>> & i2i,                             
                          std::map<std::pair<int, int>, std::vector<double>> & coord,  
                          std::map<int, std::vector<double>> & coord_globalIndex, 
                          std::vector<int>& i2Tag,                                      
                          std::map<int, int> & map_t2l,                                
                          std::map<int, int> & map_sT2l,                               
                          std::map<int,std::set<int>> & map_II,                         
                          std::map<int,std::set<int>> & map_IGamma,                    
                          std::map<int,std::set<int>> & map_GammaGamma,
                          std::map<int,std::set<int>> & map_IGamma_noDuplicate,
                          std::map<int,std::set<int>> & map_GammaGamma_noDuplicate,                 
                          std::map<int,std::set<int>> & map_GammaGamma_W_Nbr,  
                          std::map<int,std::set<int>> & map_GammaNbr_Nbr,  
                          std::map<int,std::set<int>> & map_GammaNbr_Nbr_noDuplicate,       
                          std::map<int,std::map<int,std::set<int>>> & map_GammaNbr,
                          std::set<int> & interface_extra_dofs)     
{
 	getInnerInterfaceDofsForeachSubdomain(fse,  
                                        GetGlobalCoordinate(),
                                        material, 
                                        e2i, 
                                        i2e,
                                        i2t,
                                        i2i,
                                        coord,
                                        coord_globalIndex,
                                        i2Tag,
                                        map_IGamma);

  markedIndicesOnInterfacesForeachSubdomain(fse,
                                            GetGlobalCoordinate(), 
                                            material,
                                            coord,
                                            e2i,
                                            e2e,
                                            i2i,
                                            map_GammaGamma,
                                            map_GammaGamma_W_Nbr,
                                            map_GammaNbr);



  std::set<int> arr_extra_set(arr_extra.begin(), arr_extra.end());

  for ( const auto &IGamma : map_IGamma ) {
    int tag = IGamma.first;
    std::vector<int> igamma(IGamma.second.begin(), IGamma.second.end());
    std::vector<int> gamma(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

    std::set<int> I;
    std::set_difference(IGamma.second.begin(), IGamma.second.end(), 
                        map_GammaGamma[tag].begin(), map_GammaGamma[tag].end(),
                        std::inserter(I, I.end()));

    map_II[IGamma.first] = I;
  
    // auto it_extra = std::find(arr_extra.begin(), arr_extra.end(), tag);
    // if (it_extra == arr_extra.end()) {  
    //   // std::cout << "~~~~~ ~~~~~ ~~~~~ " << tag << std::endl;
    //   map_GammaGamma_noDuplicate[tag] = map_GammaGamma[tag];
    // }
    // if(true) std::cout <<"tag: " <<  IGamma.first << " inner + gamma: " << IGamma.second.size() << " , gamma: "<<map_GammaGamma[IGamma.first].size() <<" , inner: " << I.size() << std::endl;
  }


  {
    std::map<int,std::set<int>>::iterator it;
    for (it=map_IGamma.begin(); it!=map_IGamma.end(); ++it){
      map_t2l[it->first] = it->second.size();
    }
  }

  // GET only neighbrs without duplication the extra cellular
  // question, I should make them be,ong to one subdomain???
  for ( const auto &Gamma : map_GammaGamma ) 
  {
    int tag = Gamma.first;
    std::set<int> difference;

    // Use std::set_difference to find the difference between set1 and set2
    std::set_difference(map_GammaGamma_W_Nbr[tag].begin(), map_GammaGamma_W_Nbr[tag].end(),
                        map_GammaGamma[tag].begin(), map_GammaGamma[tag].end(),
                        std::inserter(difference, difference.begin()));
    map_GammaNbr_Nbr[tag] = difference;
  }

  // GET only neighbrs without duplication the extra cellular
  // question, I should make them be,ong to one subdomain???
  for ( const auto &Gamma : map_IGamma_noDuplicate ) 
  {
    int tag = Gamma.first;
    std::set<int> difference;

    // Use std::set_difference to find the difference between set1 and set2
    std::set_difference(map_GammaGamma_W_Nbr[tag].begin(), map_GammaGamma_W_Nbr[tag].end(),
                        map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end(),
                        std::inserter(difference, difference.begin()));
    map_GammaNbr_Nbr_noDuplicate[tag] = difference;
  }


  // REMOVE THE INTERFACES FROM ONE OF THE EXTERNAL, 
  // IT MEANS ONLY ONE EXTRA CELLULAR SUBDOMAIN WILL HAVE THE DOFS ON THE COMMON INTERFACES BETWEEN DIFFERENT EXTRACELLULAR SUBDOMAINS

  std::cout << " remove teh extracellular interfaces " <<std::endl;
  map_IGamma_noDuplicate = map_IGamma;
  map_GammaGamma_noDuplicate = map_GammaGamma;
  for (int i = 0; i < i2t.size(); ++i)
  {
    std::set<int> s = i2t[i];
    std::set<int> s_ii = i2i[i];
    std::set<int>::iterator itr;
    std::set<int>::iterator itr_s_ii;
    std::set<int>::iterator itr_g_nbr;

    if(s.size()>1){
      int count = 0;
      for (itr = s.begin(); itr != s.end(); itr++) 
      {
        if(count==0) {
          i2Tag[i] = *itr;
          interface_extra_dofs.insert(i);
        }
        if(count>0){
          map_IGamma_noDuplicate[*itr].erase(i);
          map_GammaGamma_noDuplicate[*itr].erase(i);
          
        
          
          // for (itr_s_ii = s_ii.begin(); itr_s_ii != s_ii.end(); itr_s_ii++) 
          // {
          //   itr_g_nbr = map_GammaGamma_W_Nbr[*itr].find(*itr_s_ii);
          //   if(itr_g_nbr!=map_GammaGamma_W_Nbr[*itr].end())
          //     map_GammaGamma_W_Nbr[*itr].erase(*itr_s_ii);
          // }
        }
        count++;
      }
    }else{
      std::vector<int> v(s.begin(),s.end());
      i2Tag[i] = v[0];
    }
  }

  int n_count = 0;
  {
    std::map<int, int>::iterator it;

    for (it = map_t2l.begin(); it != map_t2l.end(); it++)
    {
      map_sT2l[n_count] = it->second; 
      if(true){
      std::cout   << "pre( "<<it->first    // tag
          << ':'
          << it->second << " ) -> "   // length of dosf for this tag 
          << "next( "<< n_count    // new tag
          << ':'
          << it->second << " ) "   // length of dosf for this tag 
          << std::endl;       
      } 
      n_count++;
    }
  }

  if(true)
  {
    std::cout << "  i2Tag "<< std::endl;
    for (int i = 0; i < i2Tag.size(); ++i)
    {
        std::cout << i << " -> " <<i2Tag[i]<< std::endl;
    }

    std::cout << "  i2t "<< std::endl;
    for (int i = 0; i < i2t.size(); ++i)
    {
        std::cout << i << " -> " <<i2t[i].size()<< std::endl;
    }

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

    std::cout << "interface_extra_dofs" << interface_extra_dofs.size()<< ":\n";
    {
    std::set<int>::iterator it;
    for (it = interface_extra_dofs.begin(); it != interface_extra_dofs.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
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

    std::cout << "igamma_noDuplicate"<< ":\n";
    for ( const auto &IGamma : map_IGamma_noDuplicate ) {
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

    std::cout << "map_GammaGamma_noDuplicate"<< ":\n";
    for ( const auto &GammaGamma : map_GammaGamma_noDuplicate ) {
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

		std::cout << "gamma with nbr"<< ":\n";
		for ( const auto &Gamma_W_Nbr : map_GammaGamma_W_Nbr ){
			std::cout << Gamma_W_Nbr.first << ": ";
			std::set<int>::iterator itr;
			for (itr = Gamma_W_Nbr.second.begin(); itr != Gamma_W_Nbr.second.end(); itr++) 
		  {
		    std::cout << *itr << " ";
		  }
		  std::cout << "\n";
		}	

		std::cout << "nbr"<< ":\n";
		for ( const auto &gammaNbr : map_GammaNbr ){
			std::cout << gammaNbr.first << ":";
			for ( const auto &nbr : gammaNbr.second ){
				std::cout << nbr.first << " (";
				std::set<int>::iterator itr;
				for (itr = nbr.second.begin(); itr != nbr.second.end(); itr++) 
			  {
			    std::cout << *itr << " ";
			  }
			  std::cout << ") ";
			}
			std::cout << "\n";
		}

    std::cout << " only Gamma nbrs(diff) " <<std::endl;
    for ( const auto &Gammanbr : map_GammaNbr_Nbr ) 
    {
      int tag = Gammanbr.first;
      std::set<int>::iterator itr;
      for (itr = Gammanbr.second.begin(); itr != Gammanbr.second.end(); itr++) 
      {
        std::cout << *itr << " ";
      }
      std::cout << "\n";  
    }

		std::cout << "t2l"<< ":\n";
		for ( const auto &t2l : map_t2l ) {
			std::cout << t2l.first << ": " << t2l.second << "\n";
		}
  }
}

void local2GlobalMapSubdomain(std::vector<int> I, std::vector<int>  gamma, 
                              std::unordered_map<int, int> &local2Global_subdomain)
{
  int index = 0;
  for (int i = 0; i < I.size(); ++i)
  {
    local2Global_subdomain.insert(pair<int, int>(index,I[i])); 
    index++;  
  }

  for (int i = 0; i < gamma.size(); ++i)
  {
    local2Global_subdomain.insert(pair<int, int>(index,gamma[i]));     
    index++;  
  }

	// for ( const auto &gamma_nbr : gamma_nbr_subdomain ) {
	// 	int tag = gamma_nbr.first;
	// 	std::set<int> nbr = gamma_nbr.second;
	// 	std::set<int>::iterator it;

  //   for (it = nbr.begin(); it != nbr.end(); it++)
  //   {
  //   	local2Global_subdomain.insert(pair<int, int>(index,*it)); 
  //     index++; 
  //   }
	// }
}

void global2LocalMapSubdomain(std::vector<int> I, std::vector<int>  gamma, 
                              std::unordered_map<int, int> &global2Local_subdomain)
{
  int index = 0;
  for (int i = 0; i < I.size(); ++i)
  {
    global2Local_subdomain.insert(pair<int, int>(I[i], index)); 
    index++;  
  }

  for (int i = 0; i < gamma.size(); ++i)
  {
    global2Local_subdomain.insert(pair<int, int>(gamma[i], index));     
    index++;  
  }

	// for ( const auto &gamma_nbr : gamma_nbr_subdomain ) {
	// 	int tag = gamma_nbr.first;
	// 	std::set<int> nbr = gamma_nbr.second;
	// 	std::set<int>::iterator it;

  //   for (it = nbr.begin(); it != nbr.end(); it++)
  //   {
  //   	global2Local_subdomain.insert(pair<int, int>(index,*it)); 
  //     index++; 
  //   }
	// }
}

void globalIndicesSubdomain(std::vector<int> I, std::vector<int>  gamma, 
                            std::vector<int> & globalIndices_subdomain)
{
  int index = 0;
  for (int i = 0; i < I.size(); ++i)
  {
    globalIndices_subdomain.push_back(I[i]);
    index++;  
  }

  for (int i = 0; i < gamma.size(); ++i)
  {
    globalIndices_subdomain.push_back(gamma[i]); 
    index++;  
  }

	// for ( const auto &gamma_nbr : gamma_nbr_subdomain ) {
	// 	int tag = gamma_nbr.first;
	// 	std::set<int> nbr = gamma_nbr.second;
	// 	std::set<int>::iterator it;
	// 	std::vector<int> nbr_vec(nbr.begin(), nbr.end());

	// 	for (int j = 0; j < nbr_vec.size(); ++j)
  //   {
  //     globalIndices_subdomain.push_back(nbr_vec[j]);
  //     index++; 
  //   }
	// }
}


void subdomain_indices( std::vector<int> sequenceOfTags, 
                        std::map<int,std::set<int>> & map_II,
                        std::map<int,std::set<int>> & map_GammaGamma_W_Nbr, 
                        std::map<int,std::unordered_map<int, int>> & local2Global,
                        std::map<int,std::unordered_map<int, int>> & global2Local,
                        std::map<int,std::vector<int>> & globalIndices)
{
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
		int tag =  sequenceOfTags[index];
		std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
		std::vector<int> Interface(map_GammaGamma_W_Nbr[tag].begin(), map_GammaGamma_W_Nbr[tag].end());

    std::unordered_map<int, int> local2Global_subIdx;
    std::unordered_map<int, int> global2Local_subIdx;
    std::vector<int> globalIndices_subIdx;

    local2GlobalMapSubdomain(Interior,Interface,local2Global_subIdx);
    global2LocalMapSubdomain(Interior,Interface,global2Local_subIdx);
    globalIndicesSubdomain(Interior,Interface,globalIndices_subIdx);
    local2Global[tag] = local2Global_subIdx;
    global2Local[tag] = global2Local_subIdx;
    globalIndices[tag] = globalIndices_subIdx;

    if(true){
      std::cout <<" size of globalIndices_subIdx " << globalIndices_subIdx.size() << std::endl; 
      for (int i = 0; i < globalIndices_subIdx.size(); ++i)
      {
        std::cout << globalIndices_subIdx[i] << " ";
      }
      std::cout << "\n";
      std::cout <<"============================================================================="  << std::endl; 
    } 
	}
}

void map_kaskade2petcs(std::vector<int> sequenceOfTags, 
                       std::map<int,std::set<int>> & map_II,
                       std::map<int,std::set<int>> & map_GammaGamma_noDuplicate,
                       std::map<int, int> & map_indices)
{
  int counter = 0;
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
   	
   	std::vector<int> I_vec(map_II[tag].begin(),map_II[tag].end());
   	std::vector<int> gamma_vec(map_GammaGamma_noDuplicate[tag].begin(),map_GammaGamma_noDuplicate[tag].end());

    for (int i = 0; i < I_vec.size(); ++i)
    {

      std::pair<int,int> pairs;
      pairs.first = I_vec[i];
      pairs.second = tag;

      map_indices[I_vec[i]] = counter;
      std::cout << I_vec[i] << ": "  << counter <<std::endl;
      counter++;
    }

    for (int i = 0; i < gamma_vec.size(); ++i)
    {
      std::pair<int,int> pairs;
      pairs.first = gamma_vec[i];
      pairs.second = tag;

      map_indices[gamma_vec[i]] = counter;
      std::cout << gamma_vec[i] << ": "  << counter <<std::endl;
      counter++;
    }
  }
}

void removeInnerIndices_i2i(std::vector<std::set<int>> & i2i)
{
  for (int i = 0; i < i2i.size(); ++i)
  {
    std::set<int> s = i2i[i];
    // if(s.size()==1) {
    //   i2i.erase(i2i.begin()+i);
    // }
  }
}

void compute_sharedDofsKaskade( std::vector<int> sequenceOfTags, 
                                std::map<std::pair<int, int>, int> map_indices, 
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma, 
                                std::map<int,std::map<int,std::set<int>>> map_GammaNbr, 
                                bool write_to_file,
                                std::string matlab_dir,  
                                std::vector<std::vector<LocalDof>> & sharedDofsKaskade)
{
  // construct the sequance of the submatrices in BDDC in kaskade
  int n_subdomains = map_II.size();
  std::map<int,std::vector<int>> sequanceOfsubdomainsKaskade;
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
    std::vector<int> I(map_II[tag].begin(),map_II[tag].end());
    std::vector<int> gamma(map_GammaGamma[tag].begin(),map_GammaGamma[tag].end());
    std::map<int,std::set<int>> gamma_nbr_subdomain(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
    int size = I.size()+gamma.size();
    for ( const auto &gamma_nbr : gamma_nbr_subdomain ) 
    {
     size += gamma_nbr.second.size(); 
    }


    std::vector<int> vec_kaskadeIndex(size);
    int count = 0;
    for (int i = 0; i < I.size(); ++i){
      std::pair<int,int> pairs;
      pairs.first = I[i];
      pairs.second = tag;

      vec_kaskadeIndex[count] = map_indices[pairs]; 
      count++;     
    }

    for (int i = 0; i < gamma.size(); ++i){
      std::pair<int,int> pairs;
      pairs.first = gamma[i];
      pairs.second = tag;

      vec_kaskadeIndex[count] = map_indices[pairs];    
      count++;  
    }
    
    for ( const auto &gamma_nbr : gamma_nbr_subdomain ) 
    {
      int tag_nbr = gamma_nbr.first;
      std::set<int> nbr = gamma_nbr.second;
      std::set<int>::iterator it;

      for (it = nbr.begin(); it != nbr.end(); it++)
      {
        std::pair<int,int> pairs;
        pairs.first = *it;
        pairs.second = tag_nbr;

        vec_kaskadeIndex[count] = map_indices[pairs]; 
        count++; 
      }
    }

    sequanceOfsubdomainsKaskade[tag] = vec_kaskadeIndex;
  }
   

  typedef std::tuple<int,int,int> i3tuple;
  std::map<int, std::vector<i3tuple>> MapSharedDofsKaskadeTuple;
  for ( const auto &II : map_II ) 
  {
    int tag =  II.first;
    std::vector<int> tmp = sequanceOfsubdomainsKaskade[tag];

    for (int i = 0; i < tmp.size(); ++i)
    {
      auto it = MapSharedDofsKaskadeTuple.find(tmp[i]);
      if (it != MapSharedDofsKaskadeTuple.end()) {  
        std::vector<i3tuple>& values = it->second;
        values.push_back(i3tuple(tag,i, tmp[i]));
      }else{
        MapSharedDofsKaskadeTuple[tmp[i]] = {i3tuple(tag,i, tmp[i])};
      }
    }
  }

  std::vector<std::vector<LocalDof>> sharedDofsKaskadeAll;
  {
    double precision = 16;
    std::string fname = matlab_dir+"/sharedDofsKaskade.txt";
    std::ofstream f(fname.c_str());
    f.precision(precision);

    for (const auto& entry : MapSharedDofsKaskadeTuple) 
    {
      int key = entry.first;
      const std::vector<i3tuple>& values = entry.second;
      {
        std::vector<LocalDof> tmp;
        for (const auto& value : values) {
          tmp.push_back({std::get<0>(value),std::get<1>(value)});
        }
        sharedDofsKaskadeAll.push_back(tmp);
      }

      if(values.size()>1){
        std::vector<LocalDof> tmp;
        for (const auto& value : values) {
          tmp.push_back({std::get<0>(value),std::get<1>(value)});
        }
        sharedDofsKaskade.push_back(tmp);
      }

      if(values.size()>1 and write_to_file){
        f << key << "-> ";
        std::cout << key << "-> ";
        for (const auto& value : values) {
          f <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
         std::cout <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
        }
        f << "\n";
        std::cout <<"\n";
      }
    }
  } 
}

// void compute_sharedDofsKaskade_moreExtraCells( std::vector<int> sequenceOfTags, 
//                                 std::map<std::pair<int, int>, int> map_indices, 
//                                 std::map<int,std::set<int>> map_II,
//                                 std::map<int,std::set<int>> map_GammaGamma, 
//                                 std::map<int,std::map<int,std::set<int>>> map_GammaNbr, 
//                                 bool write_to_file,
//                                 std::string matlab_dir,  
//                                 std::vector<std::vector<LocalDof>> & sharedDofsKaskade)
// {
//   // construct the sequance of the submatrices in BDDC in kaskade
//   int n_subdomains = map_II.size();
//   std::map<int,std::vector<int>> sequanceOfsubdomainsKaskade;
//   for (int index = 0; index < sequenceOfTags.size(); ++index)
//   { 
//     int tag =  sequenceOfTags[index];
//     std::vector<int> I(map_II[tag].begin(),map_II[tag].end());
//     std::vector<int> gamma(map_GammaGamma[tag].begin(),map_GammaGamma[tag].end());
//     std::map<int,std::set<int>> gamma_nbr_subdomain(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
//     int size = I.size()+gamma.size();
//     for ( const auto &gamma_nbr : gamma_nbr_subdomain ) 
//     {
//      size += gamma_nbr.second.size(); 
//     }


//     std::vector<int> vec_kaskadeIndex(size);
//     int count = 0;
//     for (int i = 0; i < I.size(); ++i)
//     {
//       std::pair<int,int> pairs;
//       pairs.first = I[i];
//       pairs.second = tag;

//       vec_kaskadeIndex[count] = map_indices[pairs]; 
//       count++;     
//     }

//     for (int i = 0; i < gamma.size(); ++i)
//     {
//       std::pair<int,int> pairs;
//       pairs.first = gamma[i];
//       pairs.second = tag;

//       vec_kaskadeIndex[count] = map_indices[pairs];    
//       count++;  
//     }
    
//     for ( const auto &gamma_nbr : gamma_nbr_subdomain ) 
//     {
//       int tag_nbr = gamma_nbr.first;
//       std::set<int> nbr = gamma_nbr.second;
//       std::set<int>::iterator it;

//       for (it = nbr.begin(); it != nbr.end(); it++)
//       {
//         std::pair<int,int> pairs;
//         pairs.first = *it;
//         pairs.second = tag_nbr;

//         vec_kaskadeIndex[count] = map_indices[pairs]; 
//         count++; 
//       }
//     }

//     sequanceOfsubdomainsKaskade[tag] = vec_kaskadeIndex;
//   }
   

//   typedef std::tuple<int,int,int> i3tuple;
//   std::map<int, std::vector<i3tuple>> MapSharedDofsKaskadeTuple;
//   for ( const auto &II : map_II ) 
//   {
//     int tag =  II.first;
//     std::vector<int> tmp = sequanceOfsubdomainsKaskade[tag];

//     for (int i = 0; i < tmp.size(); ++i)
//     {
//       auto it = MapSharedDofsKaskadeTuple.find(tmp[i]);
//       if (it != MapSharedDofsKaskadeTuple.end()) {  
//         std::vector<i3tuple>& values = it->second;
//         values.push_back(i3tuple(tag,i, tmp[i]));
//       }else{
//         MapSharedDofsKaskadeTuple[tmp[i]] = {i3tuple(tag,i, tmp[i])};
//       }
//     }
//   }

//   std::vector<std::vector<LocalDof>> sharedDofsKaskadeAll;
//   {
//     double precision = 16;
//     std::string fname = matlab_dir+"/sharedDofsKaskade.txt";
//     std::ofstream f(fname.c_str());
//     f.precision(precision);

//     for (const auto& entry : MapSharedDofsKaskadeTuple) 
//     {
//       int key = entry.first;
//       const std::vector<i3tuple>& values = entry.second;
//       {
//         std::vector<LocalDof> tmp;
//         for (const auto& value : values) {
//           tmp.push_back({std::get<0>(value),std::get<1>(value)});
//         }
//         sharedDofsKaskadeAll.push_back(tmp);
//       }

//       if(values.size()>1){
//         std::vector<LocalDof> tmp;
//         for (const auto& value : values) {
//           tmp.push_back({std::get<0>(value),std::get<1>(value)});
//         }
//         sharedDofsKaskade.push_back(tmp);
//       }

//       if(values.size()>1 and write_to_file){
//         f << key << "-> ";
//         std::cout << key << "-> ";
//         for (const auto& value : values) {
//           f <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
//          std::cout <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
//         }
//         f << "\n";
//         std::cout <<"\n";
//       }
//     }
//   } 
// }

void compute_sharedDofsKaskade_moreExtraCells( std::vector<int> sequenceOfTags, 
                                std::map<int, int> map_indices, 
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma, 
                                std::map<int,std::set<int>> map_GammaGamma_W_Nbr,
                                std::map<int,std::set<int>> map_GammaNbr_Nbr,
                                bool write_to_file,
                                std::string matlab_dir,  
                                std::vector<std::vector<LocalDof>> & sharedDofsKaskade)
{
  // construct the sequance of the submatrices in BDDC in kaskade
  std::map<int,std::vector<int>> sequanceOfsubdomainsKaskade;
  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];

    std::vector<int> I(map_II[tag].begin(),map_II[tag].end());
    for (int i = 0; i < I.size(); ++i)
    {
      sequanceOfsubdomainsKaskade[tag].push_back(map_indices[I[i]]);
    }
    std::vector<int> gamma(map_GammaGamma[tag].begin(),map_GammaGamma[tag].end());;
    for (int i = 0; i < gamma.size(); ++i)
    {
      sequanceOfsubdomainsKaskade[tag].push_back(map_indices[gamma[i]]);
    }

    std::vector<int> gammanbr(map_GammaNbr_Nbr[tag].begin(),map_GammaNbr_Nbr[tag].end());;
    for (int i = 0; i < gammanbr.size(); ++i)
    {
      sequanceOfsubdomainsKaskade[tag].push_back(map_indices[gammanbr[i]]);
    }
  }

  typedef std::tuple<int,int,int> i3tuple;
  std::map<int, std::vector<i3tuple>> MapSharedDofsKaskadeTuple;

  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::vector<int> tmp = sequanceOfsubdomainsKaskade[tag];

    for (int i = 0; i < tmp.size(); ++i)
    {
      auto it = MapSharedDofsKaskadeTuple.find(tmp[i]);
      if (it != MapSharedDofsKaskadeTuple.end()) {  
        std::vector<i3tuple>& values = it->second;
        values.push_back(i3tuple(subIdx,i, tmp[i]));
      }else{
        MapSharedDofsKaskadeTuple[tmp[i]] = {i3tuple(subIdx,i, tmp[i])};
      }
    }
  }

  std::vector<std::vector<LocalDof>> sharedDofsKaskadeAll;
  {
    double precision = 16;
    std::string fname = matlab_dir+"/sharedDofsKaskade.txt";
    std::ofstream f(fname.c_str());
    f.precision(precision);

    for (const auto& entry : MapSharedDofsKaskadeTuple) {
      int key = entry.first;
      const std::vector<i3tuple>& values = entry.second;
      {
        std::vector<LocalDof> tmp;
        for (const auto& value : values) {
          tmp.push_back({std::get<0>(value),std::get<1>(value)});
        }
        sharedDofsKaskadeAll.push_back(tmp);
      }

      if(values.size()>1){
        std::vector<LocalDof> tmp;
        for (const auto& value : values) {
          tmp.push_back({std::get<0>(value),std::get<1>(value)});
        }
        sharedDofsKaskade.push_back(tmp);
      }
      if(values.size()>1 and write_to_file){
        f << key << "-> ";
        for (const auto& value : values) {
          f <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
        }
        f << "\n";
      }
    }
  } 
  // // ------------------------------------------------------------------------------------------

  // for (int i = 0; i < sharedDofsKaskade.size(); ++i)
  // {
  //   for (const auto& value : sharedDofsKaskade[i]) {
  //     std::cout << "("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
  //   }
  //   std::cout << "\n";
  // }
}

template<class FSElement, class Material>
void write_Dirichlet_and_coordinates( FSElement& fse, 
                                      Material const & material, 
                                      std::vector<std::vector<int>> e2i,
                                      std::map<int, int> map_indices,
                                      std::map<std::pair<int, int>, std::vector<double>> & coord,
                                      std::map<int, std::vector<double>> & coord_globalIndex,
                                      int dof_size,
                                      int dim,
                                      bool write_to_file,
                                      std::string matlab_dir)
{


  std::set<int> dofsDirichlet;
  markedIndicesForDirichlet(fse,
                            GetGlobalCoordinate(), 
                            material,
                            e2i,
                            dofsDirichlet);


  std::cout << "DrichletNodes"<<std::endl;
  if(write_to_file)
  {
    double precision = 16;
    std::string fname = matlab_dir+"/DrichletNodes.m";
    std::ofstream f(fname.c_str());
    f.precision(precision);
    std::set<int>::iterator it;
    for (it = dofsDirichlet.begin(); it != dofsDirichlet.end(); ++it) {
      f << map_indices[*it] << " ";
    }
    f << "\n";
  }
  
  if(write_to_file)
  {
    std::vector<std::vector<double>> indexCoordinates_petsc(dof_size);
    double precision = 16;
    std::string fname = matlab_dir+"/coordinates.txt";
    std::ofstream f(fname.c_str());
    f.precision(precision);
    for (int i = 0; i < coord_globalIndex.size(); ++i)
      indexCoordinates_petsc[map_indices[i]] = coord_globalIndex[i];

    for (int j = 0; j < indexCoordinates_petsc.size(); ++j){
      if(dim==2) 
        f << j << ":"<< indexCoordinates_petsc[j][0] << " "<< indexCoordinates_petsc[j][1] << " "<< 0.0 << "\n";
      if(dim==3) 
        f << j << ":"<<indexCoordinates_petsc[j][0] << " "<< indexCoordinates_petsc[j][1] << " "<< indexCoordinates_petsc[j][2] << "\n";
    }
  }

}

void computed_sequenceOfTags(std::map<int, int> map_t2l, std::map<int,std::map<int,std::set<int>>> map_GammaNbr, std::vector<int> & sequenceOfTags, std::map<int,int> & startingIndexOfTag, std::map<int,int> & map_nT2oT, std::map<int,std::vector<int>> & sequenceOfsubdomains)
{
  int start = 0;
  int index = 0;
  int nTag = 0;
  for ( const auto &t2l : map_t2l ) 
  {
    int tag =  t2l.first;
    sequenceOfTags[index] = tag;
    startingIndexOfTag[tag] = start;
    map_nT2oT[nTag] = tag;

    start+=map_t2l[tag];
    index++;
    nTag++;
  }

  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
    std::vector<int> sequenceOfsubdomains_subIndex;
    std::cout << "tag: " << tag << ": ";
    for ( const auto &gamma_nbr : nbrs ) {
      std::cout << gamma_nbr.first << " ";
      sequenceOfsubdomains_subIndex.push_back(gamma_nbr.first);
    }
    std::cout << "\n";
    sequenceOfsubdomains[tag] = sequenceOfsubdomains_subIndex;
  } 

}

template<class Vector>
void petsc_structure_rhs( std::vector<int> sequenceOfTags, 
                          std::map<int,int> startingIndexOfTag,
                          std::map<int,int> map_indices,
                          std::map<int,std::set<int>> map_II,
                          std::map<int,std::set<int>> map_GammaGamma_noDuplicate, 
                          Vector b_,
                          Vector &bs_)
{
  bs_ *= 0;
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];

    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
    int x = startingIndexOfTag[tag];
    for (int i = 0; i < Interior.size(); ++i)
    {
      bs_[map_indices[i]] = b_[Interior[i]];
      // bs_[i+x] = b_[Interior[i]];
    }

    for (int i = 0; i < interface.size(); ++i)
    {
      bs_[map_indices[i]] = b_[interface[i]];
      // bs_[i+x+Interior.size()] = b_[interface[i]];
    }
  }
}


template<class Vector>
void petsc_structure_rhs_subdomain_petsc( std::vector<int> sequenceOfTags, 
                                std::map<int,int> startingIndexOfTag,
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                                std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate,
                                Vector b_,
                                std::map<int, int> map_indices,
                                std::vector<std::vector<LocalDof>> sharedDofsAll,
                                std::vector<Vector> &Fs)
{

  std::cout << " RHS based on petsc Sructure" <<std::endl;
  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
    std::vector<int> gammaNbr(map_GammaNbr_Nbr_noDuplicate[tag].begin(), map_GammaNbr_Nbr_noDuplicate[tag].end());

    Vector Fs_subIdx(b_.size()); 
    double coef = 2;

    for (int indx = 0; indx < Interior.size(); ++indx)
    {
      int pos = Interior[indx];
      int pos_map = map_indices[pos];
      Fs_subIdx[pos_map] = b_[pos];
    }

    for (int indx = 0; indx < interface.size(); ++indx)
    {
      int pos = interface[indx];
      int pos_map = map_indices[pos];
      // Fs_subIdx[pos_map] = coef*b_[pos];
      int alpha  = sharedDofsAll[pos].size();
      Fs_subIdx[pos_map] = (1.0/alpha)*b_[pos];
    }

    for (int indx = 0; indx < gammaNbr.size(); ++indx)
    {
      int pos = gammaNbr[indx];
      int pos_map = map_indices[pos];
      // Fs_subIdx[pos_map] = coef*b_[pos];
      int alpha  = sharedDofsAll[pos].size();
      Fs_subIdx[pos_map] = (1.0/alpha)*b_[pos];
    }

    Fs.push_back(Fs_subIdx);
  }
}

// ------------------------------------------------------------------------------------
// construct submatrices
// ------------------------------------------------------------------------------------
template<class Matrix_, class Matrix>
void insertMatrixBlock(Matrix_ A_block, double coef, int x, int y, Matrix &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      double val = *ca;
      // std::cout << "k+x:  "<< k+x << " l+y:  "<< l+y << std::endl;
      As_[k+x][l+y] = val*coef; 
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlock(Matrix A_block, int x, int y, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      As_[k+x][l+y] = *ca; 
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlockNbr(Matrix A_block, double coef, int x, int y, Matrix_ &As_)
{  
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      //if(*ca!=0.0)
      {
        As_[k+x][l+y] = coef*(*ca); 
        As_[l+y][k+x] = coef*(*ca); 
      }
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlockNbr(Matrix A_block, int x, int y, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      As_[k+x][l+y] = *ca; 
      As_[l+y][k+x] = *ca; 
    }
  }
}

template<class Matrix>
void exctract_petsc_mass_blocks(std::vector<int> sequenceOfTags, 
                                std::map<int,int> startingIndexOfTag,
                                std::map<int,std::vector<int>> sequenceOfsubdomains,
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma,
                                std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                                Matrix A_,
                                Matrix M_,
                                std::map<int,Matrix> massSubmatrices_subIndx,
                                int subIdx,
                                Matrix &Ms)
{
  int tag = sequenceOfTags[subIdx];
  std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
  std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
  std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());

  int start = startingIndexOfTag[tag];
  std::cout << "tag: "<< tag << " start: " << start << std::endl;

  double coef = 0.5;
  auto M_GAMMAGAMMA_block = M_(Interface,Interface);
  insertMatrixBlock(M_GAMMAGAMMA_block, coef, start+Interior.size(), start+Interior.size(), Ms);


  std::vector<int> sequanceOfsubdomains_subIdx(sequenceOfsubdomains[tag].begin(), sequenceOfsubdomains[tag].end());

  for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
  {
    int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
    std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
    std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());
 
    int x = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();
    int y = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();

    Matrix massSubmatrices_nbr = massSubmatrices_subIndx[nbr_Indx];
    auto K_GAMMAGAMMA_nbr_nbr_block = massSubmatrices_nbr(Interface_nbr,Interface_nbr);
    insertMatrixBlock(K_GAMMAGAMMA_nbr_nbr_block, coef, x, y, Ms);
  }

  // cross blocks on column for subIdx 
  for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
  {
    int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
    std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
    std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());

    int x = startingIndexOfTag[tag] + Interior.size();
    int y = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();

    auto K_GAMMAGAMMA_nbr_block = A_(Interface,Interface_nbr);
    insertMatrixBlockNbr(K_GAMMAGAMMA_nbr_block, coef, x, y, Ms);
  }
}

template<class Matrix, class Matrix_>
void exctract_petsc_stiffness_blocks( std::vector<int> sequenceOfTags, 
                                      std::map<int,int> startingIndexOfTag,
                                      std::map<int,std::set<int>> map_II,
                                      std::map<int,std::set<int>> map_GammaGamma,
                                      Matrix A_,
                                      Matrix K_,
                                      int subIdx,
                                      Matrix_ &Ks_)
{
  int tag = sequenceOfTags[subIdx];
  std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
  std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

  int start = startingIndexOfTag[tag];
  auto A_II_block = A_(Interior,Interior);
  insertMatrixBlock(A_II_block, start, start, Ks_);

  auto A_IGAMMA_block = A_(Interior,Interface);
  insertMatrixBlock(A_IGAMMA_block, start, start+Interior.size(), Ks_);

  auto A_GAMMAI_block = A_(Interface,Interior);
  insertMatrixBlock(A_GAMMAI_block, start+Interior.size(), start, Ks_);

  auto K_GAMMAGAMMA_block = K_(Interface,Interface);
  insertMatrixBlock(K_GAMMAGAMMA_block, start+Interior.size(), start+Interior.size(), Ks_);
}

// -------------------------------------------------------------------------------------
// reorder matrix A with respect to petsc structure
// -------------------------------------------------------------------------------------
template<class Matrix>
void petsc_structure_Matrix(std::vector<int> arr_extra, 
                            std::vector<int> sequenceOfTags, 
                            std::map<int,int> startingIndexOfTag,
                            std::map<int,std::set<int>> map_II,
                            std::map<int,std::set<int>> map_GammaGamma,
                            Matrix A_,
                            Matrix &As_)
{ 
  for (int subIdx = 0; subIdx < map_II.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

    int x = startingIndexOfTag[tag];
    int y = startingIndexOfTag[tag];

    auto A_Interior_block = A_(Interior,Interior);
    insertMatrixBlock(A_Interior_block, x, y, As_);

    auto A_Interior_Interface_block = A_(Interior,interface);
    insertMatrixBlock(A_Interior_Interface_block, x, y+Interior.size(), As_);

    auto A_Interface_Interior_block = A_(interface,Interior);
    insertMatrixBlock(A_Interface_Interior_block, x+Interior.size(), y, As_);
   
    auto A_Interface_Interface_block = A_(interface,interface);
    insertMatrixBlock(A_Interface_Interface_block, x+Interior.size(), y+Interior.size(), As_);

    for (int nbr = 0; nbr < sequenceOfTags.size(); ++nbr)
    {
      int tag_nbr = sequenceOfTags[nbr];

      bool extra_interface = true;
      if( (std::find(arr_extra.begin(), arr_extra.end(), tag_nbr) != arr_extra.end() and
         std::find(arr_extra.begin(), arr_extra.end(), tag) != arr_extra.end()) ){
        std::cout << "extra_interface: "<<"(" << tag << ","<<tag_nbr << ")" << std::endl;
        extra_interface = false;
      } 
        

      if(extra_interface)
      {
        std::cout<< "(" << tag << ","<<tag_nbr << ")" << "inner: "<< map_II[tag_nbr].size() << " gamma: " << map_GammaGamma[tag_nbr].size()<<std::endl;

        std::vector<int> Interior_nbr(map_II[tag_nbr].begin(), map_II[tag_nbr].end());
        std::vector<int> interface_nbr(map_GammaGamma[tag_nbr].begin(), map_GammaGamma[tag_nbr].end());

        x = startingIndexOfTag[tag] + Interior.size();
        y = startingIndexOfTag[tag_nbr] + Interior_nbr.size();
        auto A_Interface_Interface_nbr_block = A_(interface,interface_nbr);
        insertMatrixBlockNbr(A_Interface_Interface_nbr_block, x, y, As_); 
      }
    }
  }
}

template<class Grid, class Functional, class VariableSet, class Spaces, class Matrix, class Vector>
typename VariableSet::VariableSet  construct_submatrices_petsc( std::vector<int> arr_extra_set, 
                                                          std::map<int,int> map_nT2oT,
                                                          GridManager<Grid>& gridManager,
                                                          Functional& F,
                                                          VariableSet const& variableSet, 
                                                          Spaces const& spaces,
                                                          Grid const& grid, 
                                                          typename VariableSet::VariableSet u,
                                                          double  dt,
                                                          std::vector<int> sequenceOfTags, 
                                                          std::map<int,int> startingIndexOfTag,
                                                          std::map<int,std::set<int>> map_II,
                                                          std::map<int,std::set<int>> map_GammaGamma,
                                                          std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                                                          std::map<int,std::vector<int>> sequenceOfsubdomains,
                                                          std::map<std::pair<int,int>, int> map_indices,
                                                          std::vector<std::set<int>> i2t, 
                                                          Matrix A_,
                                                          Matrix K_,
                                                          Matrix M_,
                                                          Vector rhs_petsc_test,
                                                          int nDofs,
                                                          int assemblyThreads,
                                                          bool write_to_file,
                                                          std::string matlab_dir,
                                                          std::vector<Vector> &Fs_petcs,
                                                          std::vector<Vector> &weights,
                                                          std::vector<Matrix> &subMatrices,
                                                          std::vector<Matrix> &subMatrices_M,
                                                          std::vector<Matrix> &subMatrices_K)
{
  // ------------------------------------------------------------------------------------ 
  // construct sparsity patterns
  // ------------------------------------------------------------------------------------ 
  NumaCRSPatternCreator<> creator(nDofs,nDofs,false);
  int counter_elements = 0;
  int count_row = 0;
  for (int r=0; r<A_.N(); ++r)
  {
    auto row  = A_[r];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const c = ca.index();
      std::vector<int> tags_r(i2t[r].begin(), i2t[r].end());
      std::vector<int> tags_c(i2t[c].begin(), i2t[c].end());

      if(tags_r.size()==2 and tags_c.size()==2){
        std::pair<int,int> pairs_r1;
        pairs_r1.first = r;
        pairs_r1.second = tags_r[0];

        std::pair<int,int> pairs_c1;
        pairs_c1.first = c;
        pairs_c1.second = tags_c[0];

        int row_indx = map_indices[pairs_r1];
        int col_indx = map_indices[pairs_c1];

        if(r==2 and c ==9){
          std::cout <<"(tags_r.size()==2 and tags_c.size()==2)" <<"\n";
        }
        creator.addElement(row_indx,col_indx); 
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        // std::cout << "(" <<r << ", " << c << "), (" <<row_indx << ", " << col_indx << ") = "<< *ca <<  "\n";


        std::pair<int,int> pairs_r2;
        pairs_r2.first = r;
        pairs_r2.second = tags_r[1];

        std::pair<int,int> pairs_c2;
        pairs_c2.first = c;
        pairs_c2.second = tags_c[1];

        row_indx = map_indices[pairs_r2];
        col_indx = map_indices[pairs_c2];

        creator.addElement(row_indx,col_indx); 
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        // std::cout << "(" <<r << ", " << c << "), (" <<row_indx << ", " << col_indx << ") = "<< *ca <<  "\n"; 

        std::pair<int,int> pairs_r3;
        pairs_r3.first = r;
        pairs_r3.second = tags_r[0];

        std::pair<int,int> pairs_c3;
        pairs_c3.first = c;
        pairs_c3.second = tags_c[1];

        row_indx = map_indices[pairs_r3];
        col_indx = map_indices[pairs_c3];
        creator.addElement(row_indx,col_indx); 
        creator.addElement(col_indx,row_indx);

        std::pair<int,int> pairs_r4;
        pairs_r4.first = r;
        pairs_r4.second = tags_r[1];

        std::pair<int,int> pairs_c4;
        pairs_c4.first = c;
        pairs_c4.second = tags_c[0];

        row_indx = map_indices[pairs_r4];
        col_indx = map_indices[pairs_c4];
        creator.addElement(row_indx,col_indx); 
        creator.addElement(col_indx,row_indx);
        if(row_indx==5 and col_indx ==8){
          std::cout <<"(tags_r.size()==2 and tags_c.size()==2)" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }

      }else if(tags_r.size()==1 and tags_c.size()==1){ // other subdomains either inner or myocytes
        std::pair<int,int> pairs_r1;
        pairs_r1.first = r;
        pairs_r1.second = tags_r[0];

        std::pair<int,int> pairs_c1;
        pairs_c1.first = c;
        pairs_c1.second = tags_c[0];

        int row_indx = map_indices[pairs_r1];
        int col_indx = map_indices[pairs_c1];

        // if(r==2 and c ==9){
        //   std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        // }

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";

        std::pair<int,int> pairs_r2;
        pairs_r2.first = r;
        pairs_r2.second = tags_r[0];

        std::pair<int,int> pairs_c2;
        pairs_c2.first = c;
        pairs_c2.second = tags_c[0];

        row_indx = map_indices[pairs_r2];
        col_indx = map_indices[pairs_c2];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }

        //std::cout <<row_indx << ", " << col_indx <<"\n";
        // --------------------------------------------------------------------------------------------------------

        pairs_r1.first = r;
        pairs_r1.second = tags_c[0];

        pairs_c1.first = c;
        pairs_c1.second = tags_r[0];

        row_indx = map_indices[pairs_r1];
        col_indx = map_indices[pairs_c1];

        // if(r==2 and c ==9){
        //   std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        // }


        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }

        pairs_r2.first = r;
        pairs_r2.second = tags_c[0];

        pairs_c2.first = c;
        pairs_c2.second = tags_r[0];

        row_indx = map_indices[pairs_r2];
        col_indx = map_indices[pairs_c2];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }

        // --------------------------------------------------------------------------------------------------------

        pairs_r1.first = r;
        pairs_r1.second = tags_c[0];

        pairs_c1.first = c;
        pairs_c1.second = tags_r[0];

        row_indx = map_indices[pairs_r1];
        col_indx = map_indices[pairs_c1];

        // if(r==2 and c ==9){
        //   std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        // }

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }
        pairs_r2.first = r;
        pairs_r2.second = tags_c[0];

        pairs_c2.first = c;
        pairs_c2.second = tags_r[0];

        row_indx = map_indices[pairs_r2];
        col_indx = map_indices[pairs_c2];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);

        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()==1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }
        // --------------------------------------------------------------------------------------------------------
        // std::cout << "(" <<r << ", " << c << "), (" <<row_indx << ", " << col_indx << ") = "<< *ca <<  "\n";
      }else if(tags_c.size()!=1 and tags_r.size()==1){ // row is inner
        std::pair<int,int> pairs_r1;
        pairs_r1.first = r;
        pairs_r1.second = tags_c[0];

        std::pair<int,int> pairs_c1;
        pairs_c1.first = c;
        pairs_c1.second = tags_c[0];

        // if(r==2 and c ==9){
        //   std::cout <<"tags_c.size()!=1 and tags_r.size()==1" <<"\n";
        // }

        int row_indx = map_indices[pairs_r1];
        int col_indx = map_indices[pairs_c1];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_c.size()!=1 and tags_r.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }

        std::pair<int,int> pairs_r2;
        pairs_r2.first = r;
        pairs_r2.second = tags_c[1];

        std::pair<int,int> pairs_c2;
        pairs_c2.first = c;
        pairs_c2.second = tags_c[1];

        row_indx = map_indices[pairs_r2];
        col_indx = map_indices[pairs_c2];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx); 
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_c.size()!=1 and tags_r.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }


        std::pair<int,int> pairs_r3;
        pairs_r3.first = r;
        pairs_r3.second = tags_r[0];

        std::pair<int,int> pairs_c3;
        pairs_c3.first = c;
        pairs_c3.second = tags_c[0];

        row_indx = map_indices[pairs_r3];
        col_indx = map_indices[pairs_c3];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx); 


        std::pair<int,int> pairs_r4;
        pairs_r4.first = r;
        pairs_r4.second = tags_r[0];

        std::pair<int,int> pairs_c4;
        pairs_c4.first = c;
        pairs_c4.second = tags_c[1];

        row_indx = map_indices[pairs_r4];
        col_indx = map_indices[pairs_c4];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx); 

        //std::cout <<row_indx << ", " << col_indx <<"\n";
        // std::cout << "(" <<r << ", " << c << "), (" <<row_indx << ", " << col_indx << ") = "<< *ca <<  "\n";
      }else if(tags_r.size()!=1 and tags_c.size()==1){ // row is interface
        std::pair<int,int> pairs_r1;
        pairs_r1.first = r;
        pairs_r1.second = tags_r[0];

        std::pair<int,int> pairs_c1;
        pairs_c1.first = c;
        pairs_c1.second = tags_r[0];

        int row_indx = map_indices[pairs_r1];
        int col_indx = map_indices[pairs_c1];

        // if(r==2 and c ==9){
        //   std::cout <<"tags_r.size()!=1 and tags_c.size()==1" <<"\n";
        // }

        creator.addElement(row_indx,col_indx); 
        creator.addElement(col_indx,row_indx);
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()!=1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }
        std::pair<int,int> pairs_r2;
        pairs_r2.first = r;
        pairs_r2.second = tags_r[1];

        std::pair<int,int> pairs_c2;
        pairs_c2.first = c;
        pairs_c2.second = tags_r[1];

        row_indx = map_indices[pairs_r2];
        col_indx = map_indices[pairs_c2];

        creator.addElement(row_indx,col_indx);
        creator.addElement(col_indx,row_indx); 
        //std::cout <<row_indx << ", " << col_indx <<"\n";
        // std::cout << "(" <<r << ", " << c << "), (" <<row_indx << ", " << col_indx << ") = "<< *ca <<  "\n";

        if(row_indx==5 and col_indx ==8){
          std::cout <<"tags_r.size()!=1 and tags_c.size()==1" << "r= "<< r << " ->  "<< row_indx << " c=  "<< c << " ->  "<< col_indx <<"\n";
        }
      }
    }
  }

  // ------------------------------------------------------------------------------------ 
  // construct submatrices
  // ------------------------------------------------------------------------------------ 
  std::vector<Matrix> Ms;// for debugging
  std::vector<Matrix> Ks;// for debugging
  Matrix M_sum(creator);// for debugging
  Matrix K_sum(creator);// for debugging

  typedef SemiLinearizationAtInner<SemiImplicitEulerStep<Functional> >  SemiLinearization;
  typedef VariationalFunctionalAssembler<SemiLinearization> Assembler;
  Assembler assembler(spaces);

  // ------------------------------------------------------------------------------------
  // construct mass and stiffness matrix from semi-implicit structure
  // ------------------------------------------------------------------------------------  
  auto du(u);
  for (int subIdx=0; subIdx<sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx]; 
    du *= 0;
    std::map<int,Matrix> massSubmatrices_nbr;
    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());

    F.Mass_stiff(1);
    F.set_mass_submatrix(true);
    SemiImplicitEulerStep<Functional>  eqM(&F,dt);
    eqM.setTau(0);
    for ( const auto & gamma_nbr: nbrs ){
      int row = gamma_nbr.first;
      int col = tag;
      std::string path = std::to_string(col) + "_" + std::to_string(row);

      F.set_row_col_subdomain(row,col);
      assembler.assemble(SemiLinearization(eqM,u,u,du), Assembler::MATRIX,assemblyThreads); 
      Matrix sub_M_ = assembler.template get<Matrix>(false);
      // if(write_to_file) writeToMatlab(assembler,matlab_dir+"/subMatrixM_"+path, "M"); 
      massSubmatrices_nbr[row] = sub_M_;
    }

    Matrix subMatrix(creator);
    std::string path = std::to_string(subIdx+1);

    Vector Fs_petcs_sub =  Fs_petcs[subIdx];

    // construct mass blocks
    Matrix M_sub(creator);
    exctract_petsc_mass_blocks(sequenceOfTags, startingIndexOfTag, sequenceOfsubdomains,
                               map_II,map_GammaGamma,map_GammaNbr,A_,M_,massSubmatrices_nbr,subIdx, M_sub);
    Ms.push_back(M_sub);
    // writeToMatlab(M_sub,rhs_petsc_test,"M_petsc_"+path);
    M_sum+=M_sub;
    
    // -------------------------------------
    // added mass blocks to subMatrix
    // -------------------------------------
    subMatrix+=M_sub;

    // construct stiffness + mass on the current subdomain blocks
    Matrix K_sub(creator);
    exctract_petsc_stiffness_blocks(sequenceOfTags, startingIndexOfTag, map_II,map_GammaGamma, A_,K_,subIdx,K_sub);
    Ks.push_back(K_sub);
    // writeToMatlab(K_sub,rhs_petsc_test,"K_petsc_"+path);
    K_sum+=K_sub;

    // -------------------------------------
    // added stiffness blocks to subMatrix
    // -------------------------------------
    subMatrix+=K_sub;

    subMatrices.push_back(subMatrix);
    subMatrices_M.push_back(M_sub);
    subMatrices_K.push_back(K_sub);
    if(write_to_file) writeToMatlabPath(subMatrix,Fs_petcs_sub,"resultBDDC"+path,matlab_dir, false);
  }

  // writeToMatlab(M_sum,rhs_petsc_test,"M_sum");
  // writeToMatlab(K_sum,rhs_petsc_test,"K_sum");

  M_sum+=K_sum; 
  // writeToMatlabPath(M_sum,rhs_petsc_test,"sum_submatrices",matlab_dir);

  Matrix A_petsc(creator);
  petsc_structure_Matrix(arr_extra_set,sequenceOfTags, startingIndexOfTag, map_II, map_GammaGamma, A_, A_petsc);  
  // if(write_to_file) writeToMatlabPath(A_petsc,rhs_petsc_test,"A_petsc",matlab_dir);
  if(write_to_file) writeToMatlabPath(A_petsc,rhs_petsc_test,"resultBDDC",matlab_dir, true);

  // -------------------------------------
  // compute the weight for each subdomain to update the rhs
  // -------------------------------------
  // - we need to go through the interfaces  
  //  - of the current and their neighbors   
  // - update the rhs for each subdomain
  // - generate the submatrices with the writeToMatlabPath again
  // -------------------------------------
  // then 
  // -------------------------------------

  for (int subIdx=0; subIdx<sequenceOfTags.size(); ++subIdx)
  {
    // std::cout << "============================" <<std::endl;
    // std::cout << "subIdx, weights " << subIdx <<std::endl;
    // std::cout << "============================" <<std::endl;
    std::string path = std::to_string(subIdx+1);
    Matrix subMatrix(creator);
    subMatrix = subMatrices[subIdx];
    Vector Fs_petcs_sub =  rhs_petsc_test;
    Vector weights_sub(subMatrix.N()); 
    // optimize it by iterating only on the interfaces
    for (int k = 0; k < subMatrix.N(); ++k)
    {
      double subvalue = subMatrix[k][k];
      double originalvalue = A_petsc[k][k];
      double weight = (subvalue/originalvalue);
      Fs_petcs_sub[k] = Fs_petcs_sub[k]*weight;
      weights_sub[k] = weight;
      // std::cout << k << " : "<< weight <<std::endl;
    }
    // std::cout << "============================" <<std::endl;
    Fs_petcs[subIdx] = Fs_petcs_sub;
    weights.push_back(weights_sub);
    // if(write_to_file) writeToMatlabPath(subMatrix,Fs_petcs_sub,"resultBDDCNew"+path,matlab_dir);
  }

  return u;
}

void generate_Interror_and_Interfaces_indices(std::vector<int> sequenceOfTags,
                                              std::map<int,std::set<int>> map_II,
                                              std::map<int,std::set<int>> map_GammaGamma,
                                              std::map<int,std::set<int>> map_GammaGamma_W_Nbr,
                                              std::map<std::pair<int,int>, int> map_indices, 
                                              std::string matlab_dir)
{
  double precision = 16;
  for (int i = 0; i < sequenceOfTags.size(); ++i)
  {
    int tag = sequenceOfTags[i]; 
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());

    std::string Ii = matlab_dir+"/BDDC_I"+ std::to_string(i+1) + ".txt";
    std::ofstream Iii(Ii.c_str());
    Iii.precision(precision);

    for (int j = 0; j < Interior.size(); ++j)
    {
      std::pair<int,int> pairs;
      pairs.first = Interior[j];
      pairs.second = tag;

      Iii << map_indices[pairs]+1 << " ";   
    }
  }


  std::string gammaAll = matlab_dir+"/BDDC_Gamma.txt";
  std::ofstream GammaAll(gammaAll.c_str());
  GammaAll.precision(precision);

  for (int i = 0; i < sequenceOfTags.size(); ++i)
  {
    int tag = sequenceOfTags[i];
    std::vector<int> interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
    std::vector<int> gamma_W_Nbr(map_GammaGamma_W_Nbr[tag].begin(), map_GammaGamma_W_Nbr[tag].end());

    std::string gammai = matlab_dir+"/BDDC_Gamma"+ std::to_string(i+1) + ".txt";
    std::ofstream Gammaii(gammai.c_str());
    Gammaii.precision(precision);

    std::string gammai_w_nbr = matlab_dir+"/BDDC_Gamma_w_nbr"+ std::to_string(i+1) + ".txt";
    std::ofstream Gammaii_w_nbr(gammai_w_nbr.c_str());
    Gammaii_w_nbr.precision(precision);

    for (int j = 0; j < gamma_W_Nbr.size(); ++j)
    {
      std::pair<int,int> pairs;
      pairs.first = gamma_W_Nbr[j];
      pairs.second = tag;

      Gammaii_w_nbr << map_indices[pairs]+1 << " "; 
    }

    for (int j = 0; j < interface.size(); ++j)
    {
      std::pair<int,int> pairs;
      pairs.first = interface[j];
      pairs.second = tag;

      Gammaii << map_indices[pairs]+1 << " "; 
      GammaAll << map_indices[pairs]+1 << " ";   
    }
  }
}

template<class Matrix, class Vector>
void construct_As(std::vector<int> sequenceOfTags, 
                  std::map<int,int> startingIndexOfTag,
                  std::map<int,std::set<int>> map_II,
                  std::map<int,std::set<int>> map_GammaGamma,
                  std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                  std::map<int,std::vector<int>> sequanceOfsubdomains,
                  std::map<std::pair<int,int>, int> map_indices_petsc,
                  std::vector<std::set<int>> i2t, 
                  Vector rhs_petsc_test,
                  std::vector<Matrix> subMatrices,
                  std::vector<Matrix> subMatrices_M,
                  std::vector<Matrix> subMatrices_K,
                  std::vector<Matrix> &subMatrices_kaskade,
                  std::vector<Matrix> &subMatrices_kaskade_M,
                  std::vector<Matrix> &subMatrices_kaskade_K)
{

  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx]; 
    std::map<int,std::set<int>> gamma_nbrs_subIdx(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
    std::vector<int> sequanceOfsubdomains_subIdx =  sequanceOfsubdomains[tag];

    Matrix subMatrix  = subMatrices[subIdx];
    Matrix subMatrix_M  = subMatrices_M[subIdx];
    Matrix subMatrix_K  = subMatrices_K[subIdx];

    std::map<std::pair<int,int>, int> map_indices_kaskade;
    {
      int counter_kasakde = 0;
      int start = startingIndexOfTag[tag];

      for (int k = 0; k < map_II[tag].size(); ++k)
      {
        std::pair<int,int> pairs;
        pairs.first = start+k;
        pairs.second = tag;

        map_indices_kaskade[pairs] = counter_kasakde;
        // std::cout << start+k+1 << ":" << counter_kasakde+1 << "\n";
        counter_kasakde++;
      }   
      int shift = map_II[tag].size();
      for (int k = 0; k < map_GammaGamma[tag].size(); ++k)
      {
        std::pair<int,int> pairs;
        pairs.first = start+shift+k;
        pairs.second = tag;

        map_indices_kaskade[pairs] = counter_kasakde;
        // std::cout << start+shift+k+1 << ":" << counter_kasakde+1 << "\n";
        counter_kasakde++;
      }    

      std::vector<int> sequanceOfsubdomains_subIdx =  sequanceOfsubdomains[tag];
      for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
      {
        int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
        std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
        std::vector<int> Interface_nbr_only(gamma_nbrs_subIdx[nbr_Indx].begin(), gamma_nbrs_subIdx[nbr_Indx].end());
    
        int x = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();

        for (int k = 0; k < Interface_nbr_only.size(); ++k)
        {
          std::pair<int,int> pairs;
          pairs.first = Interface_nbr_only[k];
          pairs.second = nbr_Indx;
          
          std::pair<int,int> pairs1;
          pairs1.first = map_indices_petsc[pairs];
          pairs1.second = nbr_Indx;

          map_indices_kaskade[pairs1] = counter_kasakde;    
          //std::cout << map_indices_petsc[Interface_nbr_only[k]] +1 << ":" << counter_kasakde+1<< " ??? " << Interface_nbr_only[k]+1 <<":"<<map_indices_petsc[Interface_nbr_only[k]] +1<<"\n";
          counter_kasakde++;      
        }
      }
    }
    // std::cout <<"\n";
    // // block for subIdx located at 0, 0 of submatrices
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // modified sub matrices which includes only nbr indices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    int Dofs_kaskade = Interior.size() + Interface.size();   
    for (int nbr = 0; nbr < gamma_nbrs_subIdx.size(); ++nbr)
    {
      Dofs_kaskade+=gamma_nbrs_subIdx[nbr].size();
    }
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // creator for kaskade structure, the shrinked version
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    NumaCRSPatternCreator<> creator_kaskade(Dofs_kaskade,Dofs_kaskade,false);
    for (int r=0; r<subMatrix.N(); ++r)
    {
      auto row  = subMatrix[r];

      for (auto ca=row.begin(); ca!=row.end(); ++ca)
      {
        int const c = ca.index();
        // std::cout << "tags_r: " << subMatrix.N() << ": " << row.size() << ", -> " << r << ", -> " << i2t[r].size()<< ", -> " << i2t[c].size()<<std::endl;
        // std::vector<int> tags_r(i2t[r].begin(), i2t[r].end());
        // std::cout << "tags_r: " << tags_r.size() <<std::endl;
        // std::vector<int> tags_c(i2t[c].begin(), i2t[c].end());

        // int row_indx = 0.0;//;map_indices_kaskade[k];
        // int col_indx = 0.0;//map_indices_kaskade[l];
        // std::cout << "tags_r: " << row_indx << ", "<< col_indx <<std::endl;
    //     // creator_kaskade.addElement(row_indx,col_indx);  
      }
    }

    // // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // // fill the sub matrices for kaskade format
    // // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // Matrix subMatrix_kaskade_shrinked(creator_kaskade);  
    // Matrix subMatrix_kaskade_shrinked_M(creator_kaskade);  
    // Matrix subMatrix_kaskade_shrinked_K(creator_kaskade);    
    // {
    //   // block for subIdx located at 0, 0 of submatrices
    //   std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    //   std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
      
    //   // IG for subIdx
    //   int start = startingIndexOfTag[tag];
    //   int size = Interior.size()+Interface.size();
    //   std::vector<int> IG(size); // vector with size ints.
    //   std::iota (std::begin(IG), std::end(IG), start); // Fill with start, 1, ..., size.
    //   // std::cout <<"subIdx:" << subIdx << ", start: " << start << " IG[0]" << IG[0] << " IG[size-1]" << IG[size-1] << std::endl;
    //   auto IGAMMA_block = subMatrix(IG,IG);
    //   insertMatrixBlock(IGAMMA_block, 0, 0, subMatrix_kaskade_shrinked);

    //   auto IGAMMA_block_M = subMatrix_M(IG,IG);
    //   insertMatrixBlock(IGAMMA_block_M, 0, 0, subMatrix_kaskade_shrinked_M);

    //   auto IGAMMA_block_K = subMatrix_K(IG,IG);
    //   insertMatrixBlock(IGAMMA_block_K, 0, 0, subMatrix_kaskade_shrinked_K);
    // }
    // // std::cout << "  rest of diagonal  " <<"\n";
    // // rest of diagonal
    // {     
    //   std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    //   std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

    //   int size = Interior.size() + Interface.size();
    //   int x = size;
    //   int y = size;

    //   for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
    //   {

    //     std::vector<int> Interface_nbr_only(gamma_nbrs_subIdx[nbr].begin(), gamma_nbrs_subIdx[nbr].end()); 
    //     std::vector<int> mapped_Interface_nbr_only(Interface_nbr_only.size()); // vector with size ints.
    //     for (int i = 0; i < Interface_nbr_only.size(); ++i)
    //     {
    //       mapped_Interface_nbr_only[i] = map_indices_petsc[Interface_nbr_only[i]];      
    //     }

    //     int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
    //     std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
    //     std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());
    //     // std::cout << x << " " << y <<"\n";
    //     auto GAMMAGAMMA_diagonal_block = subMatrix(mapped_Interface_nbr_only,mapped_Interface_nbr_only);
    //     insertMatrixBlock(GAMMAGAMMA_diagonal_block, x, y, subMatrix_kaskade_shrinked);
      
    //     auto GAMMAGAMMA_diagonal_block_M = subMatrix_M(mapped_Interface_nbr_only,mapped_Interface_nbr_only);
    //     insertMatrixBlock(GAMMAGAMMA_diagonal_block_M, x, y, subMatrix_kaskade_shrinked_M);

    //     x+= Interface_nbr_only.size();
    //     y+= Interface_nbr_only.size();
    //   }
    // }
    
    // // std::cout << "  cross blocks on column for subIdx  " <<"\n";
    // int x = Interior.size();
    // int y = Interior.size() + Interface.size();

    // int start_x = x;
    // int start_y = y;
    // // cross blocks on column for subIdx 
    // for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
    // {
    //   std::vector<int> Interface_nbr_only(gamma_nbrs_subIdx[nbr].begin(), gamma_nbrs_subIdx[nbr].end()); 
    //   std::vector<int> mapped_Interface_nbr_only(Interface_nbr_only.size()); // vector with size ints.
    //   for (int i = 0; i < Interface_nbr_only.size(); ++i)
    //   {
    //     mapped_Interface_nbr_only[i] = map_indices_petsc[Interface_nbr_only[i]];      
    //   }

    //   int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
    //   std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
    //   std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());

    //   int x_nbr = Interior.size();
    //   // std::cout << start_x << " " << start_y <<"\n";
    //   std::vector<int> G_subIdx(Interface.size()); // vector with size ints.
    //   std::iota (std::begin(G_subIdx), std::end(G_subIdx), startingIndexOfTag[tag] + x_nbr); // Fill with start, 1, ..., size.

    //   auto GAMMAGAMMA_cross_block = subMatrix(G_subIdx, mapped_Interface_nbr_only);
    //   insertMatrixBlockNbr(GAMMAGAMMA_cross_block, start_x, start_y, subMatrix_kaskade_shrinked);

    //   auto GAMMAGAMMA_cross_block_M = subMatrix_M(G_subIdx, mapped_Interface_nbr_only);
    //   insertMatrixBlockNbr(GAMMAGAMMA_cross_block_M, start_x, start_y, subMatrix_kaskade_shrinked_M);

    //   start_y+=mapped_Interface_nbr_only.size();
    // }
    // // std::cout<< "\n\n\n";
    // // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // // save sub_matrices for kaskade format
    // // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // subMatrices_kaskade.push_back(subMatrix_kaskade_shrinked);
    // subMatrices_kaskade_M.push_back(subMatrix_kaskade_shrinked_M);
    // subMatrices_kaskade_K.push_back(subMatrix_kaskade_shrinked_K);
    // std::string path = std::to_string(subIdx);   
  }
}

template<class Vector>
void construct_Fs(std::vector<int> sequenceOfTags, 
                  std::map<int,int> startingIndexOfTag,
                  std::map<int,std::set<int>> map_II,
                  std::map<int,std::set<int>> map_GammaGamma,
                  std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                  std::map<int,std::vector<int>> sequanceOfsubdomains,
                  Vector rhs_petsc_test,
                  std::vector<std::vector<LocalDof>> sharedDofs, 
                  std::vector<Vector> weights,
                  std::map<int, int> map_indices,
                  std::vector<Vector> &Fs)
{

  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {  
    int tag = sequenceOfTags[subIdx];
    std::map<int,std::set<int>> gamma_nbrs_subIdx(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
    std::vector<int> sequanceOfsubdomains_subIdx =  sequanceOfsubdomains[subIdx];
    double coef = 0.5;

    // // block for subIdx located at 0, 0 of submatrices
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // modified sub matrices which inlcludes only nbr indices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    int Dofs_kaskade = Interior.size() + Interface.size();   
    for (int nbr = 0; nbr < gamma_nbrs_subIdx.size(); ++nbr)
    {
      Dofs_kaskade+=gamma_nbrs_subIdx[nbr].size();
    }

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // fill the sub matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    Vector Fs_subIdx(Dofs_kaskade);  
    {
      // block for subIdx located at 0, 0 of submatrices
      std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
      std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
      
      // IG for subIdx
      int start = startingIndexOfTag[subIdx];
      std::vector<int> I(Interior.size()); // vector with size ints.
      std::iota (std::begin(I), std::end(I), start); // Fill with start, 1, ..., size.

      for (int i = 0; i < Interior.size(); ++i)
      {
        Fs_subIdx[i] = rhs_petsc_test[Interior[i]]; // start from 0 index for kaskade structure
      }

      std::vector<int> G(Interface.size()); // vector with size ints.
      std::iota (std::begin(G), std::end(G), Interior.size()); // Fill with start, 1, ..., size.

      for (int i = 0; i < Interface.size(); ++i)
      {
        // if(sharedDofs[Interface[i]].size()>1)
        //     coef = 1.0/sharedDofs[Interface[i]].size();
        coef = weights[subIdx][map_indices[Interface[i]]];
        Fs_subIdx[i+Interior.size()] = coef*rhs_petsc_test[Interface[i]]; // start from 0 index for kaskade structure
      }
    }

    // rest of diagonal
    {     
      std::vector<int> Interior(map_II[subIdx].begin(), map_II[subIdx].end());
      std::vector<int> Interface(map_GammaGamma[subIdx].begin(), map_GammaGamma[subIdx].end());

      int size = Interior.size() + Interface.size();
      int x = size;
      int y = size;

      for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
      {

        std::vector<int> Interface_nbr_only(gamma_nbrs_subIdx[nbr].begin(), gamma_nbrs_subIdx[nbr].end()); 
      
        for (int i = 0; i < Interface_nbr_only.size(); ++i)
        {
          // if(sharedDofs[Interface_nbr_only[i]].size()>1)
          //   coef = 1.0/sharedDofs[Interface_nbr_only[i]].size();
          coef = weights[subIdx][map_indices[Interface_nbr_only[i]]];
          Fs_subIdx[y + i] = coef*rhs_petsc_test[Interface_nbr_only[i]]; // start from y index for kaskade structure
        }
      
        x+= Interface_nbr_only.size();
        y+= Interface_nbr_only.size();
      }
    }
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // save sub_matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    Fs.push_back(Fs_subIdx);
  }
}

#endif