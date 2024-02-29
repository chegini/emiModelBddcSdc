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

struct vec {
    size_t N;
    std::vector<double> data;

    vec(size_t N) : N(N), data(N) {};

    vec(std::ifstream &in) {
        in.read(reinterpret_cast<char*>(&N), sizeof(size_t));
        data.resize(N);
        in.read(reinterpret_cast<char*>(data.data()), N * sizeof(double));
    }

    void write(std::ofstream &out) {
        out.write(reinterpret_cast<char*>(&N), sizeof(size_t));
        out.write(reinterpret_cast<char*>(data.data()), N * sizeof(double));
    }
};

struct csr_matrix {
    size_t N, nnz;
    std::vector<size_t> row_ptrs;
    std::vector<size_t> col_idxs;
    std::vector<double> values;

    csr_matrix(size_t N, size_t nnz) : N(N), nnz(nnz), row_ptrs(N + 1, 0), col_idxs(nnz), values(nnz)
    {}

    csr_matrix(std::ifstream &in) {
        in.read(reinterpret_cast<char*>(&N), sizeof(size_t));
        in.read(reinterpret_cast<char*>(&nnz), sizeof(size_t));
        row_ptrs.resize(N + 1);
        col_idxs.resize(nnz);
        values.resize(nnz);
        in.read(reinterpret_cast<char*>(row_ptrs.data()), (N + 1) * sizeof(size_t));
        in.read(reinterpret_cast<char*>(col_idxs.data()), nnz * sizeof(size_t));
        in.read(reinterpret_cast<char*>(values.data()), nnz * sizeof(double));
    }

    void write(std::ofstream &out) {
        out.write(reinterpret_cast<char*>(&N), sizeof(size_t));
        out.write(reinterpret_cast<char*>(&nnz), sizeof(size_t));
        out.write(reinterpret_cast<char*>(row_ptrs.data()), (N + 1) * sizeof(size_t));
        out.write(reinterpret_cast<char*>(col_idxs.data()), nnz * sizeof(size_t));
        out.write(reinterpret_cast<char*>(values.data()), nnz * sizeof(double));
    }
};

template <class Entry, class Index, class VEntry>
void writeToMatlabPath(NumaBCRSMatrix<Entry,Index> const& A, Dune::BlockVector<VEntry> const& b, std::string const& basename, std::string const& path, bool gen_rhs, int precision=16)
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
                          std::map<int,std::set<int>> & map_GammaGamma,
                          std::map<int,std::set<int>> & map_GammaGamma_W_Nbr,
                          std::map<int,std::map<int,std::set<int>>> & map_GammaNbr)
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
                                            map_GammaGamma,
                                            map_GammaGamma_W_Nbr,
                                            map_GammaNbr);



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

		std::cout << "t2l"<< ":\n";
		for ( const auto &t2l : map_t2l ) {
			std::cout << t2l.first << ": " << t2l.second << "\n";
		}
  }
}

void local2GlobalMapSubdomain(std::vector<int> I, std::vector<int>  gamma, 
															std::map<int,std::set<int>> gamma_nbr_subdomain,
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

	for ( const auto &gamma_nbr : gamma_nbr_subdomain ) {
		int tag = gamma_nbr.first;
		std::set<int> nbr = gamma_nbr.second;
		std::set<int>::iterator it;

    for (it = nbr.begin(); it != nbr.end(); it++)
    {
    	local2Global_subdomain.insert(pair<int, int>(index,*it)); 
      index++; 
    }
	}
}

void global2LocalMapSubdomain(std::vector<int> I, std::vector<int>  gamma, 
															std::map<int,std::set<int>> gamma_nbr_subdomain,
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

	for ( const auto &gamma_nbr : gamma_nbr_subdomain ) {
		int tag = gamma_nbr.first;
		std::set<int> nbr = gamma_nbr.second;
		std::set<int>::iterator it;

    for (it = nbr.begin(); it != nbr.end(); it++)
    {
    	global2Local_subdomain.insert(pair<int, int>(index,*it)); 
      index++; 
    }
	}
}

void globalIndicesSubdomain(std::vector<int> I, std::vector<int>  gamma, 
														std::map<int,std::set<int>> gamma_nbr_subdomain,
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

	for ( const auto &gamma_nbr : gamma_nbr_subdomain ) {
		int tag = gamma_nbr.first;
		std::set<int> nbr = gamma_nbr.second;
		std::set<int>::iterator it;
		std::vector<int> nbr_vec(nbr.begin(), nbr.end());

		for (int j = 0; j < nbr_vec.size(); ++j)
    {
      globalIndices_subdomain.push_back(nbr_vec[j]);
      index++; 
    }
	}
}

void subdomain_indices( std::map<int,std::set<int>> & map_II,
                        std::map<int,std::set<int>> & map_GammaGamma,
                        std::map<int,std::map<int,std::set<int>>> & map_GammaNbr,
                        std::map<int,std::unordered_map<int, int>> & local2Global,
                        std::map<int,std::unordered_map<int, int>> & global2Local,
                        std::map<int,std::vector<int>> & globalIndices)
{
	for ( const auto &II : map_II ) {
		int tag =  II.first;
		std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
		std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());

    std::unordered_map<int, int> local2Global_subIdx;
    std::unordered_map<int, int> global2Local_subIdx;
    std::vector<int> globalIndices_subIdx;

    local2GlobalMapSubdomain(Interior,Interface,nbrs,local2Global_subIdx);
    global2LocalMapSubdomain(Interior,Interface,nbrs,global2Local_subIdx);
    globalIndicesSubdomain(Interior,Interface,nbrs,globalIndices_subIdx);
    local2Global[tag] = local2Global_subIdx;
    global2Local[tag] = global2Local_subIdx;
    globalIndices[tag] = globalIndices_subIdx;
	}
}

void map_kaskade2petcs(std::map<int,std::set<int>> & map_II,
                               std::map<int,std::set<int>> & map_GammaGamma,
                               std::map<int, int> & map_indices,
                               std::map<int, int> & map_i2sub)
{
  int counter = 0;
	for ( const auto &II : map_II ) {
		int tag =  II.first;
   	
   	std::vector<int> I_vec(map_II[tag].begin(),map_II[tag].end());
   	std::vector<int> gamma_vec(map_GammaGamma[tag].begin(),map_GammaGamma[tag].end());

    for (int i = 0; i < I_vec.size(); ++i)
    {
      map_indices[I_vec[i]] = counter;
      map_i2sub[I_vec[i]] = tag;
      counter++;
    }

    for (int i = 0; i < gamma_vec.size(); ++i)
    {
      map_indices[gamma_vec[i]] = counter;
      map_i2sub[gamma_vec[i]] = tag;
      counter++;
    }
  }

  std::vector<int> sequanceIndices_petsc;
  if(false){
    for ( const auto &II : map_II ) {
      int tag =  II.first;
      std::vector<int> I_vec(map_II[tag].begin(),map_II[tag].end());
      std::vector<int> gamma_vec(map_GammaGamma[tag].begin(),map_GammaGamma[tag].end());

      for (int i = 0; i < I_vec.size(); ++i)
      {
        sequanceIndices_petsc.push_back(map_indices[I_vec[i]]);
        std::cout <<I_vec[i] << " " << map_indices[I_vec[i]] << "\n" ;
      }

      for (int i = 0; i < gamma_vec.size(); ++i)
      {
        sequanceIndices_petsc.push_back(map_indices[gamma_vec[i]]);
        std::cout <<gamma_vec[i] << " " << map_indices[gamma_vec[i]] << "\n" ;
      }
      std::cout << "\n" ;
    }    
  }
}

void removeInnerIndices_i2i(std::vector<std::set<int>> & i2i)
{
  for (int i = 0; i < i2i.size(); ++i)
  {
    std::set<int> s = i2i[i];
    if(s.size()==1) {
      i2i.erase(i2i.begin()+i);
    }
  }
}

#endif