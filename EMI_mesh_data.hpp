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

void subdomain_indices( std::vector<int> sequenceOfTags, 
                        std::map<int,std::set<int>> & map_II,
                        std::map<int,std::set<int>> & map_GammaGamma,
                        std::map<int,std::map<int,std::set<int>>> & map_GammaNbr,
                        std::map<int,std::unordered_map<int, int>> & local2Global,
                        std::map<int,std::unordered_map<int, int>> & global2Local,
                        std::map<int,std::vector<int>> & globalIndices)
{
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
		int tag =  sequenceOfTags[index];
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

void map_kaskade2petcs(std::vector<int> sequenceOfTags, 
                       std::map<int,std::set<int>> & map_II,
                       std::map<int,std::set<int>> & map_GammaGamma,
                       std::map<int, int> & map_indices,
                       std::map<int, int> & map_i2sub)
{
  int counter = 0;
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
   	
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

void compute_sharedDofsKaskade( std::vector<int> sequenceOfTags, 
                                std::map<int, int> map_indices, 
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
      vec_kaskadeIndex[count] = map_indices[I[i]]; 
      count++;     
    }

    for (int i = 0; i < gamma.size(); ++i){
      vec_kaskadeIndex[count] = map_indices[gamma[i]];    
      count++;  
    }
    
    for ( const auto &gamma_nbr : gamma_nbr_subdomain ) 
    {
      int tag_nbr = gamma_nbr.first;
      std::set<int> nbr = gamma_nbr.second;
      std::set<int>::iterator it;

      for (it = nbr.begin(); it != nbr.end(); it++)
      {
        vec_kaskadeIndex[count] = map_indices[*it]; 
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
        for (const auto& value : values) {
          f <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
          // std::cout <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
        }
        f << "\n";
        // std::cout <<"\n";
      }
    }
  } 
}
template<class FSElement>
void write_Dirichlet_and_coordinates( FSElement& fse, 
                                      std::vector<std::vector<int>> e2i,
                                      std::map<int, int> map_indices,
                                      std::vector<std::vector<double>> icoord,
                                      int dof_size,
                                      int dim,
                                      bool write_to_file,
                                      std::string matlab_dir)
{


  std::set<int> dofsDirichlet;
  markedIndicesForDirichlet(fse,
                            GetGlobalCoordinate(), 
                            e2i,
                            dofsDirichlet);



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
    for (int i = 0; i < icoord.size(); ++i)
      indexCoordinates_petsc[map_indices[i]] = icoord[i];

    for (int j = 0; j < indexCoordinates_petsc.size(); ++j){
      if(dim==2) 
        f <<indexCoordinates_petsc[j][0] << " "<< indexCoordinates_petsc[j][1] << " "<< 0.0 << " ";
      if(dim==3) 
        f <<indexCoordinates_petsc[j][0] << " "<< indexCoordinates_petsc[j][1] << " "<< indexCoordinates_petsc[j][2] << " ";
    }
  }
}

void computed_sequenceOfTags(std::map<int, int> map_t2l, std::map<int,std::map<int,std::set<int>>> map_GammaNbr, std::vector<int> & sequenceOfTags, std::vector<int> & startingIndexOfTag, std::map<int,int> & map_nT2oT, std::map<int,std::vector<int>> & sequanceOfsubdomains)
{
  int start = 0;
  int index = 0;
  int nTag = 0;
  for ( const auto &t2l : map_t2l ) 
  {
    int tag =  t2l.first;
    sequenceOfTags[index] = tag;
    startingIndexOfTag[index] = start;
    map_nT2oT[nTag] = tag;

    start+=map_t2l[tag];
    index++;
    nTag++;
  }

  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
    std::vector<int> sequanceOfsubdomains_subIndex;
    for ( const auto &gamma_nbr : nbrs ) {
      sequanceOfsubdomains_subIndex.push_back(gamma_nbr.first);
    }

    sequanceOfsubdomains[tag] = sequanceOfsubdomains_subIndex;
  } 

}

template<class Vector>
void petsc_structure_rhs( std::vector<int> sequenceOfTags, 
                          std::vector<int> startingIndexOfTag,
                          std::map<int,std::set<int>> map_II,
                          std::map<int,std::set<int>> map_GammaGamma, 
                          Vector b_,
                          Vector &bs_)
{
  bs_ *= 0;
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
    std::cout << "tag " << tag <<std::endl;
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());


    int x = startingIndexOfTag[index];
    for (int i = 0; i < Interior.size(); ++i)
    {
      bs_[i+x] = b_[Interior[i]];
    }

    for (int i = 0; i < interface.size(); ++i)
    {
      bs_[i+x+Interior.size()] = b_[interface[i]];
    }
  }
}


template<class Vector>
void petsc_structure_rhs_petsc( std::vector<int> sequenceOfTags, 
                                std::vector<int> startingIndexOfTag,
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma,
                                std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                                Vector b_,
                                std::map<int, int> map_indices,
                                std::vector<std::vector<LocalDof>> sharedDofsAll,
                                std::vector<Vector> &Fs)
{

  std::cout << " construct RHS for each subdomain based on petsc Sructure" <<std::endl;
  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());

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

    for ( const auto &gamma_nbr : nbrs ) 
    {
      int tag = gamma_nbr.first;
      std::set<int> nbr = gamma_nbr.second;
      std::set<int>::iterator it;
      std::vector<int> nbr_vec(nbr.begin(), nbr.end());

      for (int j = 0; j < nbr_vec.size(); ++j)
      {
        int pos = nbr_vec[j];
        int pos_map = map_indices[pos];
        // Fs_subIdx[pos_map] = coef*b_[pos];
        int alpha  = sharedDofsAll[pos].size();
        Fs_subIdx[pos_map] = (1.0/alpha)*b_[pos];
      }
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
      As_[k+x][l+y] = coef*(*ca); 
      As_[l+y][k+x] = coef*(*ca); 
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
                                std::vector<int> startingIndexOfTag,
                                std::map<int,std::vector<int>> sequanceOfsubdomains,
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma,
                                std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                                Matrix A_,
                                Matrix M_,
                                std::vector<Matrix> massSubmatrices_subIndx,
                                int subIdx,
                                Matrix &Ms)
{
  int tag = sequenceOfTags[subIdx];
  std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
  std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
  std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
  
  
  int start = startingIndexOfTag[subIdx];
  double coef = 0.5;
  auto K_GAMMAGAMMA_block = M_(Interface,Interface);
  insertMatrixBlock(K_GAMMAGAMMA_block, coef, start+Interior.size(), start+Interior.size(), Ms);


  std::vector<int> sequanceOfsubdomains_subIdx(sequanceOfsubdomains[subIdx].begin(), sequanceOfsubdomains[subIdx].end());
  for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
  {
    int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
    std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
    std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());
 
    int x = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();
    int y = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();

    Matrix massSubmatrices_nbr = massSubmatrices_subIndx[nbr];
    auto K_GAMMAGAMMA_nbr_nbr_block = massSubmatrices_nbr(Interface_nbr,Interface_nbr);
    insertMatrixBlock(K_GAMMAGAMMA_nbr_nbr_block, coef, x, y, Ms);
  }

  // cross blocks on column for subIdx 
  for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
  {
    int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
    std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
    std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());

    int x = startingIndexOfTag[subIdx] + Interior.size();
    int y = startingIndexOfTag[nbr_Indx] + Interior_nbr.size();

    auto K_GAMMAGAMMA_nbr_block = A_(Interface,Interface_nbr);
    insertMatrixBlockNbr(K_GAMMAGAMMA_nbr_block, coef, x, y, Ms);
  }
}

template<class Matrix, class Matrix_>
void exctract_petsc_stiffness_blocks( std::vector<int> sequenceOfTags, 
                                      std::vector<int> startingIndexOfTag,
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
void petsc_structure_Matrix(std::vector<int> sequenceOfTags, 
                            std::vector<int> startingIndexOfTag,
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

    for (int nbr = 0; nbr < map_II.size(); ++nbr)
    {
      int tag_nbr = sequenceOfTags[nbr];
      if(tag!=tag_nbr){
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
typename VariableSet::VariableSet  construct_submatrices_petsc( std::map<int,int> map_nT2oT,
                                                          GridManager<Grid>& gridManager,
                                                          Functional& F,
                                                          VariableSet const& variableSet, 
                                                          Spaces const& spaces,
                                                          Grid const& grid, 
                                                          typename VariableSet::VariableSet u,
                                                          double  dt,
                                                          std::vector<int> sequenceOfTags, 
                                                          std::vector<int> startingIndexOfTag,
                                                          std::map<int,std::set<int>> map_II,
                                                          std::map<int,std::set<int>> map_GammaGamma,
                                                          std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                                                          std::map<int,std::vector<int>> sequanceOfsubdomains,
                                                          std::map<int, int> map_indices, 
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
  for (int k=0; k<A_.N(); ++k)
  {
    auto row  = A_[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      int row_indx = map_indices[k];
      int col_indx = map_indices[l];
      counter_elements++;
      creator.addElement(row_indx,col_indx);  
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
    std::vector<Matrix> massSubmatrices_nbr;
    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());

    F.Mass_stiff(1);
    F.set_mass_submatrix(true);
    SemiImplicitEulerStep<Functional>  eqM(&F,dt);
    eqM.setTau(0);

    for ( const auto & gamma_nbr: nbrs ){
      int row = gamma_nbr.first;
      int col = subIdx;

      std::string path = std::to_string(col) + "_" + std::to_string(row);

      F.set_row_col_subdomain(map_nT2oT[row],map_nT2oT[col]);
      assembler.assemble(SemiLinearization(eqM,u,u,du), Assembler::MATRIX,assemblyThreads); 
      Matrix sub_M_ = assembler.template get<Matrix>(false);
      // if(write_to_file) writeToMatlab(assembler,matlab_dir+"/subMatrixM_"+path, "M"); 
      massSubmatrices_nbr.push_back(sub_M_);
    }


    Matrix subMatrix(creator);
    std::string path = std::to_string(subIdx+1);

    Vector Fs_petcs_sub =  Fs_petcs[subIdx];

    // construct mass blocks
    Matrix M_sub(creator);
    exctract_petsc_mass_blocks(sequenceOfTags, startingIndexOfTag, sequanceOfsubdomains,
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
  petsc_structure_Matrix(sequenceOfTags, startingIndexOfTag, map_II, map_GammaGamma, A_, A_petsc);  
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
                                              std::map<int, int> map_indices, 
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
      Iii << map_indices[Interior[j]]+1 << " ";   
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
      Gammaii_w_nbr << map_indices[gamma_W_Nbr[j]]+1 << " "; 
    }

    for (int j = 0; j < interface.size(); ++j)
    {
      Gammaii << map_indices[interface[j]]+1 << " "; 
      GammaAll << map_indices[interface[j]]+1 << " ";   
    }
  }
}

template<class Matrix, class Vector>
void construct_As(std::vector<int> sequenceOfTags, 
                  std::vector<int> startingIndexOfTag,
                  std::map<int,std::set<int>> map_II,
                  std::map<int,std::set<int>> map_GammaGamma,
                  std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                  std::map<int,std::vector<int>> sequanceOfsubdomains,
                  std::map<int, int> map_indices_petsc,
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

    std::map<int, int> map_indices_kaskade;
    {
      int counter_kasakde = 0;
      int start = startingIndexOfTag[subIdx];

      for (int k = 0; k < map_II[tag].size(); ++k)
      {
        map_indices_kaskade[start+k] = counter_kasakde;
        // std::cout << start+k+1 << ":" << counter_kasakde+1 << "\n";
        counter_kasakde++;
      }   
      int shift = map_II[tag].size();
      for (int k = 0; k < map_GammaGamma[tag].size(); ++k)
      {
        map_indices_kaskade[start+shift+k] = counter_kasakde;
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
          map_indices_kaskade[map_indices_petsc[Interface_nbr_only[k]]] = counter_kasakde;    
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
    for (int k=0; k<subMatrix.N(); ++k)
    {
      auto row  = subMatrix[k];
      for (auto ca=row.begin(); ca!=row.end(); ++ca)
      {
        int const l = ca.index();
        int row_indx = map_indices_kaskade[k];
        int col_indx = map_indices_kaskade[l];
        creator_kaskade.addElement(row_indx,col_indx);  
      }
    }

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // fill the sub matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    Matrix subMatrix_kaskade_shrinked(creator_kaskade);  
    Matrix subMatrix_kaskade_shrinked_M(creator_kaskade);  
    Matrix subMatrix_kaskade_shrinked_K(creator_kaskade);    
    {
      // block for subIdx located at 0, 0 of submatrices
      std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
      std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());
      
      // IG for subIdx
      int start = startingIndexOfTag[tag];
      int size = Interior.size()+Interface.size();
      std::vector<int> IG(size); // vector with size ints.
      std::iota (std::begin(IG), std::end(IG), start); // Fill with start, 1, ..., size.
      // std::cout <<"subIdx:" << subIdx << ", start: " << start << " IG[0]" << IG[0] << " IG[size-1]" << IG[size-1] << std::endl;
      auto IGAMMA_block = subMatrix(IG,IG);
      insertMatrixBlock(IGAMMA_block, 0, 0, subMatrix_kaskade_shrinked);

      auto IGAMMA_block_M = subMatrix_M(IG,IG);
      insertMatrixBlock(IGAMMA_block_M, 0, 0, subMatrix_kaskade_shrinked_M);

      auto IGAMMA_block_K = subMatrix_K(IG,IG);
      insertMatrixBlock(IGAMMA_block_K, 0, 0, subMatrix_kaskade_shrinked_K);
    }
    // std::cout << "  rest of diagonal  " <<"\n";
    // rest of diagonal
    {     
      std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
      std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

      int size = Interior.size() + Interface.size();
      int x = size;
      int y = size;

      for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
      {

        std::vector<int> Interface_nbr_only(gamma_nbrs_subIdx[nbr].begin(), gamma_nbrs_subIdx[nbr].end()); 
        std::vector<int> mapped_Interface_nbr_only(Interface_nbr_only.size()); // vector with size ints.
        for (int i = 0; i < Interface_nbr_only.size(); ++i)
        {
          mapped_Interface_nbr_only[i] = map_indices_petsc[Interface_nbr_only[i]];      
        }

        int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
        std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
        std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());
        // std::cout << x << " " << y <<"\n";
        auto GAMMAGAMMA_diagonal_block = subMatrix(mapped_Interface_nbr_only,mapped_Interface_nbr_only);
        insertMatrixBlock(GAMMAGAMMA_diagonal_block, x, y, subMatrix_kaskade_shrinked);
      
        auto GAMMAGAMMA_diagonal_block_M = subMatrix_M(mapped_Interface_nbr_only,mapped_Interface_nbr_only);
        insertMatrixBlock(GAMMAGAMMA_diagonal_block_M, x, y, subMatrix_kaskade_shrinked_M);

        x+= Interface_nbr_only.size();
        y+= Interface_nbr_only.size();
      }
    }
    
    // std::cout << "  cross blocks on column for subIdx  " <<"\n";
    int x = Interior.size();
    int y = Interior.size() + Interface.size();

    int start_x = x;
    int start_y = y;
    // cross blocks on column for subIdx 
    for (int nbr = 0; nbr < sequanceOfsubdomains_subIdx.size(); ++nbr)
    {
      std::vector<int> Interface_nbr_only(gamma_nbrs_subIdx[nbr].begin(), gamma_nbrs_subIdx[nbr].end()); 
      std::vector<int> mapped_Interface_nbr_only(Interface_nbr_only.size()); // vector with size ints.
      for (int i = 0; i < Interface_nbr_only.size(); ++i)
      {
        mapped_Interface_nbr_only[i] = map_indices_petsc[Interface_nbr_only[i]];      
      }

      int nbr_Indx = sequanceOfsubdomains_subIdx[nbr];
      std::vector<int> Interior_nbr(map_II[nbr_Indx].begin(), map_II[nbr_Indx].end());
      std::vector<int> Interface_nbr(map_GammaGamma[nbr_Indx].begin(), map_GammaGamma[nbr_Indx].end());

      int x_nbr = Interior.size();
      // std::cout << start_x << " " << start_y <<"\n";
      std::vector<int> G_subIdx(Interface.size()); // vector with size ints.
      std::iota (std::begin(G_subIdx), std::end(G_subIdx), startingIndexOfTag[tag] + x_nbr); // Fill with start, 1, ..., size.

      auto GAMMAGAMMA_cross_block = subMatrix(G_subIdx, mapped_Interface_nbr_only);
      insertMatrixBlockNbr(GAMMAGAMMA_cross_block, start_x, start_y, subMatrix_kaskade_shrinked);

      auto GAMMAGAMMA_cross_block_M = subMatrix_M(G_subIdx, mapped_Interface_nbr_only);
      insertMatrixBlockNbr(GAMMAGAMMA_cross_block_M, start_x, start_y, subMatrix_kaskade_shrinked_M);

      start_y+=mapped_Interface_nbr_only.size();
    }
    // std::cout<< "\n\n\n";
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // save sub_matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    subMatrices_kaskade.push_back(subMatrix_kaskade_shrinked);
    subMatrices_kaskade_M.push_back(subMatrix_kaskade_shrinked_M);
    subMatrices_kaskade_K.push_back(subMatrix_kaskade_shrinked_K);
    std::string path = std::to_string(subIdx);   
  }
}

template<class Vector>
void construct_Fs(std::vector<int> sequenceOfTags, 
                  std::vector<int> startingIndexOfTag,
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