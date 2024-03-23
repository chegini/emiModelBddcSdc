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
                          std::set<int>& tags,                                    
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
                          std::set<int> & interface_extra_dofs,
                          std::map<int,std::vector<int>> & sequenceOfsubdomains)     
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
    tags.insert(tag);
    std::vector<int> igamma(IGamma.second.begin(), IGamma.second.end());
    std::vector<int> gamma(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

    std::set<int> I;
    std::set_difference(IGamma.second.begin(), IGamma.second.end(), 
                        map_GammaGamma[tag].begin(), map_GammaGamma[tag].end(),
                        std::inserter(I, I.end()));

    map_II[IGamma.first] = I;
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

  // REMOVE THE INTERFACES FROM ONE OF THE EXTERNAL, 
  // IT MEANS ONLY ONE EXTRA CELLULAR SUBDOMAIN WILL HAVE THE DOFS ON THE COMMON INTERFACES BETWEEN DIFFERENT EXTRACELLULAR SUBDOMAINS
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
        }
        count++;
      }
    }else{
      std::vector<int> v(s.begin(),s.end());
      i2Tag[i] = v[0];
    }
  }

  // GET only neighbrs without duplication the extra cellular
  // question, I should make them be,ong to one subdomain???
  for ( const auto &Gamma : map_IGamma_noDuplicate ) 
  {
    int tag = Gamma.first;
    std::set<int> difference;
    std::set<int>::iterator itr;
    std::set<int> tags;
    // Use std::set_difference to find the difference between set1 and set2
    std::set_difference(map_GammaGamma_W_Nbr[tag].begin(), map_GammaGamma_W_Nbr[tag].end(),
                        map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end(),
                        std::inserter(difference, difference.begin()));
 
    for (itr = difference.begin(); itr != difference.end(); itr++)
    {
      tags.insert(i2Tag[*itr]);
    }
    std::vector<int> tags_vec(tags.begin(), tags.end());
    sequenceOfsubdomains[tag] = tags_vec;
    map_GammaNbr_Nbr_noDuplicate[tag] = difference;
  }

  int n_count = 0;
  {
    std::map<int, int>::iterator it;

    for (it = map_t2l.begin(); it != map_t2l.end(); it++)
    {
      map_sT2l[n_count] = it->second; 
      if(false){
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

  if(false)
  {
    std::cout << "\ni2Tag "<< std::endl;
    for (int i = 0; i < i2Tag.size(); ++i)
    {
        std::cout << i << " -> " <<i2Tag[i]<< std::endl;
    }

    std::cout << "\ni2t "<< std::endl;
    for (int i = 0; i < i2t.size(); ++i)
    {
        std::cout << i << " -> " <<i2t[i].size()<< std::endl;
    }

    std::cout << "\n";
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

    std::cout << "\ninterface_extra_dofs" << interface_extra_dofs.size()<< ":\n";
    {
    std::set<int>::iterator it;
    for (it = interface_extra_dofs.begin(); it != interface_extra_dofs.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
    }
  
    std::cout << "\nigamma"<< ":\n";
    for ( const auto &IGamma : map_IGamma ) {
    std::cout << IGamma.first << ": ";

    std::set<int>::iterator it;
    for (it = IGamma.second.begin(); it != IGamma.second.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
    }

    std::cout << "\nigamma_noDuplicate"<< ":\n";
    for ( const auto &IGamma : map_IGamma_noDuplicate ) {
    std::cout << IGamma.first << ": ";

    std::set<int>::iterator it;
    for (it = IGamma.second.begin(); it != IGamma.second.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
    }

    std::cout << "\ngamma"<< ":\n";
    for ( const auto &GammaGamma : map_GammaGamma ) {
    std::cout << GammaGamma.first << ": ";

    std::set<int>::iterator it;
    for (it = GammaGamma.second.begin(); it != GammaGamma.second.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
    }

    std::cout << "\nmap_GammaGamma_noDuplicate"<< ":\n";
    for ( const auto &GammaGamma : map_GammaGamma_noDuplicate ) {
    std::cout << GammaGamma.first << ": ";

    std::set<int>::iterator it;
    for (it = GammaGamma.second.begin(); it != GammaGamma.second.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
    }

    std::cout << "\nI"<< ":\n";
    for ( const auto &II : map_II ) {
    std::cout << II.first << ": ";

    std::set<int>::iterator it;
    for (it = II.second.begin(); it != II.second.end(); ++it) {
       std::cout << *it << " ";
    }
    std::cout << "\n";
    }

    std::cout << "\ngamma with nbr"<< ":\n";
    for ( const auto &Gamma_W_Nbr : map_GammaGamma_W_Nbr ){
      std::cout << Gamma_W_Nbr.first << ": ";
      std::set<int>::iterator itr;
      for (itr = Gamma_W_Nbr.second.begin(); itr != Gamma_W_Nbr.second.end(); itr++) 
      {
        std::cout << *itr << " ";
      }
      std::cout << "\n";
    } 

    std::cout << "\nnbr"<< ":\n";
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

    std::cout << "\nonly Gamma nbrs(diff) " <<std::endl;
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

    std::cout << "\nt2l"<< ":\n";
    for ( const auto &t2l : map_t2l ) {
      std::cout << t2l.first << ": " << t2l.second << "\n";
    }
  }
}

void local2GlobalMapSubdomain(std::vector<int> I, std::vector<int>  gamma, std::vector<int> Interface_nbr, 
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

  // for (int i = 0; i < Interface_nbr.size(); ++i)
  // {
  //   local2Global_subdomain.insert(pair<int, int>(index,Interface_nbr[i]));     
  //   index++;  
  // }
}

void global2LocalMapSubdomain(std::vector<int> I, std::vector<int>  gamma, std::vector<int> Interface_nbr, 
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

  // for (int i = 0; i < Interface_nbr.size(); ++i)
  // {
  //   global2Local_subdomain.insert(pair<int, int>(Interface_nbr[i], index));     
  //   index++;  
  // }
}

void globalIndicesSubdomain(std::vector<int> I, std::vector<int>  gamma, std::vector<int> Interface_nbr, 
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

  // for (int i = 0; i < Interface_nbr.size(); ++i)
  // {
  //   globalIndices_subdomain.push_back(Interface_nbr[i]); 
  //   index++;  
  // }
}


void subdomain_indices( std::vector<int> sequenceOfTags, 
                        std::map<int,std::set<int>> & map_II,
                        std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                        std::map<int,std::set<int>> & map_GammaNbr_Nbr_noDuplicate, 
                        std::map<int,std::unordered_map<int, int>> & local2Global,
                        std::map<int,std::unordered_map<int, int>> & global2Local,
                        std::map<int,std::vector<int>> & globalIndices)
{
  for (int index = 0; index < sequenceOfTags.size(); ++index)
  { 
    int tag =  sequenceOfTags[index];
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> Interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
    std::vector<int> Interface_nbr(map_GammaNbr_Nbr_noDuplicate[tag].begin(), map_GammaNbr_Nbr_noDuplicate[tag].end());

    std::unordered_map<int, int> local2Global_subIdx;
    std::unordered_map<int, int> global2Local_subIdx;
    std::vector<int> globalIndices_subIdx;

    local2GlobalMapSubdomain(Interior,Interface,Interface_nbr,local2Global_subIdx);
    global2LocalMapSubdomain(Interior,Interface,Interface_nbr,global2Local_subIdx);
    globalIndicesSubdomain(Interior,Interface,Interface_nbr,globalIndices_subIdx);
    local2Global[tag] = local2Global_subIdx;
    global2Local[tag] = global2Local_subIdx;
    globalIndices[tag] = globalIndices_subIdx;

    if(false){
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
      // std::cout << I_vec[i] << ": "  << counter <<std::endl;
      counter++;
    }

    for (int i = 0; i < gamma_vec.size(); ++i)
    {
      std::pair<int,int> pairs;
      pairs.first = gamma_vec[i];
      pairs.second = tag;

      map_indices[gamma_vec[i]] = counter;
      // std::cout << gamma_vec[i] << ": "  << counter <<std::endl;
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
                                      std::vector<int> arr_extra,
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
                            arr_extra,
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
    for (int i = 0; i < coord_globalIndex.size(); ++i)
      indexCoordinates_petsc[map_indices[i]] = coord_globalIndex[i];

    for (int j = 0; j < indexCoordinates_petsc.size(); ++j){
      if(dim==2) 
        f << indexCoordinates_petsc[j][0] << " "<< indexCoordinates_petsc[j][1] << " "<< 0.0 << "\n";
      if(dim==3) 
        f <<indexCoordinates_petsc[j][0] << " "<< indexCoordinates_petsc[j][1] << " "<< indexCoordinates_petsc[j][2] << "\n";
    }
  }

}

void computed_sequenceOfTags(std::map<int, int> map_t2l, 
                             std::map<int,std::set<int>> map_IGamma_noDuplicate,
                             std::map<int,std::map<int,std::set<int>>> map_GammaNbr, 
                             std::vector<int> & sequenceOfTags, 
                             std::map<int,int> & startingIndexOfTag, 
                             std::map<int,int> & map_nT2oT)
{  
  int start = 0;
  int index = 0;
  int nTag = 0;
  for ( const auto &IGamma : map_IGamma_noDuplicate ) {
    int tag = IGamma.first;
    sequenceOfTags[index] = tag;
    startingIndexOfTag[tag] = start;
    // std::cout << "tag: " << tag << " start: "<< start <<std::endl;
    start+=IGamma.second.size() ;
    map_nT2oT[nTag] = tag;
    nTag++;
    index++;
  }

  // for (int index = 0; index < sequenceOfTags.size(); ++index)
  // { 
  //   int tag =  sequenceOfTags[index];
  //   std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
  //   std::vector<int> sequenceOfsubdomains_subIndex;
  //   // std::cout << "tag: " << tag << ": ";
  //   for ( const auto &gamma_nbr : nbrs ) {
  //     // std::cout << gamma_nbr.first << " ";
  //     sequenceOfsubdomains_subIndex.push_back(gamma_nbr.first);
  //   }
  //   // std::cout << "\n";
  //   sequenceOfsubdomains[tag] = sequenceOfsubdomains_subIndex;
  // } 
}

void marked_corners(std::vector<int> arr_extra, std::vector<int> sequenceOfTags, std::map<int,int> map_indices, std::vector<std::set<int>> i2t, 
                    std::set<std::set<int>> i2iSet, std::map<int,std::set<int>> map_GammaNbr_Nbr,  std::string matlab_dir, std::map<int,bool> & map_markCorners)
{

  std::vector<std::set<int>> i2t_all = i2t;
  bool more_extras = arr_extra.size();
  std::set<int> arr_extra_set(arr_extra.begin(), arr_extra.end());
  std::set<int>::iterator itr;
  std::set<int>::iterator itr_tags;

  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];

    std::vector<int> gammaNbr(map_GammaNbr_Nbr[tag].begin(), map_GammaNbr_Nbr[tag].end());
    for (int i = 0; i < gammaNbr.size(); ++i)
    {
      std::set<int> s = i2t_all[gammaNbr[i]]; 
      itr = s.find(tag);
      if (itr == s.end())
      {
        s.insert(tag);
        i2t_all[gammaNbr[i]] = s;
      }
    }
  }

  for (int i = 0; i < i2t_all.size(); ++i)
  {
    std::set<int> tags = i2t_all[i];
    int count_extra = 0;
    int count_intra = 0;
    // std::cout << i << ": ";
    for (itr_tags = tags.begin(); itr_tags != tags.end(); ++itr_tags) {
      itr = arr_extra_set.find(*itr_tags );
      if (itr != arr_extra_set.end())
      {
        count_extra++;
      }else{
        count_intra++;
      }

      // std::cout << *itr_tags << " ";
    }
    if(count_intra>0 and count_extra>1)
    {
      map_markCorners[i] = true;
    }
  }
 
}

template<class Vector>
void petsc_structure_rhs( std::vector<int> sequenceOfTags,
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
    for (int i = 0; i < Interior.size(); ++i)
    {
      bs_[map_indices[Interior[i]]] = b_[Interior[i]];
    }

    for (int i = 0; i < interface.size(); ++i)
    {
      bs_[map_indices[interface[i]]] = b_[interface[i]];
    }
  }
}


template<class Vector>
void petsc_structure_rhs_subdomain_petsc( std::vector<int> sequenceOfTags, 
                                std::map<int,std::set<int>> map_II,
                                std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                                Vector b_,
                                std::map<int, int> map_indices,
                                std::vector<std::vector<LocalDof>> sharedDofsAll,
                                std::vector<Vector> &Fs)
{
  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());

    Vector Fs_subIdx(b_.size()); 
    double coef = 0.5;

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
      // int alpha  = sharedDofsAll[pos].size();
      Fs_subIdx[pos_map] = coef*b_[pos];
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
void insertMatrixBlock_extra_exclude(Matrix A_block, std::map<int,int> map_indices ,std::vector<int> i2Tag, int tag, std::vector<int> Interior, std::vector<int> Interface, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      if(i2Tag[Interior[k]] != tag and i2Tag[Interface[l]] != tag) {
        continue;
      }else{
        As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = *ca; 
      }
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlock_extra(Matrix A_block, std::map<int,int> map_indices ,std::vector<int> Interior, std::vector<int> Interface, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();

      As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = *ca; 
    }
  }
}


template<class Matrix, class Matrix_>
void insertMatrixBlock_extra(Matrix A_block, std::map<int,int> map_indices , double coef, std::vector<int> Interior, std::vector<int> Interface, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();

      As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = coef*(*ca); 
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlock_mass_extra(Matrix A_block, std::map<int,int> map_indices , double coef, std::vector<int> Interior, std::vector<int> Interface, std::map<int,bool> map_markCorners, std::vector<int> i2Tag, int tag ,Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    for (auto ca=row.begin(); ca!=row.end(); ++ca)
    {
      int const l = ca.index();
      if(Interior[k]==Interface[l] and map_markCorners[Interface[l]] and tag!=i2Tag[Interface[l]]){
        // std::cout << Interior[k]<< " &&&& " << Interface[l] << " tag: " << tag << " i2Tag[Interface[l]] " << i2Tag[Interface[l]] <<" coef*(*ca)" << coef*(*ca) << std::endl;
      }
      else 
        As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = coef*(*ca); 
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlock_extra(Matrix A_block, std::map<int,int> map_indices, std::vector<int> i2Tag, int tag, std::vector<int> Interior, std::vector<int> Interface, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    {
      for (auto ca=row.begin(); ca!=row.end(); ++ca)
      {
        int const l = ca.index();
        if(i2Tag[Interior[k]] == tag and i2Tag[Interface[l]] == tag) {
          As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = *ca;
        }else if((i2Tag[Interior[k]] == tag and i2Tag[Interface[l]] != tag) or
                 (i2Tag[Interior[k]] != tag and i2Tag[Interface[l]] == tag) ){
          As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = *ca; 
        }
      }
    }
  }
}

template<class Matrix, class Matrix_>
void insertMatrixBlock_extra(Matrix A_block, std::map<int,int> map_indices, std::vector<int> i2Tag, int tag, double coef, std::vector<int> Interior, std::vector<int> Interface, Matrix_ &As_)
{
  for (int k = 0; k < A_block.N(); ++k)
  {
    auto row  = A_block[k];
    {
      for (auto ca=row.begin(); ca!=row.end(); ++ca)
      {
        int const l = ca.index();
        if(i2Tag[Interior[k]] == tag and i2Tag[Interface[l]] == tag) {
          As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = *ca;
        }else if((i2Tag[Interior[k]] == tag and i2Tag[Interface[l]] != tag) or
                 (i2Tag[Interior[k]] != tag and i2Tag[Interface[l]] == tag) ){
          As_[map_indices[Interior[k]]][map_indices[Interface[l]]] = coef*(*ca); 
        }
      }
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


template<class Matrix, class Matrix_>
void exctract_petsc_stiffness_blocks_moreExtracellular( std::vector<int> sequenceOfTags, 
                                            std::map<int, int> map_indices,
                                            std::map<int,std::set<int>> map_II,
                                            std::map<int,std::set<int>> map_GammaGamma,
                                            Matrix K_,
                                            std::vector<int> i2Tag,
                                            int subIdx,
                                            Matrix_ &Ks_)
{
  int tag = sequenceOfTags[subIdx];
  std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
  std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

  auto A_II_block = K_(Interior,Interior);
  insertMatrixBlock_extra(A_II_block, map_indices, Interior, Interior, Ks_);

  auto A_IGAMMA_block = K_(Interior,Interface);
  insertMatrixBlock_extra(A_IGAMMA_block, map_indices, Interior, Interface, Ks_);

  auto A_GAMMAI_block = K_(Interface,Interior);
  insertMatrixBlock_extra(A_GAMMAI_block, map_indices, Interface, Interior, Ks_);

  auto K_GAMMAGAMMA_block = K_(Interface,Interface);
  insertMatrixBlock_extra(K_GAMMAGAMMA_block, map_indices, Interface, Interface, Ks_);
  // insertMatrixBlock_extra_exclude(K_GAMMAGAMMA_block, map_indices, i2Tag, tag, Interface, Interface, Ks_);
}

template<class Matrix, class Matrix_>
void exctract_petsc_stiffness_blocks( std::vector<int> sequenceOfTags, 
                                      std::map<int,int> startingIndexOfTag,
                                      std::map<int, int> map_indices,
                                      std::map<int,std::set<int>> map_II,
                                      std::map<int,std::set<int>> map_GammaGamma,
                                      Matrix K_,
                                      int subIdx,
                                      Matrix_ &Ks_)
{
  int tag = sequenceOfTags[subIdx];
  std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
  std::vector<int> Interface(map_GammaGamma[tag].begin(), map_GammaGamma[tag].end());

  int start = startingIndexOfTag[tag];
  auto A_II_block = K_(Interior,Interior);
  insertMatrixBlock(A_II_block, start, start, Ks_);

  auto A_IGAMMA_block = K_(Interior,Interface);
  insertMatrixBlock(A_IGAMMA_block, start, start+Interior.size(), Ks_);

  auto A_GAMMAI_block = K_(Interface,Interior);
  insertMatrixBlock(A_GAMMAI_block, start+Interior.size(), start, Ks_);

  auto K_GAMMAGAMMA_block = K_(Interface,Interface);
  insertMatrixBlock(K_GAMMAGAMMA_block, start+Interior.size(), start+Interior.size(), Ks_);
}

// -------------------------------------------------------------------------------------
// reorder matrix A with respect to petsc structure
// -------------------------------------------------------------------------------------
template<class Matrix>
void petsc_structure_Matrix( std::vector<int> arr_extra, 
                                  std::vector<int> sequenceOfTags, 
                                  std::map<int,int> startingIndexOfTag,
                                  std::map<int,std::set<int>> map_II,
                                  std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                                  Matrix A_,
                                  Matrix &As_)
{ 
  bool more_extras = arr_extra.size();
  std::set<int> arr_extra_set(arr_extra.begin(), arr_extra.end());
  std::set<int>::iterator itr;
  std::set<int>::iterator itr_nbr;
  for (int subIdx = 0; subIdx < map_II.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx];
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());

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

    for (int nbr = 0; nbr <  map_II.size(); ++nbr)
    {
      int tag_nbr =  sequenceOfTags[nbr]; 
      if(tag!=tag_nbr){
        std::vector<int> Interior_nbr(map_II[tag_nbr].begin(), map_II[tag_nbr].end());
        std::vector<int> interface_nbr(map_GammaGamma_noDuplicate[tag_nbr].begin(), map_GammaGamma_noDuplicate[tag_nbr].end());

        x = startingIndexOfTag[tag] + Interior.size();
        y = startingIndexOfTag[tag_nbr] + Interior_nbr.size();
        auto A_Interface_Interface_nbr_block = A_(interface,interface_nbr);
        insertMatrixBlockNbr(A_Interface_Interface_nbr_block, x, y, As_); 
      }
    }

    if(more_extras){
      for (int nbr = 0; nbr <  map_II.size(); ++nbr)
      {
        int tag_nbr =  sequenceOfTags[nbr]; 
        bool extra_cells = false;
        itr = arr_extra_set.find(tag);
        itr_nbr = arr_extra_set.find(tag_nbr);
        if(itr!= arr_extra_set.end() and itr_nbr!= arr_extra_set.end())
        {
          extra_cells = true;
        }

        if(tag!=tag_nbr and extra_cells){
          std::vector<int> Interior_nbr(map_II[tag_nbr].begin(), map_II[tag_nbr].end());
          std::vector<int> interface_nbr(map_GammaGamma_noDuplicate[tag_nbr].begin(), map_GammaGamma_noDuplicate[tag_nbr].end());

          x = startingIndexOfTag[tag] + Interior.size();
          y = startingIndexOfTag[tag_nbr];
          auto A_Interface_Interior_nbr_block = A_(interface,Interior_nbr);
          insertMatrixBlockNbr(A_Interface_Interior_nbr_block, x, y, As_); 

          x = startingIndexOfTag[tag_nbr];
          y = startingIndexOfTag[tag]+ Interior.size();
          A_Interface_Interior_nbr_block = A_(Interior_nbr,interface);
          insertMatrixBlockNbr(A_Interface_Interior_nbr_block, x, y, As_); 
        }
      }
    }
  }
}

template<class Grid, class Functional, class CellFilter, class VariableSet, class Spaces, class Matrix, class Vector>
typename VariableSet::VariableSet  construct_submatrices_petsc( std::vector<int> arr_extra, 
                                                          std::map<int,int> map_nT2oT,
                                                          GridManager<Grid>& gridManager,
                                                          Functional& F,
                                                          CellFilter & Cellfltr,
                                                          VariableSet const& variableSet, 
                                                          Spaces const& spaces,
                                                          Grid const& grid, 
                                                          typename VariableSet::VariableSet u,
                                                          double  dt,
                                                          std::vector<int> sequenceOfTags, 
                                                          std::map<int,int> startingIndexOfTag,
                                                          std::map<int,std::set<int>> map_II,
                                                          std::map<int,std::set<int>> map_GammaGamma,
                                                          std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                                                          std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                                                          std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate,
                                                          std::map<int,std::vector<int>> sequenceOfsubdomains,
                                                          std::map<int, int> map_indices,
                                                          std::map<int,bool> map_markCorners,
                                                          std::set<int> cells_set,
                                                          std::set<int> tags,
                                                          std::vector<int> i2Tag,
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
                                                          std::map<int, Vector> &weights,
                                                          std::map<int, Matrix> &subMatrices,
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

  // original lhs in petsc format
  Matrix A_petsc(creator);
  {
    for (int k = 0; k < A_.N(); ++k)
    {
      auto row_ = A_[k];
      for (auto ca=row_.begin(); ca!=row_.end(); ++ca)
      {
        int const l = ca.index();
        A_petsc[map_indices[k]][map_indices[l]] = (*ca); 
      }
    }
  }

  writeToMatlabPath(A_petsc,rhs_petsc_test,"resultBDDC",matlab_dir, true);

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

  std::set<int> arr_extra_set(arr_extra.begin(), arr_extra.end());
  std::set<int>::iterator itr;

  // ------------------------------------------------------------------------------------
  // construct mass and stiffness matrix from semi-implicit structure
  // ------------------------------------------------------------------------------------  

  auto du(u);
  for (int subIdx=0; subIdx<sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx]; 
    std::string path = std::to_string(subIdx);
    du *= 0;

    std::map<int,std::set<int>> nbrs(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
    double coef = 0.5;
    std::vector<int> nbr_tags = sequenceOfsubdomains[tag];
    F.Mass_stiff(1);
    F.set_mass_submatrix(true);
    SemiImplicitEulerStep<Functional>  eqM(&F,dt);
    eqM.setTau(0);
    Matrix subMatrix(creator);

    // mass
    Matrix subMatrix_mass(creator);
    Vector Fs_petcs_sub =  Fs_petcs[subIdx];

    // construct mass
    for ( const auto & gamma_nbr: nbrs )
    {
      Cellfltr.set_cells(cells_set);
      Cellfltr.set_tags(tags);
      int row = tag;
      int col = gamma_nbr.first;
      {
        Matrix M_sub(creator);
        F.set_row_col_subdomain(row,col);
        assembler.assemble(SemiLinearization(eqM,u,u,du), Assembler::MATRIX,assemblyThreads); 
        Matrix sub_M_ = assembler.template get<Matrix>(false);
        for (int k = 0; k < sub_M_.N(); ++k)
        {
          auto row_ = sub_M_[k];
          for (auto ca=row_.begin(); ca!=row_.end(); ++ca)
          {
            int const l = ca.index();
            M_sub[map_indices[k]][map_indices[l]] = coef*(*ca); 
          }
        }
        subMatrix_mass+=M_sub;
      }

      {
        Matrix M_sub(creator);
        F.set_row_col_subdomain(col,row);
        assembler.assemble(SemiLinearization(eqM,u,u,du), Assembler::MATRIX,assemblyThreads); 
        Matrix sub_M_ = assembler.template get<Matrix>(false);

        for (int k = 0; k < sub_M_.N(); ++k)
        {
          auto row_ = sub_M_[k];
          for (auto ca=row_.begin(); ca!=row_.end(); ++ca)
          {
            int const l = ca.index();
            M_sub[map_indices[k]][map_indices[l]] = coef*(*ca); 
          }
        }
        subMatrix_mass+=M_sub;
      }
    }

    subMatrix+=subMatrix_mass;
    // if(write_to_file) writeToMatlabPath(subMatrix_mass,Fs_petcs_sub,"mass"+path,matlab_dir, false);


   // construct stiffness
    Matrix subMatrix_stiffness(creator);
    {
      SemiImplicitEulerStep<Functional>  eqK(&F,dt);
      F.Mass_stiff(0);
      F.set_mass_submatrix(false);
      eqK.setTau(1);
    
      std::set<int> tags_target;
      tags_target.insert(tag);
      Cellfltr.set_tags(tags_target);
      Cellfltr.select_based_on_tag(true);
      assembler.template assemble<AssemblyDetail::TakeAllBlocks,CellFilter>(SemiLinearization(eqK,u,u,du),Cellfltr,Assembler::MATRIX|Assembler::RHS,assemblyThreads);  
      K_ = assembler.template get<Matrix>(false);
      K_*=(-dt);
      exctract_petsc_stiffness_blocks_moreExtracellular(sequenceOfTags, map_indices, map_II, map_GammaGamma, K_, i2Tag, subIdx,subMatrix_stiffness); 
      subMatrix+=subMatrix_stiffness;
      // if(write_to_file) writeToMatlabPath(subMatrix_stiffness,Fs_petcs_sub,"stiffness"+path,matlab_dir, false);
      Cellfltr.select_based_on_tag(false);
    }
   
    subMatrices[tag] = subMatrix;
    subMatrices_M.push_back(subMatrix_mass);
    subMatrices_K.push_back(subMatrix_stiffness);
    if(write_to_file) writeToMatlabPath(subMatrix,Fs_petcs_sub,"resultBDDC"+path,matlab_dir, false);
  }

  // -------------------------------------
  // compute the weight for each subdomain to update the rhs
  // -------------------------------------
  // - we need to go through the interfaces  
  //  - of the current and their neighbors   
  // - update the rhs for each subdomain
  // - generate the submatrices with the writeToMatlabPath again
  // -------------------------------------

  for (int subIdx=0; subIdx<sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx]; 
    std::string path = std::to_string(subIdx+1);
    Matrix subMatrix(creator);
    subMatrix = subMatrices[tag];
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
    }
    Fs_petcs[subIdx] = Fs_petcs_sub;
    weights[tag] = weights_sub;
    // std::cout<< " ===================================================== "<<std::endl;
    // if(write_to_file) writeToMatlabPath(subMatrix,Fs_petcs_sub,"resultBDDCNew"+path,matlab_dir);
  }

  return u;
}

void generate_Interror_and_Interfaces_indices(std::vector<int> sequenceOfTags,
                                              std::map<int,std::set<int>> map_II,
                                              std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                                              std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate,
                                              std::map<int, int> map_indices, 
                                              std::string matlab_dir)
{
  double precision = 16;
  for (int i = 0; i < sequenceOfTags.size(); ++i)
  {
    int tag = sequenceOfTags[i]; 
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());

    std::string Ii = matlab_dir+"/BDDC_I"+ std::to_string(i) + ".txt";
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
    std::vector<int> interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
    std::vector<int> gamma_W_Nbr(map_GammaNbr_Nbr_noDuplicate[tag].begin(), map_GammaNbr_Nbr_noDuplicate[tag].end());

    std::string gammai = matlab_dir+"/BDDC_Gamma"+ std::to_string(i) + ".txt";
    std::ofstream Gammaii(gammai.c_str());
    Gammaii.precision(precision);

    std::string gammai_w_nbr = matlab_dir+"/BDDC_Gamma_w_nbr"+ std::to_string(i) + ".txt";
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
void construct_As(std::vector<int> arr_extra, std::vector<int> sequenceOfTags, 
                  std::map<int,int> startingIndexOfTag,
                  std::map<int,std::set<int>> map_II,
                  std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                  std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate,
                  Vector rhs_kaskade,
                  std::map<int,std::vector<int>> sequanceOfsubdomains,
                  std::map<int, Vector> weights,
                  std::map<int, int> map_indices_petsc,
                  std::map<int,std::map<int,std::set<int>>> map_GammaNbr,
                  std::vector<std::set<int>> i2t, 
                  std::string matlab_dir, bool write_to_file,
                  std::map<int,Matrix> subMatrices,
                  std::vector<Matrix> subMatrices_M,
                  std::vector<Matrix> subMatrices_K,
                  std::map<int,Matrix> &subMatrices_kaskade,
                  std::vector<Matrix> &subMatrices_kaskade_M,
                  std::vector<Matrix> &subMatrices_kaskade_K,
                  std::map<int,Vector> &Fs,
                  std::map<int,std::map<int,int>> &map_kaskadeToPetscAll,
                  std::map<int,std::map<int,int>> &map_indices_kaskadeAll,
                  std::vector<std::vector<LocalDof>> &sharedDofsKaskade)
{

  std::map<int,std::vector<int>> sequanceOfsubdomainsKaskade;

  std::set<int> arr_extra_set(arr_extra.begin(), arr_extra.end()); 
  std::set<int>::iterator itr_extra; 
  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {
    int tag = sequenceOfTags[subIdx]; 
    int count = 0;

    bool extra = false;
    itr_extra =arr_extra_set.find(tag);
    if(itr_extra!=arr_extra_set.end())
      extra = true;


    Matrix subMatrix  = subMatrices[subIdx];
    Matrix subMatrix_M  = subMatrices_M[subIdx];
    Matrix subMatrix_K  = subMatrices_K[subIdx];

    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> Interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
    int size = Interior.size()+Interface.size();
    

    int counter_kasakde = 0;
    std::map<int, int> map_indices_kaskade;
    std::map<int, int> map_kaskadeToPetsc;
    {
      for (int k = 0; k < Interior.size(); ++k)
      {
        //map_indices_kaskade[counter_kasakde] = map_indices_petsc[Interior[k]];
        map_indices_kaskade[map_indices_petsc[Interior[k]]] = counter_kasakde;
        map_kaskadeToPetsc[counter_kasakde] = map_indices_petsc[Interior[k]];
        if(false) std::cout <<map_indices_petsc[Interior[k]] + 1 << ":" << counter_kasakde+1 << "\n";
        counter_kasakde++;
      }   

      for (int k = 0; k < Interface.size(); ++k)
      {
        // map_indices_kaskade[counter_kasakde] = map_indices_petsc[Interface[k]];
        map_indices_kaskade[map_indices_petsc[Interface[k]]] = counter_kasakde;
        map_kaskadeToPetsc[counter_kasakde] = map_indices_petsc[Interface[k]];
        if(false) std::cout <<map_indices_petsc[Interface[k]] + 1 << ":" << counter_kasakde+1 << "\n";
        counter_kasakde++;
      }  
    
      std::map<int,std::set<int>> Nbrs= map_GammaNbr[tag];

      for ( const auto &nbr : Nbrs)
      {
        int tag_nbr = nbr.first;
        std::set<int>::iterator itr;

        bool extra_second = false;
        itr_extra =arr_extra_set.find(tag_nbr);
        if(itr_extra!=arr_extra_set.end())
          extra_second = true;

        std::vector<int> Interior_nbr(map_II[tag_nbr].begin(), map_II[tag_nbr].end());
        std::vector<int> Interface_nbr(map_GammaGamma_noDuplicate[tag_nbr].begin(), map_GammaGamma_noDuplicate[tag_nbr].end());
        int size_nbr = Interior_nbr.size()+Interface_nbr.size();
        if(extra and extra_second)
        {
          int start = startingIndexOfTag[tag_nbr];
          
          std::vector<int> IG(size_nbr); // vector with size ints.
          std::iota (std::begin(IG), std::end(IG), start); // Fill with start, 1, ..., size.
          auto block_nbr = subMatrix(IG,IG);

          std::set<int> indices_nbr;       
          for (int k = 0; k < block_nbr.N(); ++k)
          {
            auto row  = block_nbr[k];
            for (auto ca=row.begin(); ca!=row.end(); ++ca)
            {
              int const l = ca.index();
              double val = *ca;
              if(val> 1e-16 or val < -1e-16){
                indices_nbr.insert(k+start);
              }
            }
          }

          std::set<int>::iterator it;
          for (it = indices_nbr.begin(); it != indices_nbr.end(); ++it) 
          {
            map_indices_kaskade[*it] = counter_kasakde;
            map_kaskadeToPetsc[counter_kasakde] = *it;
            if(false) std::cout  << *it + 1 << ":" << counter_kasakde+1 << "\n";
            counter_kasakde++;
          }
          if(indices_nbr.size()>0)
            size+=size_nbr;
        }
        else
        {
          int start = startingIndexOfTag[tag_nbr];
          int size_nbr = Interior_nbr.size()+Interface_nbr.size();

          for (itr = nbr.second.begin(); itr != nbr.second.end(); itr++) 
          {
            map_indices_kaskade[map_indices_petsc[*itr]] = counter_kasakde;
            map_kaskadeToPetsc[counter_kasakde] = map_indices_petsc[*itr];
            if(false) std::cout <<map_indices_petsc[*itr] + 1 << ":" << counter_kasakde+1 << "\n";
            counter_kasakde++;
          }
          size+=size_nbr;
        }
      }
      std::cout << "\n";
    }

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // creator for kaskade structure, the shrinked version
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    NumaCRSPatternCreator<> creator_kaskade(counter_kasakde,counter_kasakde,false);
    for (int k=0; k<subMatrix.N(); ++k)
    {
      auto row  = subMatrix[k];
      for (auto ca=row.begin(); ca!=row.end(); ++ca)
      {
        int const l = ca.index();
        int row_indx = map_indices_kaskade[k];
        int col_indx = map_indices_kaskade[l];

        double val = *ca;
        if(val> 1e-16 or val < -1e-16){
          creator_kaskade.addElement(row_indx,col_indx);                
        }
      }
    }

    Matrix subMatrix_kaskade_shrinked(creator_kaskade); 
    for (int k=0; k<subMatrix.N(); ++k)
    {
      auto row  = subMatrix[k];
      for (auto ca=row.begin(); ca!=row.end(); ++ca)
      {
        int const l = ca.index();
        int row_indx = map_indices_kaskade[k];
        int col_indx = map_indices_kaskade[l];

        double val = *ca;
        if(val> 1e-16 or val < -1e-16){
          subMatrix_kaskade_shrinked[row_indx][col_indx] = val;                
        }
      }
    }
    subMatrices_kaskade[tag] = subMatrix_kaskade_shrinked;
    std::cout <<"tag " << tag << " : " << counter_kasakde << std::endl;

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // save sub_matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    map_kaskadeToPetscAll[tag] = map_kaskadeToPetsc;
    map_indices_kaskadeAll[tag] = map_indices_kaskade;
    Vector Fs_subIdx(counter_kasakde);  
    { 
      for (int i = 0; i < counter_kasakde; ++i)
      {
        int index = map_kaskadeToPetsc[i]; 
        sequanceOfsubdomainsKaskade[tag].push_back(index);
        int coef = weights[tag][index];
        Fs_subIdx[i] = coef*rhs_kaskade[index];
      }
    }
    Fs[subIdx] = Fs_subIdx;
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
      //if(values.size()>1 and write_to_file){
      if(write_to_file){
        f << key << "-> ";
        for (const auto& value : values) {
          f <<"("<<std::get<0>(value) << " "<< std::get<1>(value) << ") ";
        }
        f << "\n";
      }
    }
  } 
}

// template<class Matrix>
// void construct_As(std::vector<int> sequenceOfTags, 
//                   std::map<int,int> startingIndexOfTag,
//                   std::map<int,std::set<int>> map_II,
//                   std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
//                   std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate,
//                   std::map<int,std::vector<int>> sequanceOfsubdomains,
//                   std::map<int, int> map_indices_petsc,
//                   std::vector<std::set<int>> i2t, 
//                   std::vector<Matrix> subMatrices,
//                   std::vector<Matrix> subMatrices_M,
//                   std::vector<Matrix> subMatrices_K,
//                   std::vector<Matrix> &subMatrices_kaskade,
//                   std::vector<Matrix> &subMatrices_kaskade_M,
//                   std::vector<Matrix> &subMatrices_kaskade_K)
// {

//   for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
//   {
//     int tag = sequenceOfTags[subIdx]; 
//     // std::map<int,std::set<int>> gamma_nbrs_subIdx(map_GammaNbr[tag].begin(), map_GammaNbr[tag].end());
//     std::vector<int> sequanceOfsubdomains_subIdx =  sequanceOfsubdomains[tag];

//     Matrix subMatrix  = subMatrices[subIdx];
//     Matrix subMatrix_M  = subMatrices_M[subIdx];
//     Matrix subMatrix_K  = subMatrices_K[subIdx];

//     std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
//     std::vector<int> Interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
//     std::vector<int> Interface_nbr(map_GammaNbr_Nbr_noDuplicate[tag].begin(), map_GammaNbr_Nbr_noDuplicate[tag].end());

//     if(false) std::cout <<"tag: "<< tag << std::endl;
//     int start = startingIndexOfTag[tag];
//     std::map<int, int> map_indices_kaskade;
//     {
//       int counter_kasakde = 0;
  
//       for (int k = 0; k < Interior.size(); ++k)
//       {
//         map_indices_kaskade[map_indices_petsc[Interior[k]]] = counter_kasakde;
//         if(false) std::cout << Interior[k] << ":" <<map_indices_petsc[Interior[k]] + 1 << ":" << counter_kasakde+1 << "  ";
//         counter_kasakde++;
//       }   

//       for (int k = 0; k < Interface.size(); ++k)
//       {
//         map_indices_kaskade[map_indices_petsc[Interface[k]]] = counter_kasakde;
//         if(false) std::cout << Interface[k] << ":" <<map_indices_petsc[Interface[k]] + 1 << ":" << counter_kasakde+1 << "  ";
//         counter_kasakde++;
//       }    

//       for (int k = 0; k < Interface_nbr.size(); ++k)
//       {
//         map_indices_kaskade[map_indices_petsc[Interface_nbr[k]]] = counter_kasakde;
//         if(false) std::cout << Interface_nbr[k] << ":" <<map_indices_petsc[Interface_nbr[k]] + 1 << ":" << counter_kasakde+1 << "  ";
//         counter_kasakde++;
//       }  
//       if(false) std::cout << "\n";
//     }
//     if(false) std::cout <<"============================" << std::endl;    
//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     // modified sub matrices which includes only nbr indices for kaskade format
//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     int Dofs_kaskade = Interior.size() + Interface.size() + Interface_nbr.size();   

//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     // creator for kaskade structure, the shrinked version
//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     NumaCRSPatternCreator<> creator_kaskade(Dofs_kaskade,Dofs_kaskade,false);
//     for (int k=0; k<subMatrix.N(); ++k)
//     {
//       auto row  = subMatrix[k];
//       for (auto ca=row.begin(); ca!=row.end(); ++ca)
//       {
//         int const l = ca.index();
//         int row_indx = map_indices_kaskade[k];
//         int col_indx = map_indices_kaskade[l];
//         creator_kaskade.addElement(row_indx,col_indx);  
//       }
//     }

//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     // fill the sub matrices for kaskade format
//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     Matrix subMatrix_kaskade_shrinked(creator_kaskade);  
//     Matrix subMatrix_kaskade_shrinked_M(creator_kaskade);  
//     Matrix subMatrix_kaskade_shrinked_K(creator_kaskade);  
//     {
//       // block for subIdx located at 0, 0 of submatrices    
//       // IG for subIdx
//       int size = Interior.size()+Interface.size();
//       std::vector<int> IG(size); // vector with size ints.
//       std::iota (std::begin(IG), std::end(IG), start); // Fill with start, 1, ..., size.
//       // std::cout <<"subIdx:" << subIdx << ", start: " << start << " IG[0]" << IG[0] << " IG[size-1]" << IG[size-1] << std::endl;
//       auto IGAMMA_block = subMatrix(IG,IG);
//       insertMatrixBlock(IGAMMA_block, 0, 0, subMatrix_kaskade_shrinked);

//       auto IGAMMA_block_M = subMatrix_M(IG,IG);
//       insertMatrixBlock(IGAMMA_block_M, 0, 0, subMatrix_kaskade_shrinked_M);

//       auto IGAMMA_block_K = subMatrix_K(IG,IG);
//       insertMatrixBlock(IGAMMA_block_K, 0, 0, subMatrix_kaskade_shrinked_K);

//     }

//     // rest of diagonal
//     {     

//       int size = Interior.size() + Interface.size();
//       int x = size;
//       int y = size;


//       std::vector<int> mapped_Interface_nbr(Interface_nbr.size()); // vector with size ints.
//       for (int i = 0; i < Interface_nbr.size(); ++i)
//       {
//         mapped_Interface_nbr[i] = map_indices_petsc[Interface_nbr[i]];      
//       }

//       auto GAMMAGAMMA_diagonal_block = subMatrix(mapped_Interface_nbr,mapped_Interface_nbr);
//       insertMatrixBlock(GAMMAGAMMA_diagonal_block, x, y, subMatrix_kaskade_shrinked);
    
//       auto GAMMAGAMMA_diagonal_block_M = subMatrix_M(mapped_Interface_nbr,mapped_Interface_nbr);
//       insertMatrixBlock(GAMMAGAMMA_diagonal_block_M, x, y, subMatrix_kaskade_shrinked_M);
//     }
    
//     // std::cout << "  cross blocks on column for subIdx  " <<"\n";
//     int start_x = Interior.size();
//     int start_y = Interior.size() + Interface.size();
//     // cross blocks on column for subIdx 
//     {
//       std::vector<int> mapped_Interface_nbr(Interface_nbr.size()); // vector with size ints.
//       for (int i = 0; i < Interface_nbr.size(); ++i)
//       {
//         mapped_Interface_nbr[i] = map_indices_petsc[Interface_nbr[i]];      
//       }

//       int x_nbr = Interior.size();
//       std::vector<int> G_subIdx(Interface.size()); // vector with size ints.
//       std::iota (std::begin(G_subIdx), std::end(G_subIdx), startingIndexOfTag[tag] + x_nbr); // Fill with start, 1, ..., size.

//       auto GAMMAGAMMA_cross_block = subMatrix(G_subIdx, mapped_Interface_nbr);
//       insertMatrixBlockNbr(GAMMAGAMMA_cross_block, start_x, start_y, subMatrix_kaskade_shrinked);

//       auto GAMMAGAMMA_cross_block_M = subMatrix_M(G_subIdx, mapped_Interface_nbr);
//       insertMatrixBlockNbr(GAMMAGAMMA_cross_block_M, start_x, start_y, subMatrix_kaskade_shrinked_M);
//     }
//     // std::cout<< "\n\n\n";
//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     // save sub_matrices for kaskade format
//     // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//     subMatrices_kaskade.push_back(subMatrix_kaskade_shrinked);
//     subMatrices_kaskade_M.push_back(subMatrix_kaskade_shrinked_M);
//     subMatrices_kaskade_K.push_back(subMatrix_kaskade_shrinked_K);
//   }
// }

template<class Vector>
void construct_Fs(std::vector<int> sequenceOfTags, 
                  std::map<int,int> startingIndexOfTag,
                  std::map<int,std::set<int>> map_II,
                  std::map<int,std::set<int>> map_GammaGamma_noDuplicate,
                  std::map<int,std::set<int>> map_GammaNbr_Nbr_noDuplicate,
                  Vector rhs_kaskade,
                  std::vector<std::vector<LocalDof>> sharedDofs, 
                  std::map<int, Vector> weights,
                  std::map<int, int> map_indices,
                  std::vector<Vector> &Fs)
{

  for (int subIdx = 0; subIdx < sequenceOfTags.size(); ++subIdx)
  {  
    int tag = sequenceOfTags[subIdx];
    double coef;

    // // block for subIdx located at 0, 0 of submatrices
    std::vector<int> Interior(map_II[tag].begin(), map_II[tag].end());
    std::vector<int> Interface(map_GammaGamma_noDuplicate[tag].begin(), map_GammaGamma_noDuplicate[tag].end());
    std::vector<int> Interface_nbr(map_GammaNbr_Nbr_noDuplicate[tag].begin(), map_GammaNbr_Nbr_noDuplicate[tag].end());
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // modified sub matrices which inlcludes only nbr indices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    int Dofs_kaskade = Interior.size() + Interface.size() + Interface_nbr.size(); 
    if(false){
      std::cout << "========================================================================" <<std::endl;
      std::cout << "tag "<< tag << " Dofs_kaskade: " << Dofs_kaskade <<std::endl;
      std::cout << "========================================================================" <<std::endl;      
    }  

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // fill the sub matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    Vector Fs_subIdx(Dofs_kaskade);  
    {      
      for (int i = 0; i < Interior.size(); ++i)
      {
        Fs_subIdx[i] = rhs_kaskade[Interior[i]];
      }
      if(false) std::cout <<"\n";
      for (int i = 0; i < Interface.size(); ++i)
      {
        coef = weights[tag][map_indices[Interface[i]]];
        Fs_subIdx[i+Interior.size()] = coef*rhs_kaskade[Interface[i]]; 
      }
      if(false) std::cout <<"\n";
      // rest of diagonal
      int size = Interior.size() + Interface.size();
      for (int i = 0; i < Interface_nbr.size(); ++i)
      {
        coef = weights[tag][map_indices[Interface_nbr[i]]];
        Fs_subIdx[size + i] = coef*rhs_kaskade[Interface_nbr[i]];
      }
      if(false) std::cout <<"\n==============================================" <<std::endl;
    }

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // save sub_matrices for kaskade format
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    Fs.push_back(Fs_subIdx);
  }
}

#endif