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