#ifndef ap_utils_read_tpp__
#define ap_utils_read_tpp__

template <typename T>
std::vector<T> read(std::istream & ifs) {
  // read all data in to a flat sequence container
  std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  if(!ifs.good()) std::cout << "ERROR :: input stream is not good " << std::endl;
  std::string line{};
  std::vector<T> v{};
  while(std::getline(ifs, line)) {
    auto first_non_white_char = line.find_first_not_of(" \t\n");
    if(first_non_white_char != std::string::npos) {
      if('#' != line[first_non_white_char]) {
        std::istringstream iss(line);
        T tmp{};
        while(iss >> tmp) {
          v.push_back(tmp);
        }
      }
    }
  }
  return v;
}


template <typename T>
std::vector<T> read(std::string const & filename) {
  // read all data in to a flat sequence container
  std::cout << "\n" __FILE__ << " --> " << __func__ << "() --> " << __LINE__ << "\n";
  std::ifstream ifs(filename, std::ios::in);
  if(!ifs.good()) std::cout << "ERROR :: input stream is not good " << std::endl;
  std::string line{};
  std::vector<T> v{};
  while(std::getline(ifs, line)) {
    auto first_non_white_char = line.find_first_not_of(" \t\n");
    if(first_non_white_char != std::string::npos) {
      if('#' != line[first_non_white_char]) {
        std::istringstream iss(line);
        T tmp{};
        while(iss >> tmp) {
          v.push_back(tmp);
        }
      }
    }
  }
  return v;
}

#endif
