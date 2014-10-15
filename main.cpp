#include <clocale>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <chrono>
#include <ctime>

#define BOOL char

using namespace std;

template<class V>
class Matrix {
  std::vector<V> vec;
public:
  size_t dimX;
  size_t dimY;
  V& operator()(size_t i, size_t j) {
    if(i>=dimX || j>=dimY) {
      wcout << "Error in operator(): i or j >= dimX or dimY." << endl;
      exit(1);
    }
    return vec[i*dimY+j];
  }
  Matrix() {}
  Matrix(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    vec.resize(dimX*dimY);
  }
  void write(std::string path) const {
    ofstream outFILE(path, ios::out | ofstream::binary);
    outFILE.write(reinterpret_cast<const char *>(&dimX), sizeof(size_t));
    outFILE.write(reinterpret_cast<const char *>(&dimY), sizeof(size_t));
    outFILE.write(reinterpret_cast<const char *>(&vec[0]), vec.size()*sizeof(V));
  }
  void load(std::string path) {
    ifstream inFILE(path, ios::in | ifstream::binary);
    inFILE.read(reinterpret_cast<char *>(&dimX), sizeof(size_t));
    inFILE.read(reinterpret_cast<char *>(&dimY), sizeof(size_t));
    vec.resize(dimX*dimY);
    inFILE.read(reinterpret_cast<char *>(&vec[0]), vec.size()*sizeof(V));
  }
  void printValToConsole() {
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++)
        wcout << (int) this->operator ()(i,j) << "\t";
      wcout << endl;
    }
    wcout << endl;
  }
};

class MatrixB {
  vector<BOOL> b;
public:
  size_t dimX;
  size_t dimY;
  inline BOOL& operator()(size_t i, size_t j) {
    return b[i*dimY+j];
  }
  inline BOOL operator()(size_t i, size_t j) const {
    return b[i*dimY+j];
  }
  MatrixB(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    b.resize(dimX*dimY);
  }

  //operators
  MatrixB operator ^(MatrixB mat) const {
    MatrixB ans = *this;
    for(size_t i=0;i<b.size();i++)
      ans.b[i] = this->b[i] ^ mat.b[i];
    return ans;
  }
  bool operator==(MatrixB mat) const {
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        if (b[i*dimY+j] != mat.b[i*dimY+j])
          return false;
    return true;
  }

  //write-load-routines
  void write(std::string path) const {
    std::ofstream FILE(path, std::ios::out | std::ofstream::binary);
    FILE.write(reinterpret_cast<const char *>(&dimX), sizeof(size_t));
    FILE.write(reinterpret_cast<const char *>(&dimY), sizeof(size_t));
    FILE.write(reinterpret_cast<const char *>(&b[0]), sizeof(BOOL)*dimX*dimY);
  }
  void load(std::string path) {
    std::ifstream INFILE(path, std::ios::in | std::ifstream::binary);
    INFILE.read(reinterpret_cast<char *>(&dimX), sizeof(size_t));
    INFILE.read(reinterpret_cast<char *>(&dimY), sizeof(size_t));
    b.resize(dimX*dimY);
    INFILE.read(reinterpret_cast<char *>(&b[0]), sizeof(BOOL)*dimX*dimY);
  }
  void printBoolToConsole() const {
    setlocale(LC_ALL, "");
    wcout << "--" << endl;
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++)
        if(this->operator ()(i,j)) {
          wcout << L'\u2588';
        } else {
          wcout << " ";
        }
      wcout << endl;
    }
    wcout << "--" << endl;
  }

  //set-routines
  void set(bool val) {
    std::fill(b.begin(), b.end(), val);
  }
  void set(size_t i, vector<BOOL> vecJ, bool val) {
    for(size_t j=0;j<dimY;j++)
      if (vecJ[j])
        b[i*dimY+j] = val;
  }
  bool all() const {
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++) {
        if(!(b[i*dimY+j])) return false;
      }
    }
    return true;
  }
  bool anyInRow(size_t line) const {
    for(size_t j=0;j<dimY;j++)
      if (this->operator()(line,j))
        return true;
    return false;
  }
  void transpose() {
    MatrixB newMat(dimY,dimX);
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        newMat(j,i) = this->operator()(i,j);
    swap(newMat,*this);
  }
};

class PossibleCombinations {
private:
  std::vector<BOOL> vec;
  size_t nAllCombinations;
  size_t dim;
public:
  vector<BOOL> active;
  PossibleCombinations() {}
  PossibleCombinations(size_t nAllCombinations_, size_t dim_) : nAllCombinations(nAllCombinations_), dim(dim_) {
    vec.resize(nAllCombinations*dim);
    active.resize(nAllCombinations,true);
  }
  inline BOOL& operator()(size_t combination, size_t j) {
    return vec[combination*dim+j];
  }
  inline BOOL operator()(size_t combination, size_t j) const {
    return vec[combination*dim+j];
  }
  size_t getNAllCombinations() const {
    return nAllCombinations;
  }
  size_t getNActiveLines() const {
    size_t sum = 0;
    for(size_t i=0 ; i < nAllCombinations ; i++ )
      if(active[i])
        sum++;
    return sum;
  }
  size_t getDim() const {
    return dim;
  }
  void keepFirstActiveLine() {
    size_t firstActiveIdx = 0;
    for(size_t i=0 ; i < nAllCombinations ; i++ )
      if(active[i]) {
        firstActiveIdx = i;
        break;
      }
    std::fill(active.begin(), active.end(), false);
    active[firstActiveIdx] = true;
  }
  void dropFirstActiveLine() {
    size_t firstActiveIdx = 0;
    for(size_t i=0 ; i < nAllCombinations ; i++ )
      if(active[i]) {
        firstActiveIdx = i;
        break;
      }
    active[firstActiveIdx] = false;
  }
  void write(std::string path) const {
    std::ofstream FILE(path, std::ios::out | std::ofstream::binary);
    FILE.write(reinterpret_cast<const char *>(&nAllCombinations), sizeof(size_t));
    FILE.write(reinterpret_cast<const char *>(&dim), sizeof(size_t));
    FILE.write(reinterpret_cast<const char *>(&vec[0]), sizeof(BOOL)*nAllCombinations*dim);
  }
  void load(std::string path) {
    std::ifstream INFILE(path, std::ios::in | std::ifstream::binary);
    INFILE.read(reinterpret_cast<char *>(&nAllCombinations), sizeof(size_t));
    INFILE.read(reinterpret_cast<char *>(&dim), sizeof(size_t));
    active.resize(nAllCombinations,true);
    vec.resize(nAllCombinations*dim);
    INFILE.read(reinterpret_cast<char *>(&vec[0]), sizeof(BOOL)*nAllCombinations*dim);
  }
  void printBoolToConsole() const {
    setlocale(LC_ALL, "");
    wcout << "--" << endl;
    for(size_t i=0;i<nAllCombinations;i++) {
      for(size_t j=0;j<dim;j++)
        if(vec[i*dim+j]) {
          wcout << L'\u2588';
        } else {
          wcout << '.';
        }
      wcout << endl;
    }
    wcout << "--" << endl;
  }
};

size_t alg_n_k(size_t n_, size_t k_) {
  double n = n_;
  double k = k_;

  if (k > n) return 0;
  if (k > n / 2) {
    k  = n  - k;
    k_ = n_ - k_;
  }
  if (k_==0) return 1;
  double ans = 1.;
  size_t count = 0;
  for(size_t i=n_-k_+1;i<=n_;i++) {
    count++;
    ans *= (double)i/count;
  }
  return (size_t) ans;
}

void fun_cr_comps(size_t noWhiteBlocks, size_t noWhiteSquares) {
  size_t noComp1 = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1);
  size_t noComp2 = alg_n_k(noWhiteSquares-1,noWhiteBlocks-2);
  size_t noComp3 = alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
  size_t noComps = noComp1 + 2*noComp2 + noComp3;
  wcout << "Starting (" << noWhiteBlocks << " " << noWhiteSquares << ") composition, this will give " << noComps << " compositions." << endl;

  size_t i = 0;
  Matrix<uint8_t> A(noComps,noWhiteBlocks);

  if (noComp1 > 0) {
    int n = noWhiteSquares;
    int k = noWhiteBlocks;
    std::vector<int> a(k,0);
    a[0] = n-k;
    int t = n - k;
    int h = 0;
    for(size_t j=0;j<noWhiteBlocks;j++)
      A(i,j) = a[j]+1;
    i++;
    while(i<noComp1) {
      if (1 < t) h = 0;
      h++;
      t = a[h-1];
      a[h-1] = 0;
      a[0] = t - 1;
      a[h]++;
      for(size_t j=0;j<noWhiteBlocks;j++)
        A(i,j) = a[j]+1;
      i++;
    }
  }
  if (noComp2 > 0) {
    int n = noWhiteSquares;
    int k = noWhiteBlocks-1;
    std::vector<int> a(k,0);
    a[0] = n-k;
    int t = n - k;
    int h = 0;
    for(size_t j=0;j<noWhiteBlocks-1;j++)
      A(i,j) = a[j]+1;
    A(i,noWhiteBlocks-1) = 0;
    i++;
    A(i,0) = 0;
    for(size_t j=1;j<noWhiteBlocks;j++)
      A(i,j) = a[j-1]+1;
    i++;
    while(i<noComp1+2*noComp2) {
      if (1 < t) h = 0;
      h++;
      t = a[h-1];
      a[h-1] = 0;
      a[0] = t - 1;
      a[h] += 1;
      for(size_t j=0;j<noWhiteBlocks-1;j++)
        A(i,j) = a[j]+1;
      A(i,noWhiteBlocks-1) = 0;
      i++;
      A(i,0) = 0;
      for(size_t j=1;j<noWhiteBlocks;j++)
        A(i,j) = a[j-1]+1;
      i++;
    }
  }
  if (noComp3 > 0) {
    int n = noWhiteSquares;
    int k = noWhiteBlocks-2;
    std::vector<int> a(k,0);
    a[0] = n-k;
    int t = n - k;
    int h = 0;
    A(i,0) = 0;
    for(size_t j=1;j<noWhiteBlocks-1;j++)
      A(i,j) = a[j-1]+1;
    A(i,noWhiteBlocks-1) = 0;
    i++;
    while(i<noComp1+2*noComp2+noComp3) {
      if (1 < t) h = 0;
      h++;
      t = a[h-1];
      a[h-1] = 0;
      a[0] = t - 1;
      a[h] += 1;
      A(i,0) = 0;
      for(size_t j=1;j<noWhiteBlocks-1;j++)
        A(i,j) = a[j-1]+1;
      A(i,noWhiteBlocks-1) = 0;
      i++;
    }
  }
  std::stringstream ss;
  ss << "comp_" << noWhiteBlocks << "_" << noWhiteSquares << ".dat";
  A.write(ss.str());
  wcout << "Finished (" << noWhiteBlocks << " " << noWhiteSquares << ") composition." << endl;
}

bool exist (const std::string& name) {
  ifstream f(name.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }
}

MatrixB fun_logics_rec(MatrixB sol, MatrixB solSet, Matrix<PossibleCombinations>& posSol, bool& err, size_t incLevel) {
  vector<size_t> dim(2);dim[0] = sol.dimX;dim[1] = sol.dimY;

  MatrixB solSet_old(dim[1],dim[0]);solSet_old.set(true);

  vector<BOOL> allTrueInActiveColumns(max(dim[0],dim[1]),true);
  vector<BOOL> allFalseInActiveColumns(max(dim[0],dim[1]),true);

  while(true) {
    MatrixB solSetChange = solSet ^ solSet_old;
    solSet_old = solSet;
    MatrixB sol_old = sol;
    for(size_t ori=0;ori<2;ori++) {
      for(size_t line=0;line < dim[ori];line++) {
        PossibleCombinations& possibleCombinations = posSol(ori,line);
        size_t npSol = possibleCombinations.getNAllCombinations();
        size_t dimLine = dim[1-ori];
        if (solSetChange.anyInRow(line)) {
          bool allLinesDropped = true;
          for(size_t i=0 ; i < npSol ; i++) {
            if(possibleCombinations.active[i]) {
              for(size_t j=0 ; j < dimLine ; j++) {
                if (solSet(line,j) & (sol(line,j) ^ possibleCombinations(i,j))) {
                  possibleCombinations.active[i] = false;
                  break;
                } else {
                  allLinesDropped = false;
                }
              }
            }
          }
          if(allLinesDropped) {
            err = true;
            return sol;
          }
        }
        std::fill(allTrueInActiveColumns .begin(),allTrueInActiveColumns .end(),true);
        std::fill(allFalseInActiveColumns.begin(),allFalseInActiveColumns.end(),true);
        for(size_t i=0;i<npSol;i++)
          if(possibleCombinations.active[i])
            for(size_t j=0;j<dimLine;j++)
              if(possibleCombinations(i,j)) {
                allFalseInActiveColumns[j] = false;
              } else {
                allTrueInActiveColumns[j] = false;
              }
        sol.set(line,allTrueInActiveColumns,true);
        solSet.set(line,allTrueInActiveColumns,true);
        sol.set(line,allFalseInActiveColumns,false);
        solSet.set(line,allFalseInActiveColumns,true);
      }
      sol.transpose();
      solSet.transpose();
      solSetChange.transpose();
    }
    // Check if we are finished
    if(solSet.all()) {
      err = false;
      return sol;
    }

    //If the latest step brought no news, make a gues (=inception level)
    if(sol_old == sol) {
      size_t min_s = 1e6;
      size_t minOri = 2;
      size_t minLine = 100;
      for(size_t ori=0;ori<2;ori++)  {
        for(size_t line=0;line<dim[ori];line++) {
          size_t nActiveLines = posSol(ori,line).getNActiveLines();
          if (1 < nActiveLines && nActiveLines < min_s) {
            min_s = nActiveLines;
            minOri = ori;
            minLine = line;
          }
        }
      }
      Matrix<PossibleCombinations> thPosSol = posSol;
      thPosSol(minOri,minLine).keepFirstActiveLine();
      wcout << wstring(2*incLevel,L'-') << "Guess correct solution in O" << minOri << "L" << minLine << " of " << min_s << " solutions there." << endl;
      MatrixB thSol = fun_logics_rec(sol,solSet,thPosSol,err,incLevel+1);
      if(!err) {
        wcout << wstring(2*incLevel,L'-') << "Returned from inception level " << incLevel << ". Solved!" << endl;
        return thSol;
      } else {
        wcout << wstring(2*incLevel,L'-') << "Invalid solution in O" << minOri << "L" << minLine << " found and dropped, remaining " << min_s-1 << " solutions." << endl;
        posSol(minOri,minLine).dropFirstActiveLine();
      }
    }
  }
}

int main(int argc, const char* argv[]) {
  setlocale(LC_ALL, "");
  if (argc==1) {
    wcout << "Error: No input number provided! Exiting..." << endl;
    exit(1);
  }
  size_t inpNr = atoi(argv[1]);
  wcout << "Solving input number: " << inpNr << endl;
  std::vector<std::vector<uint8_t> > H;
  std::vector<std::vector<uint8_t> > V;
  std::vector<uint8_t> tmp;
  stringstream ss;
  ss << "data" << inpNr;
  ifstream myfile(ss.str(), ios::in);
  std::string s;
  getline(myfile,s);
  size_t dimY = atoi(s.c_str());
  getline(myfile,s);
  size_t dimX = atoi(s.c_str());
  for(size_t i=0;i<dimY;i++) {
    getline(myfile,s);
    std::stringstream ss(s);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    tmp.clear();
    for(auto & i : vstrings)
      tmp.push_back(atoi(i.c_str()));
    H.push_back(tmp);
  }
  for(size_t i=0;i<dimX;i++) {
    getline(myfile,s);
    std::stringstream ss(s);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    tmp.clear();
    for(auto & i : vstrings)
      tmp.push_back(atoi(i.c_str()));
    V.push_back(tmp);
  }

  //Creation
  wcout << "Creation started!" << endl;
  std::vector<size_t> dim(2);
  dim[0] = dimX;
  dim[1] = dimY;
  std::vector<std::vector<std::vector<uint8_t> >* > M;
  M.push_back(&H);
  M.push_back(&V);
  size_t maxDim = max(dimX, dimY);
  for(size_t ori=0;ori<2;ori++) {
    for(size_t line=0;line < dim[1-ori];line++) {
      std::stringstream ss;
      ss << inpNr << "/ori" << ori << "_line" << line << ".dat";
      if (!exist(ss.str())) {
        std::vector<uint8_t> black = M[ori]->operator[](line);
        uint8_t noWhiteBlocks = black.size() + 1;
        uint8_t noWhiteSquares = dim[ori] - accumulate(black.begin(),black.end(),0);
        size_t noComps = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1)+2*alg_n_k(noWhiteSquares-1,noWhiteBlocks-2)+alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
        wcout << "Creating sol O" << ori << "L" << line << ", using composition (" << (int) noWhiteBlocks << " " << (int) noWhiteSquares << ") with " << noComps << " compositions." << endl;
        std::stringstream ss2;
        ss2 << "comp_" << (size_t) noWhiteBlocks << "_" << (size_t) noWhiteSquares << ".dat";
        if (!exist(ss2.str())) fun_cr_comps(noWhiteBlocks,noWhiteSquares);
        PossibleCombinations S(noComps,dim[ori]);
        Matrix<uint8_t> A(noComps,noWhiteBlocks);
        A.load(ss2.str());
        for(size_t k=0;k<noComps;k++) {
          size_t pos = 0;
          for(size_t l=0;l<black.size();l++) {
            pos += A(k,l);
            for (size_t m = 0 ; m < black[l] ; m++) {
              S(k,pos++) = true;
            }
          }
        }
        S.write(ss.str());
        wcout << "Finished sol O" << ori << "L" << line << ", using composition (" << (int)  noWhiteBlocks << " " << (int)  noWhiteSquares << ") with " << noComps << " compositions." << endl;
      }
    }
  }
  MatrixB sol(dim[1],dim[0]);
  MatrixB solSet(dim[1],dim[0]);
  Matrix<PossibleCombinations> posSol(2,maxDim);
  Matrix<vector<BOOL> > active(2,maxDim);
  for(size_t ori=0;ori<2;ori++) {
    for(size_t line=0;line < dim[1-ori];line++) {
      std::stringstream ss;
      ss << inpNr << "/ori" << ori << "_line" << line << ".dat";
      posSol(ori,line).load(ss.str());
      active(ori,line).resize(posSol(ori,line).getNAllCombinations(),true);
    }
  }
  wcout << "Creation finished!" << endl;

  //Iteration
  wcout << "Iteration started!" << endl;
  bool err = false;
  auto t_start = std::chrono::high_resolution_clock::now();
  sol = fun_logics_rec(sol,solSet,posSol,err,0);
  auto t_end = std::chrono::high_resolution_clock::now();
  double calcTime = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
  if (err) {
    wcout << "Returned error :-(" << endl;
    return(1);
  } else {
    wcout << "Iteration finished!" << endl;
    wcout << "Took " << calcTime << " milliseconds." << endl;
    sol.printBoolToConsole();
    return 0;
  }
}
