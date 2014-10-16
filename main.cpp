#include <clocale>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <chrono>
#include <ctime>

#include "omp.h"

#include "boost/dynamic_bitset.hpp"

#include <stdio.h>
#include <unistd.h>

#define BOOL char

using namespace std;

template<class V>
class Matrix {
  std::vector<V> vec;
public:
  size_t dimX;
  size_t dimY;
  V& operator()(size_t i, size_t j) {
    return vec[i*dimY+j];
  }
  const V operator()(size_t i, size_t j) const {
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
  void printValToConsole() const {
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
  size_t dimX;
  size_t dimY;
public:
  inline BOOL& operator()(size_t i, size_t j) {
    return b[i*dimY+j];
  }
  inline BOOL operator()(size_t i, size_t j) const {
    return b[i*dimY+j];
  }
  MatrixB() {}
  MatrixB(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    b.resize(dimX*dimY);
  }
  size_t getDimX() { return dimX; }
  size_t getDimY() { return dimY; }

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
  bool none() const {
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++) {
        if(b[i*dimY+j]) return false;
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
    MatrixB newMat = *this;
    size_t tmp = dimX;
    dimX = dimY;
    dimY = tmp;
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        this->operator()(i,j) = newMat(j,i);
  }
};

class MatrixSmallB {
  vector<BOOL> b;
  size_t dimX;
  size_t dimY;
public:
  inline BOOL& operator()(size_t i, size_t j)       { return b[i*dimY+j];}
  inline BOOL  operator()(size_t i, size_t j) const { return b[i*dimY+j];}
  MatrixSmallB() {}
  MatrixSmallB(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    b.resize(dimX*dimY);
  }
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
};

class MatrixSmallB2 {
  boost::dynamic_bitset<> b;
  size_t dimX;
  size_t dimY;
public:
  inline boost::dynamic_bitset<>::reference operator()(size_t i, size_t j)       { return b[i*dimY+j];}
  inline bool                               operator()(size_t i, size_t j) const { return b[i*dimY+j];}
  MatrixSmallB2() {}
  MatrixSmallB2(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    b.resize(dimX*dimY);
  }
  void write(std::string path) const {
    std::ofstream FILE(path, std::ios::out | std::ofstream::binary);
    FILE.write(reinterpret_cast<const char *>(&dimX), sizeof(size_t));
    FILE.write(reinterpret_cast<const char *>(&dimY), sizeof(size_t));
    FILE << b;
  }
  void load(std::string path) {
    std::ifstream INFILE(path, std::ios::in | std::ifstream::binary);
    INFILE.read(reinterpret_cast<char *>(&dimX), sizeof(size_t));
    INFILE.read(reinterpret_cast<char *>(&dimY), sizeof(size_t));
    b.resize(dimX*dimY);
    INFILE >> b;
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
  if (!exist(ss.str())) {
    char cCurrentPath[FILENAME_MAX];

    if (!getcwd(cCurrentPath, sizeof(cCurrentPath)))
    {
      return errno;
    }

    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

    wcout << "The current working directory is " << cCurrentPath;
    wcout << "Data file does not exist." << endl;
    exit(1);
  }
  ifstream myfile(ss.str(), ios::in);
  std::string s;
  getline(myfile,s);
  wcout << s.c_str() << endl;
  size_t dimX = atoi(s.c_str());
  getline(myfile,s);
  size_t dimY = atoi(s.c_str());
  for(size_t i=0;i<dimX;i++) {
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
  for(size_t i=0;i<dimY;i++) {
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
    for(size_t line=0;line < dim[ori];line++) {
      std::stringstream ss;
      ss << inpNr << "/ori" << ori << "_line" << line << ".dat";
      if (!exist(ss.str())) {
        std::vector<uint8_t> black = M[ori]->operator[](line);
        uint8_t noWhiteBlocks = black.size() + 1;
        uint8_t noWhiteSquares = dim[1-ori] - accumulate(black.begin(),black.end(),0);
        size_t noComps = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1)+2*alg_n_k(noWhiteSquares-1,noWhiteBlocks-2)+alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
        wcout << "Creating sol O" << ori << "L" << line << ", using composition (" << (int) noWhiteBlocks << " " << (int) noWhiteSquares << ") with " << noComps << " compositions." << endl;
        std::stringstream ss2;
        ss2 << "comp_" << (size_t) noWhiteBlocks << "_" << (size_t) noWhiteSquares << ".dat";
        if (!exist(ss2.str())) fun_cr_comps(noWhiteBlocks,noWhiteSquares);
        MatrixSmallB2 S(noComps,dim[1-ori]);
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
  Matrix<MatrixSmallB2> posSol(2,maxDim);
  vector<Matrix<vector<BOOL> > > vActive;
  vActive.push_back(Matrix<vector<BOOL> >(2,maxDim));
  for(size_t ori=0;ori<2;ori++) {
    for(size_t line=0;line < dim[ori];line++) {
      std::vector<uint8_t> black = M[ori]->operator[](line);
      uint8_t noWhiteBlocks = black.size() + 1;
      uint8_t noWhiteSquares = dim[1-ori] - accumulate(black.begin(),black.end(),0);
      size_t noComps = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1)+2*alg_n_k(noWhiteSquares-1,noWhiteBlocks-2)+alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
      wcout << "Ori=" << ori << " Line=" << line << endl;
      std::stringstream ss;
      ss << inpNr << "/ori" << ori << "_line" << line << ".dat";
      posSol(ori,line).load(ss.str());
      vActive.back()(ori,line).resize(noComps,true);
    }
  }
  wcout << "Creation finished!" << endl;

  //Iteration
  wcout << "Iteration started!" << endl;

  vector<MatrixB> vSolution;
  vSolution.push_back(MatrixB(dimX,dimY));

  vector<MatrixB> vSolutionIsSet;
  vSolutionIsSet.push_back(MatrixB(dimX,dimY));

  //  vector<MatrixB> vSolutionIsSetOld;
  //  vSolutionIsSetOld.push_back(MatrixB(dim[1],dim[0]));
  //  vSolutionIsSetOld[0].set(true);
  MatrixB solutionIsSetOld(dimX,dimY);
  solutionIsSetOld.set(true);
  MatrixB solutionIsSetChange(dimX,dimY);

  vector<BOOL> allTrueInActiveColumns(maxDim,true);
  vector<BOOL> allFalseInActiveColumns(maxDim,true);

  struct guessTriple {
    size_t nActiveLines;
    size_t ori;
    size_t line;
    size_t firstActiveIdx;
  };

  vector<guessTriple> vGuessTriple;

  auto t_start = std::chrono::high_resolution_clock::now();
  while(true) {
    //Assign current correct pointers at the current incLevel;
    Matrix<vector<BOOL> >* active = &vActive.back();
    MatrixB* solution = &vSolution.back();
    MatrixB* solutionIsSet = &vSolutionIsSet.back();

    solutionIsSetChange = *solutionIsSet ^ solutionIsSetOld;
    solutionIsSetOld = *solutionIsSet;

    for(size_t ori=0;ori<2;ori++) {
      for(size_t line=0;line < dim[ori];line++) {
        MatrixSmallB2& possibleCombinations = posSol(ori,line);
        size_t npSol = active->operator ()(ori,line).size();
        size_t dimLine = dim[1-ori];
        if (solutionIsSetChange.anyInRow(line)) {
          bool allLinesDropped = true;
          for(size_t i=0 ; i < npSol ; i++) {
            if(active->operator ()(ori,line)[i]) {
              for(size_t j=0 ; j < dimLine ; j++) {
                if (solutionIsSet->operator ()(line,j) & (solution->operator ()(line,j) ^ possibleCombinations(i,j))) {
                  active->operator ()(ori,line)[i] = false;
                  break;
                } else {
                  allLinesDropped = false;
                }
              }
            }
          }
          if(allLinesDropped) {
            if (vActive.size() == 1) {
              wcout << "Error in data!" << endl;
              exit(1);
            }
            wcout << wstring(2*(vActive.size()-1),L'-') << "Invalid solution in O" << vGuessTriple.back().ori << "L" << vGuessTriple.back().line << " found and dropped, remaining " << ((vGuessTriple.back().nActiveLines) - 1) << " solutions." << endl;
            vActive.pop_back();
            vSolution.pop_back();
            vSolutionIsSet.pop_back();
            active = &vActive.back();
            solution = &vSolution.back();
            solutionIsSet = &vSolutionIsSet.back();
            active->operator ()(vGuessTriple.back().ori,vGuessTriple.back().line)[vGuessTriple.back().firstActiveIdx] = false;
            vGuessTriple.pop_back();
            solutionIsSetChange = *solutionIsSet ^ solutionIsSetOld;
          }
        }
        std::fill(allTrueInActiveColumns .begin(),allTrueInActiveColumns .end(),true);
        std::fill(allFalseInActiveColumns.begin(),allFalseInActiveColumns.end(),true);
        for(size_t i=0;i<npSol;i++)
          if(active->operator ()(ori,line)[i])
            for(size_t j=0;j<dimLine;j++)
              if(possibleCombinations(i,j)) {
                allFalseInActiveColumns[j] = false;
              } else {
                allTrueInActiveColumns[j] = false;
              }
        solution->set(line,allTrueInActiveColumns,true);
        solutionIsSet->set(line,allTrueInActiveColumns,true);
        solution->set(line,allFalseInActiveColumns,false);
        solutionIsSet->set(line,allFalseInActiveColumns,true);
      }
      solution->transpose();
      solutionIsSet->transpose();
      solutionIsSetChange.transpose();
    }
    solution->printBoolToConsole();
    // Check if we are finished
    if(solutionIsSet->all()) break;

    //If the latest step brought no news, make a guess
    if(solutionIsSetChange.none()) {
      guessTriple tmp;
      tmp.nActiveLines = 1e6;
      for(size_t ori=0;ori<2;ori++)  {
        for(size_t line=0;line<dim[ori];line++) {
          size_t nActiveLines = 0;
          for(size_t i=0 ; i < active->operator ()(ori,line).size() ; i++ )
            if(active->operator ()(ori,line)[i])
              nActiveLines++;
          if (1 < nActiveLines && nActiveLines < tmp.nActiveLines) {
            tmp.nActiveLines = nActiveLines;
            tmp.ori = ori;
            tmp.line = line;
          }
        }
      }
      vector<BOOL> activeCombinations = active->operator ()(tmp.ori, tmp.line);
      for(size_t i=0 ; i < activeCombinations.size() ; i++ )
        if(activeCombinations[i]) {
          tmp.firstActiveIdx = i;
          break;
        }
      vActive.push_back(vActive.back());
      vSolution.push_back(vSolution.back());
      vSolutionIsSet.push_back(vSolutionIsSet.back());
      std::fill(vActive.back()(tmp.ori,tmp.line).begin(), vActive.back()(tmp.ori,tmp.line).end(), false);
      vActive.back()(tmp.ori,tmp.line)[tmp.firstActiveIdx] = true;
      wcout << wstring(2*(vActive.size()-1),L'-') << "Guess correct solution in O" << tmp.ori << "L" << tmp.line << " of " << tmp.nActiveLines << " solutions there." << endl;
      vGuessTriple.push_back(tmp);
    }
  }
  auto t_end = std::chrono::high_resolution_clock::now();

  double calcTime = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
  wcout << "Iteration finished!" << endl;
  wcout << "Took " << calcTime << " microseconds." << endl;
  vSolution[0].printBoolToConsole();
  return 0;
}
