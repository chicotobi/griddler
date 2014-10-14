#include <clocale>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <sstream>
#include <chrono>
#include <ctime>

#include "boost/dynamic_bitset.hpp"

#include "omp.h"

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
    V v;
    vec.clear();
    while( inFILE.read(reinterpret_cast<char *>(&v), sizeof(V)))
      vec.push_back(v);
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
  boost::dynamic_bitset<> b;
public:
  size_t dimX;
  size_t dimY;
  boost::dynamic_bitset<>::reference operator()(size_t i, size_t j) {
    if(i>=dimX || j>=dimY) {
      wcout << "Error in operator(): i or j >= dimX or dimY." << endl;
      exit(1);
    }
    return b[i*dimY+j];
  }
  bool operator()(size_t i, size_t j) const {
    if(i>=dimX || j>=dimY) {
      wcout << "Error in operator(): i or j >= dimX or dimY." << endl;
      exit(1);
    }
    return b[i*dimY+j];
  }
  MatrixB() {}
  MatrixB(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    b.resize(dimX*dimY);
  }

  //operators
  MatrixB operator &(MatrixB mat) const {
    MatrixB ans = *this;
    for(size_t i=0;i<b.size();i++)
      ans.b[i] = this->b[i] & mat.b[i];
    return ans;
  }
  MatrixB operator ^(MatrixB mat) const {
    MatrixB ans = *this;
    for(size_t i=0;i<b.size();i++)
      ans.b[i] = this->b[i] ^ mat.b[i];
    return ans;
  }
  MatrixB operator!() const {
    MatrixB ans = *this;
    for(size_t i=0;i<b.size();i++)
      ans.b[i] = !this->b[i];
    return ans;
  }
  bool operator==(MatrixB mat) const {
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        if (this->operator()(i,j) != mat(i,j)) return false;
    return true;
  }

  //write-load-routines
  void write(std::string path) const {
    ofstream outFILE(path, ios::out | ofstream::binary);
    outFILE.write(reinterpret_cast<const char *>(&dimX), sizeof(size_t));
    outFILE.write(reinterpret_cast<const char *>(&dimY), sizeof(size_t));
    outFILE << b;
  }
  void load(std::string path) {
    ifstream inFILE(path, ios::in | ifstream::binary);
    inFILE.read(reinterpret_cast<char *>(&dimX), sizeof(size_t));
    inFILE.read(reinterpret_cast<char *>(&dimY), sizeof(size_t));
    b.resize(dimX*dimY);
    inFILE >> b;
  }
  void printBoolToConsole() const {
    setlocale(LC_ALL, "");
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++)
        if(this->operator ()(i,j)) {
          wcout << L'\u2588';
        } else {
          wcout << " ";
        }
      wcout << endl;
    }
    wcout << endl;
  }

  //set-routines
  void setFalse() {
    for(size_t i=0;i<b.size();i++)
      b[i] = false;
  }
  void setTrue() {
    for(size_t i=0;i<b.size();i++)
      b[i] = true;
  }
  void setTrue(size_t i, MatrixB vecJ) {
    for(size_t j=0;j<dimY;j++)
      if (vecJ(0,j))
        this->operator ()(i,j) = true;
  }
  void setFalse(size_t i, MatrixB vecJ) {
    for(size_t j=0;j<dimY;j++)
      if (vecJ(0,j))
        this->operator ()(i,j) = false;
  }

  //any-all-routines
  bool any() const {
    for(size_t i=0;i<b.size();i++)
      if (b[i]) return true;
    return false;
  }
  bool all() const {
    for(size_t i=0;i<b.size();i++)
      if (!b[i]) return false;
    return true;
  }
  MatrixB allInColumns() const {
    MatrixB ans(1,dimY);
    ans.setTrue();
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        if(!((bool)this->operator ()(i,j))) ans(0,j) = false;
    return ans;
  }
  MatrixB allInRows() const {
    MatrixB ans(dimX,1);
    ans.setTrue();
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        if(!((bool)this->operator ()(i,j))) ans(i,0) = false;
    return ans;
  }
  MatrixB anyInColumns() const {
    MatrixB ans(1,dimY);
    ans.setFalse();
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        if(this->operator ()(i,j)) ans(0,j) = true;
    return ans;
  }
  MatrixB anyInRows() const {
    MatrixB ans(dimX,1);
    ans.setFalse();
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        if(this->operator ()(i,j)) ans(i,0) = true;
    return ans;
  }
  bool anyInRow(size_t line) const {
    for(size_t j=0;j<dimY;j++)
      if (this->operator()(line,j))
        return true;
    return false;
  }

  //complex logical routines
  void dropRows(MatrixB idxDrop) {
    //Memory check
    size_t lines = dimX;
    for(size_t i=0;i<dimX;i++)
      if(idxDrop(i,0))
        lines--;
    MatrixB newMat(lines,dimY);
    lines = 0;
    for(size_t i=0;i<dimX;i++)
      if(!(idxDrop(i,0))) {
        for(size_t j=0;j<dimY;j++)
          newMat(lines,j) = this->operator()(i,j);
        lines++;
      }
    swap(newMat,*this);
  }
  void keepFirstLine() {
    MatrixB newMat(1,dimY);
    for(size_t j=0;j<dimY;j++)
      newMat(0,j) = this->operator()(0,j);
    swap(newMat,*this);
  }
  void dropFirstLine() {
    MatrixB newMat(dimX-1,dimY);
    for(size_t i=1;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        newMat(i-1,j) = this->operator()(i,j);
    swap(newMat,*this);
  }
  void transpose() {
    MatrixB newMat(dimY,dimX);
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        newMat(j,i) = this->operator()(i,j);
    swap(newMat,*this);
  }
};

class Line {
public:
  MatrixB possibleCombinations;
  MatrixB active;
  Line() {}
  boost::dynamic_bitset<>::reference operator()(size_t i, size_t j) {
    if (active(i,0)) {
      return possibleLines(i,j);
    } else {
      wcout << "Not active!" << endl;
      exit(1);
    }
  }
  bool operator()(size_t i, size_t j) const {
    if (active(i,0)) {
      return possibleLines(i,j);
    } else {
      wcout << "Not active!" << endl;
      exit(1);
    }
  }
  size_t getNumberOfLines() {
    return possibleLines.dimX;
  }
  size_t getNumberOfActiveLines() {
    size_t ans = 0;
    for(size_t i=0;i<active.dimY;i++)
      if(active(0,i))
        ans++;
    return ans;
  }
  MatrixB allTrueInActiveColumns() const {
    MatrixB ans(1,possibleLines.dimY);
    ans.setTrue();
    for(size_t i=0;i<possibleLines.dimX;i++)
      if(active(i,0))
        for(size_t j=0;j<possibleLines.dimY;j++)
          if(!(possibleLines(i,j))) ans(0,j) = false;
    return ans;
  }
  MatrixB allFalseInActiveColumns() const {
    MatrixB ans(1,possibleLines.dimY);
    ans.setTrue();
    for(size_t i=0;i<possibleLines.dimX;i++)
      if(active(i,0))
        for(size_t j=0;j<possibleLines.dimY;j++)
          if(possibleLines(i,j)) ans(0,j) = false;
    return ans;
  }
  void keepFirstActiveLine() {
    size_t firstActiveIdx = 0;
    for(size_t i=0;i<active.dimY;i++)
      if(active(0,i)) {
        firstActiveIdx = i;
        break;
      }
    active.setFalse();
    active(0,firstActiveIdx) = true;
  }
  void dropFirstActiveLine() {
    for(size_t i=0;i<active.dimY;i++)
      if(active(0,i)) {
        active(0,i) = false;
        break;
      }
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

inline bool exist (const std::string& name) {
  ifstream f(name.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }
}

MatrixB fun_logics_rec(MatrixB sol, MatrixB solSet, Matrix<Line>& posSol, bool& err, size_t incLevel) {
  vector<size_t> dim(2);dim[0] = sol.dimX;dim[1] = sol.dimY;

  MatrixB solSet_old(dim[1],dim[0]);solSet_old.setTrue();

  while(true) {
    MatrixB solSetChange = solSet & (!solSet_old);
    solSet_old = solSet;
    MatrixB sol_old = sol;
    for(size_t ori=0;ori<2;ori++) {
      for(size_t line=0;line < dim[ori];line++) {
        if (solSetChange.anyInRow(line)) {
          err = true;
          Line* posLines = &(posSol(ori,line));
          size_t npSol = posLines.getNumberOfLines();
          size_t dimLine = posLines.possibleLines.dimY;
          //Define idxDrop vector, init false
          MatrixB idxDrop(npSol,1);
          //Loop over all compositions
          for(size_t i=0;i<npSol;i++) {
            if(posSol(ori,line).active(i,0)) {
              for(size_t j=0;j < dimLine;j++) {
                if (solSet(line,j) & (sol(line,j) ^ posSol(ori,line).possibleLines(i,j))) {
                  posSol(ori,line).active(i,0) = false;
                  break;
                }
              }
              if(!(idxDrop(i,0))) {
                err = false;
              }
            }
          }
          if(err) {
            return sol;
          }
        }
        MatrixB v1 = posSol(ori,line).allTrueInActiveColumns();
        sol.setTrue(line,v1);
        solSet.setTrue(line,v1);
        MatrixB v2 = posSol(ori,line).allTrueInActiveColumns();
        sol.setFalse(line,v2);
        solSet.setTrue(line,v2);
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
      size_t minOri = 0;
      size_t minLine = 0;
      for(size_t ori=0;ori<2;ori++)  {
        for(size_t line=0;line<dim[ori];line++) {
          if (1 < posSol(ori,line).getNumberOfActiveLines() && posSol(ori,line).getNumberOfActiveLines() < min_s) {
            min_s = posSol(ori,line).getNumberOfActiveLines();
            minOri = ori;
            minLine = line;
          }
        }
      }
      Matrix<Line> thPosSol = posSol;
      thPosSol(minOri,minLine).keepFirstActiveLine();
      wstringstream ind;for(size_t i=0;i<incLevel;i++) ind << "--";
      wcout << ind.str() << "Guess correct solution in O" << minOri << "L" << minLine << " of " << min_s << " solutions there." << endl;
      MatrixB thSol = fun_logics_rec(sol,solSet,thPosSol,err,incLevel+1);
      if(!err) {
        wcout << ind.str() << "Returned from inception level " << incLevel << ". Solved!" << endl;
        return thSol;
      } else {
        wcout << ind.str() << "Invalid solution in O" << minOri << "L" << minLine << " found and dropped, remaining " << min_s-1 << " solutions." << endl;
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
        MatrixB S(noComps,dim[ori]);
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
  Matrix<Line> posSol(2,maxDim);
  for(size_t ori=0;ori<2;ori++) {
    for(size_t line=0;line < dim[1-ori];line++) {
      std::vector<uint8_t> black = M[ori]->operator[](line);
      uint8_t noWhiteBlocks = black.size() + 1;
      uint8_t noWhiteSquares = dim[ori] - accumulate(black.begin(),black.end(),0);
      std::stringstream ss;
      ss << inpNr << "/ori" << ori << "_line" << line << ".dat";
      size_t noComps = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1)+2*alg_n_k(noWhiteSquares-1,noWhiteBlocks-2)+alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
      MatrixB S(noComps,dim[ori]);
      S.load(ss.str());
      posSol(ori,line).possibleLines = S;
      posSol(ori,line).active = MatrixB(S.dimX,1);
      posSol(ori,line).active.setTrue();
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
  wcout << "Iteration finished!" << endl;
  wcout << "Took " << calcTime << " milliseconds." << endl;

  sol.printBoolToConsole();

  return 0;
}
