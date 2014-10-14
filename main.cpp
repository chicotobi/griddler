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
      cout << "Error in operator(): i or j >= dimX or dimY." << endl;
      exit(1);
    }
    return vec[i*dimY+j];
  }
  Matrix() {}
  Matrix(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    vec.resize(dimX*dimY);
  }
  Matrix operator &(Matrix mat) const {
    Matrix ans = *this;
    for(size_t i=0;i<vec.size();i++)
      ans.vec[i] = ((bool) this->vec[i]) & ((bool) mat.vec[i]);
    return ans;
  }
  Matrix operator!() const {
    Matrix ans = *this;
    for(size_t i=0;i<vec.size();i++)
      ans.vec[i] = (!(bool) this->vec[i]);
    return ans;
  }
  void transpose() {
    Matrix newMat(dimY,dimX);
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        newMat(j,i) = this->operator()(i,j);
    swap(newMat,*this);
  }
  void setFalse() {
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        this->operator ()(i,j) = 0;
  }
  void setTrue() {
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        this->operator ()(i,j) = 1;
  }
  void setTrue(size_t i, Matrix<BOOL> vecJ) {
    if(i>=dimX) {
      cout << "Error in setTrue: i>dimX" << endl;
      exit(1);
    }
    for(size_t j=0;j<dimY;j++)
      if (vecJ(0,j))
        this->operator ()(i,j) = true;
  }
  void setFalse(size_t i, Matrix<BOOL> vecJ) {
    if(i>=dimX) {
      cout << "Error in setFalse: i>dimX" << endl;
      exit(1);
    }
    for(size_t j=0;j<dimY;j++)
      if (vecJ(0,j)) this->operator ()(i,j) = false;
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
  bool any() {
    for(size_t i=0;i<vec.size();i++)
      if (vec[i]) return true;
    return false;
  }
  bool all() {
    for(size_t i=0;i<vec.size();i++)
      if (!((bool)vec[i])) return false;
    return true;
  }
  Matrix<BOOL> any(size_t ori) {
    if (ori == 1) {
      Matrix<BOOL> ans(dimX,1);
      ans.setFalse();
      for(size_t i=0;i<dimX;i++)
        for(size_t j=0;j<dimY;j++)
          if(this->operator ()(i,j)) ans(i,0) = true;
      return ans;
    } else if (ori == 0) {
      Matrix<BOOL> ans(1,dimY);
      ans.setFalse();
      for(size_t i=0;i<dimX;i++)
        for(size_t j=0;j<dimY;j++)
          if(this->operator ()(i,j)) ans(0,j) = true;
      return ans;
    }
    return Matrix<BOOL>();
  }
  Matrix<BOOL> all(size_t ori) {
    if (ori == 1) {
      Matrix<BOOL> ans(dimX,1);
      ans.setTrue();
      for(size_t i=0;i<dimX;i++)
        for(size_t j=0;j<dimY;j++)
          if(!((bool)this->operator ()(i,j))) ans(i,0) = false;
      return ans;
    } else if (ori == 0) {
      Matrix<BOOL> ans(1,dimY);
      ans.setTrue();
      for(size_t i=0;i<dimX;i++)
        for(size_t j=0;j<dimY;j++)
          if(!((bool)this->operator ()(i,j))) ans(0,j) = false;
      return ans;
    }
    return Matrix<BOOL>();
  }

  bool anyInRow(size_t line) {
    for(size_t j=0;j<dimY;j++)
      if (this->operator()(line,j))
        return true;
    return false;
  }
  std::vector<V> getRow(size_t line) {
    std::vector<V> ans;
    for(size_t j=0;j<dimY;j++)
      ans.push_back(this->operator ()(line,j));
    return ans;
  }
  void keepRows(Matrix<BOOL> idxKeep) {
    //Memory check
    size_t lines = 0;
    for(size_t i=0;i<dimX;i++)
      if(idxKeep(i,0))
        lines++;
    Matrix<BOOL> newMat(lines,dimY);
    lines = 0;
    for(size_t i=0;i<dimX;i++)
      if(idxKeep(i,0)) {
        for(size_t j=0;j<dimY;j++)
          newMat(lines,j) = this->operator()(i,j);
        lines++;
      }
    swap(newMat,*this);
  }
  void keepFirstLine() {
    Matrix<BOOL> newMat(1,dimY);
    for(size_t j=0;j<dimY;j++)
      newMat(0,j) = this->operator()(0,j);
    swap(newMat,*this);
  }
  void dropFirstLine() {
    Matrix<BOOL> newMat(dimX-1,dimY);
    for(size_t i=1;i<dimX;i++)
      for(size_t j=0;j<dimY;j++)
        newMat(i-1,j) = this->operator()(i,j);
    swap(newMat,*this);
  }
  void printBoolToConsole() {
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++)
        if(this->operator ()(i,j)) {
          //wcout << L'\u25A1';
          cout << "x";
        } else {
          //wcout << L'\u25A0';
          cout << ".";
        }
      wcout << endl;
    }
    wcout << endl;
    wcout << endl;
  }

  void printValToConsole() {
    for(size_t i=0;i<dimX;i++) {
      for(size_t j=0;j<dimY;j++)
        cout << (int) this->operator ()(i,j) << "\t";
      cout << endl;
    }
    cout << endl;
  }
  bool operator==(Matrix<BOOL> mat) {
    for(size_t i=0;i<dimX;i++)
      for(size_t j=0;j<dimY;j++) {
        bool bool1 = (bool) this->operator()(i,j);
        bool bool2 = (bool) mat(i,j);
        if (bool1 != bool2) return false;
      }
    return true;
  }
};

class MatrixB : public boost::dynamic_bitset<> {
public:
  size_t dimX;
  size_t dimY;
  boost::dynamic_bitset<>::reference operator()(size_t i, size_t j) {
    if(i>=dimX || j>=dimY) {
      cout << "Error in operator(): i or j >= dimX or dimY." << endl;
      exit(1);
    }
    return this->operator [](i*dimY+j);
  }
  MatrixB() {}
  MatrixB(size_t dimX_, size_t dimY_) : dimX(dimX_), dimY(dimY_) {
    this->resize(dimX*dimY);
  }
  void write(std::string path) const {
    ofstream outFILE(path, ios::out | ofstream::binary);
    outFILE.write(reinterpret_cast<const char *>(&dimX), sizeof(size_t));
    outFILE.write(reinterpret_cast<const char *>(&dimY), sizeof(size_t));
    outFILE << *this;
    //outFILE.write(reinterpret_cast<const char *>(&this->operator [](0), this->size()*sizeof(V));
  }
  void load(std::string path) {
    ifstream inFILE(path, ios::in | ifstream::binary);
    inFILE.read(reinterpret_cast<char *>(&dimX), sizeof(size_t));
    inFILE.read(reinterpret_cast<char *>(&dimY), sizeof(size_t));
    this->resize(dimX*dimY);
    inFILE >> *this;
    //V v;
    //vec.clear();
    //while( inFILE.read(reinterpret_cast<char *>(&v), sizeof(V)))
    //  vec.push_back(v);
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
  cout << "Starting (" << noWhiteBlocks << " " << noWhiteSquares << ") composition, this will give " << noComps << " compositions." << endl;
  size_t i = 0;
  Matrix<char> A(noComps,noWhiteBlocks);
  A.setFalse();
  if (noComp1 > 0) {
    int n = noWhiteSquares;
    int k = noWhiteBlocks;
    std::vector<int> a(k,0);
    a[0] = n-k;
    int t = n - k;
    int h = 0;
    for(int j=0;j<k;j++)
      A(i,j) = a[j]+1;
    i++;
    while(i<noComp1) {
      if (1 < t) h = 0;
      h++;
      t = a[h-1];
      a[h-1] = 0;
      a[0] = t - 1;
      a[h]++;
      for(int j=0;j<k;j++)
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
    for(int j=0;j<k;j++)
      A(i,j) = a[j]+1;
    i++;
    for(int j=0;j<k;j++)
      A(i,j+1) = a[j]+1;
    i++;
    while(i<noComp1+2*noComp2) {
      if (1 < t) h = 0;
      h++;
      t = a[h-1];
      a[h-1] = 0;
      a[0] = t - 1;
      a[h] += 1;
      for(int j=0;j<k;j++)
        A(i,j) = a[j]+1;
      i++;
      for(int j=0;j<k;j++)
        A(i,j+1) = a[j]+1;
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
    for(int j=0;j<k;j++)
      A(i,j+1) = a[j]+1;
    i++;
    while(i<noComp1+2*noComp2+noComp3) {
      if (1 < t) h = 0;
      h++;
      t = a[h-1];
      a[h-1] = 0;
      a[0] = t - 1;
      a[h] += 1;
      for(int j=0;j<k;j++)
        A(i,j+1) = a[j]+1;
      i++;
    }
  }
  cout << "Finished (" << noWhiteBlocks << " " << noWhiteSquares << ") composition." << endl;
  std::stringstream ss;
  ss << "comp_" << noWhiteBlocks << "_" << noWhiteSquares << ".dat";
  A.write(ss.str());
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

void fun_logics_rec(Matrix<char>& sol, Matrix<char>& solSet, size_t incLevel, Matrix<Matrix<BOOL> >& posSol, bool& err) {

  err = false;

  vector<size_t> dim(2);
  dim[0] = sol.dimX;
  dim[1] = sol.dimY;

  Matrix<BOOL> solSet_old(dim[1],dim[0]);
  solSet_old.setTrue();

  stringstream ind;
  for(size_t i=0;i<incLevel;i++)
    ind << "--";

  while(true) {
    Matrix<BOOL> solSetChange = solSet & (!solSet_old);
    solSet_old = solSet;
    Matrix<BOOL> sol_old = sol;
    for(size_t ori=0;ori<2;ori++) {
      for(size_t line=0;line < dim[ori];line++) {
        if (solSetChange.anyInRow(line)) {
          size_t npSol = posSol(ori,line).dimX;
          Matrix<BOOL> mat = Matrix<BOOL>(npSol,dim[1-ori]);
          mat.setFalse();
          Matrix<BOOL> posSolOriLine = posSol(ori,line);
          for(size_t j=0;j < dim[1-ori];j++) {
            if(solSet(line,j)) {
              bool solLineJ = sol(line,j);
#pragma omp parallel for
              for(size_t i=0;i < npSol ; i++) {
                mat(i,j) = solLineJ^posSolOriLine(i,j);
              }
            }
          }
          Matrix<BOOL> idxKeep = !(mat.any(1));
          if(!idxKeep.any()) {
            err = true;
            return;
          }
          posSol(ori,line).keepRows(idxKeep);
        }
        Matrix<BOOL> v1 = posSol(ori,line).all(0);
        sol.setTrue(line,v1);
        solSet.setTrue(line,v1);
        Matrix<BOOL> v2 = !(posSol(ori,line).any(0));
        sol.setFalse(line,v2);
        solSet.setTrue(line,v2);
      }
      sol.transpose();
      solSet.transpose();
      solSetChange.transpose();
    }
    sol.printBoolToConsole();
    if(solSet.all())
      return;
    if(sol_old == sol) {
      Matrix<Matrix<BOOL> > thPosSol = posSol;
      incLevel++;
      size_t min_s = 1e6;
      size_t minOri = 2;
      size_t minLine = 10000;
      for(size_t ori=0;ori<2;ori++)  {
        for(size_t line=0;line<dim[ori];line++) {
          if (1 < thPosSol(ori,line).dimX && thPosSol(ori,line).dimX < min_s) {
            min_s = thPosSol(ori,line).dimX;
            minOri = ori;
            minLine = line;
          }
        }
      }
      thPosSol(minOri,minLine).keepFirstLine();
      cout << ind << "Guess correct solution in O" << minOri << "L" << minLine << " of " << min_s << " solutions there." << endl;
      Matrix<BOOL> thSol = sol;
      Matrix<BOOL> thSolSet = solSet;
      fun_logics_rec(thSol,thSolSet,incLevel,thPosSol,err);
      if(err) {
        cout << ind << "Invalid solution in O" << minOri << "L" << minLine << " found and dropped, remaining " << min_s-1 << " solutions." << endl;
        posSol(minOri,minLine).dropFirstLine();
        incLevel--;
      } else {
        sol = thSol;
        cout << ind << "Returned from inception level " << incLevel << ". Solved!" << endl;
      }
    }
  }
}

int main() {
  setlocale(LC_ALL, "");
  size_t inpNr = 41;
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
  cout << "Creation started!" << endl;
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
        std::stringstream ss2;
        ss2 << "comp_" << (size_t) noWhiteBlocks << "_" << (size_t) noWhiteSquares << ".dat";
        if (!exist(ss2.str())) fun_cr_comps(noWhiteBlocks,noWhiteSquares);
        size_t noComps = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1)+2*alg_n_k(noWhiteSquares-1,noWhiteBlocks-2)+alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
        Matrix<BOOL> S(noComps,dim[ori]);
        S.setFalse();
        MatrixB S2(noComps,dim[ori]);
        cout << "Creating sol O" << ori << "L" << line << ", loading combination (" << (int) noWhiteBlocks << " " << (int) noWhiteSquares << ") with " << noComps << " compositions." << endl;
        Matrix<char> A(noComps,noWhiteBlocks);
        A.load(ss2.str());
        for(size_t k=0;k<noComps;k++) {
          size_t pos = 0;
          for(size_t l=0;l<black.size();l++) {
            pos += A(k,l);
            for (size_t m = 0 ; m < black[l] ; m++) {
              S(k,pos) = true;
              S2(k,pos) = true;
              pos++;
            }
          }
        }
        cout << "Finished sol O" << ori << "L" << line << ", using combination (" << (int)  noWhiteBlocks << " " << (int)  noWhiteSquares << ") with " << noComps << " compositions." << endl;

        std::stringstream ss3;
        ss3 << inpNr << "/ori" << ori << "_line" << line << ".dat11";
        S.write(ss.str());
        S2.write(ss3.str());
      }
    }
  }

  Matrix<BOOL> sol = Matrix<char>(dim[1],dim[0]);
  sol.setFalse();
  Matrix<BOOL> solSet = Matrix<char>(dim[1],dim[0]);
  solSet.setFalse();
  Matrix<Matrix<BOOL> > posSol = Matrix<Matrix<BOOL> >(2,maxDim);
  for(size_t ori=0;ori<2;ori++) {
    for(size_t line=0;line < dim[1-ori];line++) {
      std::vector<uint8_t> black = M[ori]->operator[](line);
      uint8_t noWhiteBlocks = black.size() + 1;
      uint8_t noWhiteSquares = dim[ori] - accumulate(black.begin(),black.end(),0);
      std::stringstream ss;
      ss << inpNr << "/ori" << ori << "_line" << line << ".dat";
      size_t noComps = alg_n_k(noWhiteSquares-1,noWhiteBlocks-1)+2*alg_n_k(noWhiteSquares-1,noWhiteBlocks-2)+alg_n_k(noWhiteSquares-1,noWhiteBlocks-3);
      Matrix<BOOL> S = Matrix<BOOL>(noComps,dim[ori]);
      S.load(ss.str());
      posSol(ori,line) = S;
    }
  }
  cout << "Creation finished!" << endl;

  //Iteration
  cout << "Iteration started!" << endl;
  auto t_start = std::chrono::high_resolution_clock::now();
  bool err = false;
  size_t incLevel = 0;
  fun_logics_rec(sol,solSet,incLevel,posSol,err);
  auto t_end = std::chrono::high_resolution_clock::now();
  double calcTime = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
  cout << "Iteration finished!" << endl;
  cout << "Took " << calcTime << " milliseconds." << endl;

  sol.printBoolToConsole();

  //Plot
  //solSet = true | sol
  //scr_plot

  return 0;
}
