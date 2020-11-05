// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cwctype>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <regex>
#include "hypersurface.h"






typedef struct {
  string label;
  int64_t p;
  map< Vec<int64_t>, ZZ, vi64less> f;
  Mat<int64_t> A;
  Vec<int64_t> b;
  Mat<int64_t> monomials;
  Vec<ZZ> coefficients;

  Vec<ZZ> zeta;
  Mat<ZZ> frob;
  Vec<int64_t> hodge_numbers;
  Vec<int64_t> r_vector;
} hypersurface;





typedef std::chrono::time_point<std::chrono::system_clock> SystemTime ;


void do_hypersurface(hypersurface& X, int64_t verbose=0, int64_t abs_precision=0) {
  zeta_and_frob_Fp(X.zeta, X.frob, X.hodge_numbers, X.r_vector, X.p, X.f, X.A, X.b, verbose, abs_precision);
}


// expects as input label:monomials:coefficients:A:b:p
istream& operator>>(istream &is,  hypersurface &o) {
  for(size_t i = 0; i < 6; ++i) {
    string s;
    if(not getline(is, s, ':'))
      throw_line("missing field"s);
    stringstream ss(s);

    switch(i) {
      case 0:
        {
          if(!(ss >> o.label))
            throw_line("bad label"s);
          break;
        }
      case 1:
        if(!(ss >>= o.monomials))
          throw_line("bad monomials"s);
        break;
      case 2:
        {
          if(!(ss >>= o.coefficients))
            throw_line("bad coefficients"s);
          assert_print(o.monomials.NumRows(), ==, o.coefficients.length());
          for(long i = 0; i < o.coefficients.length(); ++i)
            o.f[o.monomials[i]] = o.coefficients[i];
          break;
        }
      case 3:
        if(!(ss >>= o.A))
          throw_line("bad A"s);
        break;
      case 4:
        if(!(ss >>= o.b))
          throw_line("bad b"s);
        break;
      case 5:
        if(!(ss >> o.p))
          throw_line("bad prime"s);
        break;
      default:
        throw_line("too many fields in the line!"s);
    }
  }
  return is;
}

// outputs label:f:Ap:bP:p:hodge:zeta
ostream& operator<<(ostream &s, const hypersurface &X) {

  stringstream buffer;
  buffer << X.label << ":";
  buffer <<= X.monomials;
  buffer << ":";
  buffer <<= X.coefficients;
  buffer << ":";
  buffer <<= X.A;
  buffer << ":";
  buffer <<= X.b;
  buffer << ":";
  buffer << X.p<< ":";
  buffer <<= X.hodge_numbers;
  buffer << ":";
  buffer <<= X.zeta;
  string buffer_str = buffer.str();
  buffer_str.erase(std::remove(buffer_str.begin(), buffer_str.end(), '\n'), buffer_str.end());
  buffer_str.erase(std::remove(buffer_str.begin(), buffer_str.end(), ' '), buffer_str.end());
  //std::replace(buffer_str.begin(), buffer_str.end(), '\n', '');
  std::replace(buffer_str.begin(), buffer_str.end(), '(', '[');
  std::replace(buffer_str.begin(), buffer_str.end(), ')', ']');
  s << buffer_str;
  return s;
}


int message(int argc, char **argv) {
  fprintf(stderr, "Usage:\n %s input_filename output_filename\n", argv[0]);
  fprintf(stderr, "Where each row follows the following format:\n");
  fprintf(stderr, "\tlabel:monomials:coefficients:A:b:p\n");
  fprintf(stderr, "where,\n");
  fprintf(stderr, "\t- label, some string of text that identifies the hypersurface\n");
  fprintf(stderr, "\t- monomials, the monomials of the laurent polynomial f defining hypersurface X\n");
  fprintf(stderr, "\t- coefficients, the corresponding coefficients\n");
  fprintf(stderr, "\t- A, the matrix associated to the half space representation of the polyhedron P = convex(supp(f)), st w in k*P <=>  A * v + k *b >= 0\n");
  fprintf(stderr, "\t- b, the vector associated to the half space representation of the polyhedron P, st w in k*P <=>  AP * v + k *bP >= 0\n");
  fprintf(stderr, "\t- p, the prime p\n");
  if(argc >= 2) {
    fprintf(stderr, "Arguments given:\n");
    fprintf(stderr, "\t<intput_filename> = %s\n", argv[1]);
  }
  if(argc >= 3) {
    fprintf(stderr, "\t<output_filename> = %s\n", argv[2]);
  }
  return false;
}

int main (int argc, char**argv) {
    if(argc < 3)
      return message(argc, argv);

    ifstream input(argv[1]);
    ofstream output(argv[2]);
    string   line;
    int r = 0;
    while(std::getline(input, line)) {
      try {
      SystemTime start(std::chrono::system_clock::now());
      std::time_t startt = std::chrono::system_clock::to_time_t(start);
      cout << "Date:   \t" <<  std::put_time(std::localtime(&startt), "%F %T") << endl;

      hypersurface X;
      // read a line
      stringstream linestream(line);
      linestream >> X;
      cout << "Starting:\t"<< X.label<<endl;
      do_hypersurface(X, 0);
      output << X << endl;
      //cout <<= X.frob;
      //cout << endl;

      SystemTime end(std::chrono::system_clock::now());
      std::time_t endt = std::chrono::system_clock::to_time_t(end);
      double walltime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
      cout << "Date:   \t" << std::put_time(std::localtime(&endt), "%F %T") << endl;
      cout << "Done:   \t"<< X.label << "\ttook ";
      cout << std::setw(6) << std::setfill(' ')  << std::fixed << std::setprecision(2) << walltime/1000 << "s"<< endl << endl;
      } catch( const std::exception & ex ) {
        cerr << "Uncaught exception: " <<ex.what() << endl;
        cerr << "Moving to next input line" << endl;
        //std::abort();
      }
    }
    return r;
}

