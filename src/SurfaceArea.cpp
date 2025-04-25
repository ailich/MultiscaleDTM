#include <Rcpp.h>
using namespace Rcpp;


// Calculate area of triangle based on side lengths
// [[Rcpp::export]]
double C_TriArea (double a, double b, double c){
  double s = (a+b+c)/2;
  double out =sqrt(s*(s-a)*(s-b)*(s-c));
  return out;
}

// [[Rcpp::export]]
NumericVector C_SurfaceArea(const NumericVector& z, double x_res, double y_res, bool na_rm, size_t ni, size_t nw) {
  NumericVector out(ni, NA_REAL);
  
  double Lx2 = x_res * x_res;
  double Ly2 = y_res * y_res;
  double Ld2 = Lx2 + Ly2;
  
  const double* z_ptr = z.begin();
  
  for (size_t i = 0; i < ni; ++i) {
    const double* zw = z_ptr + i * nw;
    
    // Fast access to elevation values
    double A = zw[0], B = zw[1], C = zw[2], D = zw[3], E = zw[4], 
                                                             F = zw[5], G = zw[6], H = zw[7], I = zw[8];
    //|A|B|C|
    //|D|E|F|
    //|G|H|I|
    
    // Horizontal
    double AB = sqrt(Lx2 + pow(A - B, 2)) / 2;
    double BC = sqrt(Lx2 + pow(B - C, 2)) / 2;
    double DE = sqrt(Lx2 + pow(D - E, 2)) / 2;
    double EF = sqrt(Lx2 + pow(E - F, 2)) / 2;
    double GH = sqrt(Lx2 + pow(G - H, 2)) / 2;
    double HI = sqrt(Lx2 + pow(H - I, 2)) / 2;
    // Vertical
    double AD = sqrt(Ly2 + pow(A - D, 2)) / 2;
    double BE = sqrt(Ly2 + pow(B - E, 2)) / 2;
    double CF = sqrt(Ly2 + pow(C - F, 2)) / 2;
    double DG = sqrt(Ly2 + pow(D - G, 2)) / 2;
    double EH = sqrt(Ly2 + pow(E - H, 2)) / 2;
    double FI = sqrt(Ly2 + pow(F - I, 2)) / 2;
    
    // Diagonal
    double EA = sqrt(Ld2 + pow(E - A, 2)) / 2;
    double EC = sqrt(Ld2 + pow(E - C, 2)) / 2;
    double EG = sqrt(Ld2 + pow(E - G, 2)) / 2;
    double EI = sqrt(Ld2 + pow(E - I, 2)) / 2;
    
    // Triangles
    double tri[8] = {
      C_TriArea(EA, AB, BE),
      C_TriArea(BE, BC, EC),
      C_TriArea(AD, DE, EA),
      C_TriArea(EC, CF, EF),
      C_TriArea(DE, DG, EG),
      C_TriArea(EF, FI, EI),
      C_TriArea(EG, EH, GH),
      C_TriArea(EH, EI, HI)
    };
    
    if (na_rm) {
      double sum = 0;
      int count = 0;
      for (int t = 0; t < 8; ++t) {
        if (!std::isnan(tri[t])) {
          sum += tri[t];
          count++;
        }
      }
      if (count > 0) out[i] = (sum / count) * 8;
    } else {
      bool has_na = false;
      for (int t = 0; t < 8; ++t) {
        if (std::isnan(tri[t])) {
          has_na = true;
          break;
        }
      }
      if (!has_na) {
        out[i] = tri[0] + tri[1] + tri[2] + tri[3] + tri[4] + tri[5] + tri[6] + tri[7];
      }
    }
  }
  
  return out;
}