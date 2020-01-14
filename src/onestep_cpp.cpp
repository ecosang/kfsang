//http://www.mjdenny.com/Rcpp_Intro.html
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]


List onestep_cpp(const arma::mat& A,const arma::mat& B,const arma::mat& C,const arma::mat& D,
                const arma::mat& x,const arma::mat& u) {
  //input y_mat, u_mat, dt, A,B,C,D, x0
  int nx=A.n_rows;
  int ny=C.n_rows;
  arma::mat x_pred=arma::zeros(nx,1);
  arma::mat y_curr=arma::zeros(ny,1);
  x_pred=A*x+B*u; //x_pred;
  y_curr=C*x+D*u; // y_curr;
  
  List out(2);
  return Rcpp::List::create(Rcpp::Named("y")=y_curr,
                            Rcpp::Named("x")=x_pred
  );

}

