//http://www.mjdenny.com/Rcpp_Intro.html
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]

List nstep_cpp(const arma::mat& A,const arma::mat& B,const arma::mat& C,const arma::mat& D,
            const arma::mat& x0,const arma::mat& u_mat) {
  //input y_mat, u_mat, dt, A,B,C,D, x0
  int nt =u_mat.n_cols;
  int nu =u_mat.n_rows;
  int nx=A.n_rows;
  int ny=C.n_rows;
  double na_check=0;
  arma::mat x_mat=arma::zeros(nx,nt);
  arma::mat y_pred=arma::zeros(ny,nt);
  x_mat.col(0)=x0;
  for (int t=1; t<nt;++t){
    if(u_mat.col(t-1).has_nan()){
      x_mat.col(t)=x_mat.col(t-1);
      y_pred.col(t-1)=C*x_mat.col(t-1);
    }else{
      x_mat.col(t)=A*x_mat.col(t-1)+B*u_mat.col(t-1);//x_next;
      y_pred.col(t-1)=C*x_mat.col(t-1)+D*u_mat.col(t-1);
    }
    
  }
  if(u_mat.col(nt-1).has_nan()){
    y_pred.col(nt-1)=C*x_mat.col(nt-1)+D*u_mat.col(nt-1);
  }else{
    y_pred.col(nt-1)=C*x_mat.col(nt-1)+D*u_mat.col(nt-1);
  }
  //int ny=y_mat.n_rows;
  //int TT=y_mat.n_cols;
  //int nu =u_mat.n_rows;
  List out(2);
  return Rcpp::List::create(
	Rcpp::Named("y_pred")=y_pred,
	Rcpp::Named("x_pred")=x_mat
  );
}



/*** R
# A=nx x nx , B=nx x nu, Q=nx x nx, u_mat=nu x TT  
#  y_mat: ny x TT, C= ny x nx, D= ny x nu, R=ny x ny
 nstep_rcpp(A=matrix(2,nrow=2,ncol=2),B=matrix(3,nrow=2,ncol=1),C=matrix(0.5,nrow=3,ncol=2),D=matrix(1,nrow=3,ncol=1),
      x0=matrix(3,nrow=2,ncol=1),u_mat=matrix(1,nrow=1,ncol=5))
  */