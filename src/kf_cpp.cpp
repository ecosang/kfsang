//http://www.mjdenny.com/Rcpp_Intro.html
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]

List kf_cpp(const arma::mat& A,const arma::mat& B,const arma::mat& C,const arma::mat& D,
               const arma::mat& Q, const arma::mat& R, const arma::mat& mu0, const arma::mat& P0,
               const arma::mat& u_mat, const arma::mat& y_mat) {
  int nx =A.n_rows;
  int ny=y_mat.n_rows;
  int TT=y_mat.n_cols;
  //int nu =u_mat.n_rows;
  
  //arma::vec new_vec = arma::zeros(len);
  //arma::mat new_mat = arma::zeros(len,len);
  //arma::cube new_array = arma::zeros(len,len,len);
  arma::cube xp = arma::zeros(nx,1,TT);
  arma::cube Pp = arma::zeros(nx,nx,TT);
  arma::cube xf = arma::zeros(nx,1,TT);
  arma::cube Pf = arma::zeros(nx,nx,TT);
  arma::cube innov = arma::zeros(ny,1,TT);
  arma::cube sig = arma::zeros(ny,ny,TT);
  arma::mat sigtemp = arma::zeros(ny,ny);
  arma::mat siginv = arma::zeros(ny,ny);
  arma::mat K = arma::zeros(nx,nx);
  arma::mat sigmat = arma::zeros(ny,ny);
  // https://thecoatlessprofessor.com/programming/common-operations-with-rcpparmadillo/
  xp.slice(0)=mu0;//xp[,,0]
  Pp.slice(0)=P0;//Pp[,,0]
  arma::vec ll = arma::zeros(1);
  
  for (int t=0; t < TT; ++t) {
    
    sigtemp=C*Pp.slice(t)*trans(C)+R; // https://stackoverflow.com/questions/28465766/matrix-multiplication-using-numericmatrix-and-numericvector-in-rcpp
    sig.slice(t)=(trans(sigtemp)+sigtemp)/2;
    siginv=inv(sig.slice(t));
    K=Pp.slice(t)*trans(C)*siginv;
    innov.slice(t)=y_mat.col(t)-C*xp.slice(t)-D*u_mat.col(t);
    xf.slice(t)=xp.slice(t)+K*innov.slice(t);
    Pf.slice(t)=Pp.slice(t)-K*C*Pp.slice(t);
    sigmat=sig.slice(t);
    //http://dirk.eddelbuettel.com/papers/RcppArmadillo.pdf
    ll= ll + log(det(sigmat)) + trans(innov.slice(t))*siginv*innov.slice(t)+ny*log(2*3.14159265); //-log(likelihood) ,
    
      if (t == (TT-1)){
        
      }else{
        xp.slice(t+1)=A*xf.slice(t) + B*u_mat.col(t);
        Pp.slice(t+1)=A*Pf.slice(t)*trans(A)+Q;
      } 
  }
  ll=0.5*ll;
  List out(8);
  return Rcpp::List::create(Rcpp::Named("xp")=xp,
                            Rcpp::Named("Pp")=Pp,
                            Rcpp::Named("xf")=xf,
                            Rcpp::Named("Pf")=Pf,
                            Rcpp::Named("ll")=ll,
                            Rcpp::Named("innov")=innov,
                            Rcpp::Named("sig")=sig,
                            Rcpp::Named("ll")=ll
  );
}


/*** R
# A=nx x nx , B=nx x nu, Q=nx x nx, u_mat=nu x TT  
#  y_mat: ny x TT, C= ny x nx, D= ny x nu, R=ny x ny
kf_cpp(A=matrix(2,nrow=2,ncol=2),B=matrix(3,nrow=2,ncol=1),C=matrix(0.5,nrow=3,ncol=2),D=matrix(1,nrow=3,ncol=1),
          Q=diag(1,2),R=diag(1.5,3),mu0=matrix(3,nrow=2,ncol=1),P0=diag(4,2),
          y_mat=matrix(1,nrow=3,ncol=5),u_mat=matrix(1,nrow=1,ncol=5))
  */