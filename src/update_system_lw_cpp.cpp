#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]

List update_system_lw_cpp(const arma::vec& x0,const arma::vec& par, const arma::vec& dims, int& dt) {
  //input y_mat, u_mat, dt, A,B,C,D, x0
  
  
  //int nt =u_mat.n_cols;
  int nx=dims(0);
  int nu=dims(1);
  int ny=dims(2);
  arma::mat AA=arma::zeros(nx,nx);
  arma::mat BB=arma::zeros(nx,nu);
  // A  B
  // A0 B0
  arma::mat AB0=arma::zeros(nu,nx+nu);
  arma::mat CC=arma::zeros(ny,nx);
  arma::mat DD=arma::zeros(ny,nu);
  //arma::vec x0=arma::zeros(nx);
  
  double Te0=x0(0);
  double Ti0=x0(1);
  double Tm0=x0(2);
  double Ts0=x0(3);
  double Tn0=x0(4);
  
  double Ce=par(0);
  double Ci=par(1);
  double Cm=par(2);
  double Cs=par(3);
  
  double Rea=par(4);
  double Rie=par(5);
  double Rim=par(6);
  double Ris=par(7);
  double Rnm=par(8);
  
  double Aw=par(9);
  double kig=par(10);
  double khphtg=par(11);
  double kaux=par(12);
  double kdf=par(13);
  double khpclg=par(14);
  
  double sdte=par(15);
  double sdti=par(16);
  double sdtm=par(17);
  double sdts=par(18);
  double sdtn=par(19);
  
  AA(0,0)=-1/(Rea*Ce)-1/(Rie*Ce);
  AA(0,1)=1/(Rie*Ce);
  AA(1,0)=1/(Rie*Ci);
  //AA(1,1)=-(1/(Rie*Ci)+1/(Rim*Ci)+1/(Ris*Ci)+1/(Rin*Ci));
  AA(1,1)=-(1/(Rie*Ci)+1/(Rim*Ci)+1/(Ris*Ci));
  AA(1,2)=1/(Rim*Ci);
  AA(1,3)=1/(Ris*Ci);
  //AA(1,4)=1/(Rin*Ci);

  AA(2,1)=1/(Rim*Cm);
  AA(2,2)=-1/(Rim*Cm)-1/(Rnm*Cm);
  AA(2,4)=1/(Rnm*Cm);
  AA(3,1)=1/(Ris*Cs);
  AA(3,3)=-1/(Ris*Cs);
  
  
  BB(0,0)=1/(Rea*Ce);
  BB(1,1)=Aw/Ci;
  BB(1,2)=kig/Ci;
  BB(1,3)=khphtg/Ci;
  BB(1,4)=kaux/Ci;
  BB(1,5)=kdf/Ci;
  BB(1,6)=khpclg/Ci;
  CC(0,3)=1;
  
  arma::mat mat1= arma::join_vert(arma::join_horiz(AA,BB)*dt,AB0);
  arma::mat smat= arma::expmat(mat1);
  arma::mat Ad= smat.submat(0,0,nx-1,nx-1);//row,col, row, col
  arma::mat Bd= smat.submat(0,nx,nx-1,nx+nu-1);//row,col, row, col
  
  List out(16);
  return Rcpp::List::create(
    Rcpp::Named("Ad")=Ad,
    Rcpp::Named("Bd")=Bd,
    Rcpp::Named("Cd")=CC,
    Rcpp::Named("Dd")=DD,
    Rcpp::Named("dt")=dt,
    Rcpp::Named("dims")=dims,
    Rcpp::Named("x0")=x0,
    Rcpp::Named("sdte")=sdte,
    Rcpp::Named("sdti")=sdti,
    Rcpp::Named("sdtm")=sdtm,
    Rcpp::Named("sdts")=sdts,
    Rcpp::Named("sdtn")=sdtn,
    Rcpp::Named("khphtg")=khphtg,
    Rcpp::Named("kaux")=kaux,
    Rcpp::Named("kdf")=kdf,
    Rcpp::Named("khpclg")=khpclg
    
  );
}

