#include<RcppEigen.h>
#include "EigenQP.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;
using Eigen::MatrixXf;
using Eigen::VectorXf;

//C++ function to call the R function quadprog::solve.QP (for comparison)
//
//Call the R function quadprog::solve.QP to solve a quadratic programming problem
//
//@param Dmat,dvec,Amat See \code{\link[quadprog]{solve.QP}} for explanation.
//
//@return a list containing the solution
List quadprog_solveR(const MatrixXd & Dmat, const VectorXd & dvec, const MatrixXd & Amat)
{
  List result;
  Environment myEnv = Environment::namespace_env("quadprog");
  Function quadR = myEnv["solve.QP"];
  result = quadR(Dmat, dvec, Amat);
  return result;
}

//C++ function solving quadratic programming implemented by QuadProg++
//
//Utilized the codes by wingsit. See EigenQP.h and EigenQP.cpp for details. This function restricts the first component of the solution
//  to be 0, to avoid unidentifiability.
//
//@param Dmat,dvec,At See \code{\link[quadprog]{solve.QP}} for explanation. This function is wrapped up in such a way
//	that it works the same as quadprog::solve.QP. At means Amat transpose.
//
//@return a list containing the solution and the value
//'@export
//[[Rcpp::export]]
List quadprog_solveC2(const MatrixXd & Dmat, const VectorXd & dvec, const MatrixXd & At)
{
  //Because the solve_quadprog() function requires at least one equality constraint, I force the first
  // of the latent beta's to be 0.
  int n = dvec.size();
  int m = At.cols();

  MatrixXd G(Dmat); //because the solve_quadprog function changes the input G, so make a copy
  VectorXd g0 = (-1.)*dvec; //different definition of QuadProg++ and quadprog R package

  VectorXd x(n); //store the result of beta's

  MatrixXd CE = MatrixXd::Zero(n, 1);
  //CE(0,0) = 1.; //equality constrant to force the first beta to be 0.
  Rcout<<"CE = "<<CE<<std::endl;

  VectorXd ce0(1);
  ce0(0) = 0;

  VectorXd ci0 = VectorXd::Zero(m);

  double value = QP::solve_quadprog(G, g0, CE, ce0, At, ci0, x);

  List res = List::create(Named("solution") = x, Named("value") = value);
  return res;
}

//C++ function solving quadratic programming implemented by QuadProg++
//
//utilized the codes by wingsit. See EigenQP.h and EigenQP.cpp for details.
//
//@param Dmat,dvec See \code{\link[quadprog]{solve.QP}} for explanation. This function is wrapped up in such a way
//	that it works the same as quadprog::solve.QP. At means Amat transpose.
//@param At transpose of the partial order matrix: t(At)*beta>=0 gives the inequality constrain
//@param Bt tranpose of the equality constrain matrix: t(Bt)*beta=0 gives the equality constrain
//
//@return a list containing the solution and the value
//'@export
//[[Rcpp::export]]
List quadprog_solveC1(const MatrixXd & Dmat, const VectorXd & dvec, const MatrixXd & At, const MatrixXd & Bt)
{
  //Because the solve_quadprog() function requires at least one equality constraint, I mannually add one more
  //dimension to the variable.
  int n = dvec.size();
  int mi = At.cols(); //number of inequality constrains
  int me = Bt.cols(); //number of equality constrains

  MatrixXd G(n+1, n+1);
  G << Dmat, MatrixXd::Zero(n, 1),
       MatrixXd::Zero(1, n), 1;


  VectorXd g0(n+1);
  g0 << (-1.)*dvec, 0; //different definition of QuadProg++ and quadprog R package

  VectorXd x(n+1); //store the result of beta's

  MatrixXd CE(n+1, me);
  CE << Bt, MatrixXd::Zero(n, 1),
        MatrixXd::Ones(1, me), 1;

  VectorXd ce0 = VectorXd::Zero(me+1);

  MatrixXd CI(n+1, mi);
  CI << At,
        MatrixXd::Zero(1, mi);

  VectorXd ci0 = VectorXd::Zero(mi);

  double value = QP::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

  List res = List::create(Named("solution") = x.head(n), Named("value") = value);
  return res;
}

//C++ function solving quadratic programming implemented by QuadProg++
//
//utilized the codes by wingsit. See EigenQP.h and EigenQP.cpp for details.
//
//@param Dmat,dvec See \code{\link[quadprog]{solve.QP}} for explanation. This function is wrapped up in such a way
//	that it works the same as quadprog::solve.QP. At means Amat transpose.
//@param At transpose of the partial order matrix: t(At)*beta>=0 gives the inequality constrain
//@param Bt tranpose of the equality constrain matrix: t(Bt)*beta=0 gives the equality constrain
//
//@return a list containing the solution and the value
//'@export
//[[Rcpp::export]]
List quadprog_solveC(const MatrixXd & Dmat, const VectorXd & dvec, const MatrixXd & At)
{
  //Because the solve_quadprog() function requires at least one equality constraint, I mannually add one more
  //dimension to the variable.
  int n = dvec.size();
  int m = At.cols(); //number of inequality constrains

  MatrixXd G(n+1, n+1);
  G << Dmat, MatrixXd::Zero(n, 1),
       MatrixXd::Zero(1, n), 1;


  VectorXd g0(n+1);
  g0 << (-1.)*dvec, 0; //different definition of QuadProg++ and quadprog R package

  VectorXd x(n+1); //store the result of beta's

  MatrixXd CE(n+1, 1);
  CE << MatrixXd::Zero(n, 1),
        1;

  VectorXd ce0 = VectorXd::Zero(1);

  MatrixXd CI(n+1, m);
  CI << At,
        MatrixXd::Zero(1, m);

  VectorXd ci0 = VectorXd::Zero(m);

  double value = QP::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

  List res = List::create(Named("solution") = x.head(n), Named("value") = value);
  return res;
}

//This function perform the special summation over risk sets.
//
//@param indset the index matrix of risk set
//@param value_vec the values to be summed over. it takes value exp(eta) in our problem
//@param ind_vec the vector indicating the failure times to be summed over
double sum_risk(const MatrixXd & indset, const VectorXd & value_vec, const VectorXd & ind_vec, double p){
  double sum = 0;
  for (int i = 0; i < ind_vec.size(); i++){
    if (ind_vec(i) == 1.) sum += pow(indset.col(i).dot(value_vec), -p);
  }
  return sum;
}

//Function to calculate the negative log partial likelihood
//
//@param X design matrix generated by \code{\link{data.prep}}
//@param indset riskset matrix generated by \code{\link{data.prep}}
//@param beta vector of parameter beta's
//@param cen vector of censoring indicator
//
//@return the negative log partial likelihood
//'@export
//[[Rcpp::export]]
double logPL_C(const MatrixXd & X, const MatrixXd & indset, const VectorXd & cen, const VectorXd & beta)
{
  int n_fail = indset.cols();

  VectorXd eta = X*beta;
  VectorXd expeta = exp(eta.array());

  double sum = 0;
  for (int i = 0; i < n_fail; i++){
    VectorXd col = indset.col(i);
    sum += log(col.dot(expeta));
  }

  sum = sum - cen.dot(eta);
  return(sum);
}

//Function to perform the iterative solution of the group fused lasso penalty
//
//@param Dmat0 the original quadratic term matrix
//@param dvec0 the original linear term vector
//@param Apo the partial order matrix
//@param lambda the tuning parameter
//@param w adaptive param
//@param nb number of boundaries
//@param m number of data sets
//@param inibeta initial beta, for warm start
//@param maxiter maximum number of iteration
//@param eps tolerance
//@param trace debug
//
//@return solution beta
VectorXd iter_pen(const MatrixXd & Dmat0, const VectorXd & dvec0, const MatrixXd & Apo, const double lambda, const VectorXd & w, const int nb, const int m,
                  const VectorXd & inibeta, const int maxiter, const double eps, const bool trace = false)
{
  double delta = 1e-8; //small cutoff for group norm in the denominator
  int pqm = Apo.cols();
  if(Apo.rows() != nb*m || Dmat0.cols() != pqm) Rcout<<"Dimension inconsistent!"<<std::endl;
  VectorXd currbeta = inibeta;
  double diffbeta = 1.;
  int itt = 0;
  VectorXd oldbeta(pqm);
  while(itt < maxiter && diffbeta > eps){
    oldbeta = currbeta;

    MatrixXd Q = MatrixXd::Zero(pqm, pqm); //quadratice contribution of the numerator of the penalty
    for(int j = 0; j < nb; j++){
      VectorXd g = VectorXd::Zero(nb);
      g(j) = 1.;
      MatrixXd G = g.replicate(m,1).asDiagonal();
      MatrixXd PGP = Apo.transpose()*G*Apo;
      double norm_g = sqrt(oldbeta.transpose()*PGP*oldbeta);
      //protection from 0 denominator: if norm_g<delta, that PGP would be amplified and force corresponding boundary to 0. Besides, this
      //precision should be higher than the convergence precision.
      if (norm_g > delta)
        Q += (2/norm_g)*PGP*w(j);
      else Q += (2/(delta))*PGP*w(j);
    }

    MatrixXd Dmat = Dmat0 + lambda*Q + 0.01*eps*MatrixXd::Identity(pqm, pqm);
    List QP = quadprog_solveR(delta*Dmat, delta*dvec0, Apo.transpose()); //Rescale both Dmat and dvec to avoid overflow of quadprog::sovle.QP
    currbeta = QP["solution"];
    diffbeta = (currbeta - oldbeta).norm();

    if(trace){
      Rcpp::Rcout<<"penalty iteration "<<itt<<":"<<std::endl;
      double val = QP["value"];
      Rcpp::Rcout<<"QP value: "<<val<<std::endl;
      Rcpp::Rcout<<"diff beta: "<<diffbeta<<std::endl;
      Rcpp::Rcout<<"current beta: "<<currbeta.transpose()<<std::endl;
      Rcpp::Rcout<<std::endl;
    }

    itt += 1;
  }

  // if(itt == maxiter){
  //   Rcpp::Rcout<<"Iterative updating of penalty doesn't converge! Try increasing maxiter. "<<std::endl;
  // }

  return currbeta;
}

//function to calculate the value of the penalty term
//
//@param Apo partial order matri
//@param w adaptive param
//@param nb number of boundaries of the TN grid
//@param m number of dataset
//@param beta vector of beta
//
//@return value of the penalty term
//'@export
//[[Rcpp::export]]
double val_pen(const MatrixXd Apo, const VectorXd & w,  const int nb, const int m, const VectorXd beta){
  if (Apo.rows() != nb*m) Rcout<<"Inconsistant dimensions of the partial order matrix!"<<std::endl;

  double val = 0;
  for(int j = 0; j<nb; j++){
    VectorXd g = VectorXd::Zero(nb);
    g(j) = 1;
    MatrixXd G = g.replicate(m, 1).asDiagonal();
    double norm = sqrt(beta.transpose()*Apo.transpose()*G*Apo*beta);
    val += w(j)*norm;
  }
  return(val);
}

//function to round close values to the same
//
//@param beta vector of beta's
//@param eps tolerance of rounding
//'@export
//[[Rcpp::export]]
VectorXd roundvec(const VectorXd & beta, const double eps){
  int n = beta.size();
  VectorXd beta_round(beta);
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      if(std::abs(beta(j) - beta(i))< eps) beta_round(j) = beta_round(i);
    }
  }
  return beta_round;
}


//function to calculate the best beta's for a single point(lambda1, lambda2)
//
//@param X,indset,Apo matrices generated by the lasso.tree function in func_sol.R, representing covariates, risk set indicator and partial order
//@param nb,m number of boundaries and number of data sets
//@param cen censoring indicator vector
//@param lambda tuning parameter
//@param w adaptive param
//@param maxiter,eps,trace parameters from the lasso.tree function
//@param fullA bool variable whether calculate all elements of the Hessian, internal test only, not for users
//
//@return a vector of the solution beta
//[[Rcpp::export]]
VectorXd lasso_tree_single(const MatrixXd & X, const MatrixXd & indset, const MatrixXd & Apo, const int nb, const int m, const VectorXd & cen,
                           const double lambda, const VectorXd & w, const VectorXd & inibeta, const int maxiter, const double eps, const bool trace, const bool fullA = false)
{
  int n = X.rows();
  int p = X.cols();
  VectorXd currbeta = inibeta;
  int itt = 0;
  double diffbeta = 1.;
  VectorXd oldbeta(p);
  MatrixXd Dmat(p, p); //quadratic term coefficient in irls
  VectorXd dvec(p); //linear term coefficient in irls
  double oldval, currval, diffval; //tracking the objective values
  if(trace){
    currval = logPL_C(X, indset, cen, currbeta) + lambda*val_pen(Apo, w, nb, m, currbeta);
    diffval = 1.;
  }

  while(itt < maxiter && diffbeta > eps){
    oldbeta = currbeta;
    if(trace) oldval = currval;

    VectorXd eta = X*oldbeta;
    VectorXd expeta = exp(eta.array());

    VectorXd u(n);
    VectorXd d(n);
    for (int i = 0; i < n; i++){
      u(i) = cen(i) - expeta(i)*sum_risk(indset, expeta, indset.row(i), 1);
      d(i) = u(i) - cen(i) + exp(2*eta(i))*sum_risk(indset, expeta, indset.row(i), 2);
    }
    MatrixXd A(-1.*d.asDiagonal());

    if(fullA){
      for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++){
          VectorXd indvec((indset.row(i).array())*(indset.row(j).array()));
          A(j, i) = -2*expeta(i)*expeta(j)*sum_risk(indset, expeta, indvec, 2);
        }
        A = (A + A.transpose())/2.;
    }

    Dmat = X.transpose()*A*X;
    dvec = X.transpose()*(u + A*eta);

    currbeta = iter_pen(Dmat, dvec, Apo, lambda, w, nb, m, oldbeta, 1, eps, false);
    diffbeta = (currbeta - oldbeta).norm();

    if(trace){
      currval = logPL_C(X, indset, cen, currbeta) + lambda*val_pen(Apo, w, nb, m, currbeta);
      diffval = currval - oldval;
      Rcpp::Rcout<<"iteration "<<itt<<":"<<std::endl;
      //Rcpp::Rcout<<"old val:"<<oldval<<std::endl;
      //Rcpp::Rcout<<"curr val:"<<currval<<std::endl;
      Rcpp::Rcout<<"diff val:"<<diffval<<std::endl;
      Rcpp::Rcout<<"current objective value: "<<currval<<std::endl;
      Rcpp::Rcout<<"diff beta: "<<diffbeta<<std::endl;
      Rcpp::Rcout<<"current beta: "<<currbeta.transpose()<<std::endl;
      Rcpp::Rcout<<std::endl;
    }

    itt += 1;
  }
  if(itt == maxiter){
    Rcpp::Rcout<<"IRLS doesn't converge! Try increasing maxiter. "<<std::endl;
  }

  currbeta = roundvec(currbeta, 2*eps);
  return currbeta;
}

//Function to calculate the optimized beta's for multiple values of lambda
//
//@param X,indset,Apo matrices generated by the lasso.tree function in func_sol.R, representing covariates, risk set indicator and partial order
//@param nb,m number of boundaries and number of data sets
//@param cen censoring indicator vector
//@param lambda tuning parameter
//@param w adaptive param
//@param maxiter,eps,trace parameters from the lasso.tree function
//
//@return a NumericMatrix with the first column being values of lambda, the other columns being values of beta's
//[[Rcpp::export]]
MatrixXd lasso_tree_multi(const MatrixXd & X, const MatrixXd & indset, const MatrixXd & Apo, const int nb, const int m,
                          const VectorXd & cen, const VectorXd & lambda, const VectorXd & w, const VectorXd & inibeta, const int maxiter, const double eps, const bool trace)
{
  int nlambda = lambda.size();
  int p = X.cols();
  MatrixXd all_beta(nlambda, p);

  VectorXd beta;
  //VectorXd lastbeta = inibeta;
  for(int i = 0; i < nlambda; i++)
  {
    beta = lasso_tree_single(X, indset, Apo, nb, m, cen, lambda(i), w, inibeta, maxiter, eps, false, false);
    all_beta.row(i) = beta;
    if(trace) Rcout<<"lambda = "<<lambda(i)<<", beta's: "<<beta.transpose()<<std::endl;
    //lastbeta = beta;
  }

  return all_beta;
}
