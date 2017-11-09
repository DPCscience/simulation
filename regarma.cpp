#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <tr1/random>
#include <iostream>
#include <armadillo>
//https://ideone.com/5gdfgO

using namespace std;
using namespace arma;

int main()
{
    std::random_device rd;
    std::default_random_engine engine( rd() );
    int table[500];
    std::uniform_int_distribution<int> distr(1, 6909);
    std::generate(std::begin(table), std::end(table), [&](){ return distr(engine); });
    std::copy(std::begin(table), std::end(table), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    return 0;
}



int main(int argc, char** argv)
	{
int NoQTL = 500         # The total is 1000 but the correlation comes from 500 only
int NoLargeEffsHD = 10 
double LargeEffsRatio = 0.5 
int NoPheno = 1
int NoSmallEffsHD = NoQTL-NoLargeEffsHD
int ScaledVarLargeEffsHD = (LargeEffsRatio/NoLargeEffsHD)/((1-LargeEffsRatio)/NoSmallEffsHD)
mat smallEffsHD = randn(NoSmallEffsHD,0,1)
mat largeEffsHD = randn(NoLargeEffsHD,0,ScaledVarLargeEffsHD)
mat MAS<-join_horiz(smallEffsHD,largeEffsHD)
int NoSNP=6909 
std::random_device rd;
std::default_random_engine engine( rd() );
int table[500];
    std::uniform_int_distribution<int> distr(1, );
    std::generate(std::begin(table), std::end(table), [&](){ return distr(engine); });
    std::copy(std::begin(table), std::end(table), std::ostream_iterator<int>(std::cout, " "));

//mat Positions = sample(NoSNP,500))
//EffectSNPs=matrix(0,nrow = NoSNP,ncol=1)



//for (i in 1:NoPheno) {
 // EffectSNPs[Positions[,i],i]=as.matrix(MAS[,i])
//}


arma::vec & y, const arma::mat & X
int n = X.n_rows, k = X.n_cols;
arma::colvec coef = arma::solve(X, y); 
arma::colvec resid = y - X*coef; 
double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
arma::colvec stderrest = 
arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
##coefficients are the marker effects