//selection index
// to compile g++ selection.cpp -o selection -O2 -larmadillo
//written by Dunia Pino del Carpio 2017

//selection index
// to compile g++ selection.cpp -o selection -O2 -larmadillo

#include <iostream>
#include <armadillo>
#include <selection.h>
#include <phenotype.h>

using namespace std;
using namespace arma;

//https://stackoverflow.com/questions/5565228/quantile-functions-in-boost-c
double inverseNormal(double prob, double mean, double sd){
        boost::math::normal_distribution<>myNormal (mean, sd);
        return quantile(myNormal, prob);
}

//from http://www.cplusplus.com/forum/beginner/62864/
// Returns the probability of x, given the distribution described by mu and sigma.
double pdf(double x, double mu, double sigma)
{
  //Constants
  static const double pi = 3.14159265; 
  return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
}


arma::mat covariance (arma::mat pheno) 
    {
   	double n = pheno.n_rows-1;
	arma::mat Pmean = mean( pheno, 0 ) ;
	arma::mat Pdiff = pheno.each_row()-Pmean;
	arma::mat phecov = Pdiff.t()*Pdiff;
	double val= (1/n);
	arma::mat phenocov = val*phecov;
	cout << "covariance_calculated" << endl;
	return phenocov;
	}		

	weights.load("weights.csv",csv_ascii);//economic values
	double alpha = atof(argv[1])/10; //proportion of selection
	double qnormval = inverseNormal(1-alpha,0,1);
	double intensity = pdf(qnormval,0,1)/alpha;//selection intensity value
	
	//double intensity = 1.4;//selection intensity value
	
vec Selection::SelectionIndex(Phenotype& phenotype,vec weights){
	arma::mat phenovalues = phenotype.Value();
    arma::mat geneticvalues = phenotype.Effect(ADDITIVE);	    
	arma::mat P = covariance(phenovalues); 
    arma::mat G = covariance(geneticvalues); 
    mat A = P.i();
    arma::mat b = A*G*v;//optimal index values
    arma::mat var_T = v.t()*G*v; //Variance of the Index
	arma::mat var_I = b.t()*P*b;  //Variance of the breeding objective
	arma::mat cov_TI = b.t()*G*v; 
	arma::mat r_TI = cov_TI / sqrt(var_T * var_I);
	arma::mat response_I = intensity*sqrt(var_I);
    arma::mat response_T = intensity*r_TI*sqrt(var_T); //response considering all traits 
    ofstream file;
  	file.open ("response_T.txt");
  	file << response_T;
  	file.close();
    arma::mat res;
	res =  b%G.each_col();
	arma::mat resum = sum(res);
	arma::mat cl = intensity/sqrt(var_I);
	arma::mat geneticgain = resum.each_col()%cl; //Expected genetic gain per trait
    ofstream file;
  	file.open ("geneticgain.txt");
  	file << geneticgain;
  	file.close();

    }

/*! \brief Selection of samples as parental candidates. */
/*! \files needed :phenotype data,optimal indexes and definition of sorting direction  */

Candidates Selection::Select(Phenotype& phenotype,arma::mat sel_index,SortDir)
	arma::mat sel_index = selection.SelectionIndex()
	arma::mat phenovalues = phenotype.Value();
    arma::mat index = phenovalues.each_row()% b.t();
    arma::mat ids=sum(index,1);
    uvec id = sort_index(ids,"descend");     
    arma::mat extract= phenos.rows(id);//extract rows of original phenotypic data sorted by id
    uvec parents = id + 1;
    ofstream file;//save file of parents id sorted by weighted phenotypic traits
  	file.open ("parents.txt");
  	file << parents;
  	file.close();
    cout << parents << endl;
	}
	       
    return 0;
}
