//selection index
// to compile g++ selection.cpp -o selection -O2 -larmadillo

#include <iostream>
#include <armadillo>
#include <selection.h>
#include <phenotype.h>

using namespace std;
using namespace arma;

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
	double intensity = 1.4;//selection intensity value
	
vec Selection::SelectionIndex(Phenotype& phenotype,vec weights){
	arma::mat phenovalues = phenotype.Value();
    arma::mat geneticvalues = phenotype.Effect(ADDITIVE);	    
	arma::mat P = covariance(phenos); 
    arma::mat G = covariance(TBV); 
    mat A = P.i();
    arma::mat b = A*G*v;
    
    //arma::mat Selection::b() const {
    //return b
    //	}
    }

/*! \brief Selection of samples as parental candidates. */
/*! \files needed :phenotype data,optimal indexes and definition of sorting direction  */

Selection::Select(Phenotype& phenotype,arma::mat sel_index,SortDir)
	//b = Selection::b()
	arma::mat sel_index = selection.SelectionIndex()
	arma::mat phenovalues = phenotype.Value();
    arma::mat index = phenovalues.each_row()% b.t();
    arma::mat ids=sum(index,1);
    uvec id = sort_index(ids,"descend");     
    cout << id << endl;
	}
	    
    //arma::mat extract= phenos.rows(parents);//we can extract rows of original data with index
    
    return 0;
}

	//arma::mat var_T = v.t()*G*v; 
	//arma::mat var_I = b.t()*P*b;
	//arma::mat cov_TI = b.t()*G*v; 
	//arma::mat r_TI = cov_TI / sqrt(var_T * var_I);
	//arma::mat response_I = intensity*sqrt(var_I);
    //arma::mat response_T = intensity*r_TI*sqrt(var_T);
    //arma::mat res;
	//res =  b%G.each_col();
	//arma::mat resum = sum(res);//its ok until here compared to R
	//arma::mat cl = intensity/sqrt(var_I);
	//arma::mat res1 = resum.each_col()%cl;//missing the multiplication with intensity
