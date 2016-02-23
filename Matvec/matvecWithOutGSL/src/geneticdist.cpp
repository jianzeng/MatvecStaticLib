/*
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Library General Public
 License as published by the Free Software Foundation; either
 version 2 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Library General Public License for more details.
 
 You should have received a copy of the GNU Library General Public
 License along with this library; if not, write to the Free
 Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 MA 02111-1307, USA 
 */

#include <iostream>
#include "geneticdist.h"
#include <cstdarg>


namespace matvec {

/////////////////////////// LocusStruct class ///////////////////////

void LocusStruct::copyfrom(const LocusStruct& A)
{
   resize(A.numallele);
   qtl_ml = A.qtl_ml;
   allele_freq = A.allele_freq;
   genotypic_val_mat = A.genotypic_val_mat;
}

void LocusStruct::resize(const unsigned na)
{
   if (numallele == na) return;
   numallele = na;
   allele_freq.resize(numallele);
   genotypic_val_mat.resize(numallele,numallele);
}

/////////////////////////// ChromStruct class ///////////////////////

void ChromStruct::copyfrom(const ChromStruct& A)
{
   if (this == &A) return;
   resize(A.numloci);
   recomb_rate_mat = A.recomb_rate_mat;
}

void ChromStruct::resize(const unsigned nl)
{
   if (numloci == nl) return;
   numloci = nl;
   locus.resize(numloci);
   recomb_rate_mat.resize(numloci,numloci,0.5);
}

void ChromStruct::release(void)
{locus.resize(0);numloci = 0;}

void ChromStruct::calcLocusMask(void){
	// Authors: L. Radu Totir and Rohan L. Fernando
	// (December, 2004) 
	// Contributors: 
  locusMask.resize(numloci);
	for(unsigned i=0;i<numloci; i++){
		locusMask[i].index = i;
		locusMask[i].position = locus[i].distance;
	}
	sort(locusMask.begin(),locusMask.end());
//	for(unsigned i=0;i<numloci;i++){
//	   cout << locusMask[i].position << " for " << locusMask[i].index << endl;
//	}	
}
/////////////////////////// GeneticDist class ///////////////////////

GeneticDist::GeneticDist(const unsigned nc,...)
{
   numchrom = 0;
   chromosome = 0;
   numtrait = 1;      // temporary value
   numMarkerLoci = 0;
   resize(nc);
   va_list param_pt;
   va_start(param_pt,nc);
   unsigned nl;
   for (unsigned i=0; i<nc; i++) {
      nl = va_arg(param_pt,unsigned);
      chromosome[i].resize(nl);
   }
   va_end(param_pt);
   strcpy(distname,"GeneticDist");
}

void GeneticDist::copyfrom(const GeneticDist& A)
{
   if (this == &A) return;
   resize(A.numchrom);
   for (unsigned i=0; i<numchrom; i++) chromosome[i] = A.chromosome[i];
   strcpy(distname,A.distname);

   numtrait = A.numtrait;  // Occurs in old version
   numMarkerLoci = A.numMarkerLoci;  //LRT
}

void GeneticDist::nchrom(const unsigned nc)
{
   if (numchrom == nc) return;
   numchrom = nc;
   if (chromosome) {
     delete [] chromosome;
     chromosome=0;
   }
   if(numchrom>0){
     chromosome = new ChromStruct[numchrom];
   }
   else {
     chromosome = 0;
   }
   check_ptr(chromosome);
}

void GeneticDist::nloci(const unsigned nl0,...)
{
   chromosome[0].resize(nl0);
   va_list param_pt;
   va_start(param_pt,nl0);
   unsigned nl;
   for (int i=1;i<numchrom;i++) {
      nl = va_arg(param_pt,unsigned);
      chromosome[i].resize(nl);
   }
   va_end(param_pt);
}

//RLF
void GeneticDist::putColmNames(const unsigned c,const unsigned l,
                                char nm1[],char nm2[])
{
   if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::locus(): bad args");
   LocusStruct *L = &(chromosome[c-1].locus[l-1]);
   L->nameOfcol1 = nm1;
   L->nameOfcol2 = nm2;
}

void GeneticDist::putColmNames(const unsigned c,const unsigned l,
                                string nm1, string nm2)
{
   if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::locus(): bad args");
   LocusStruct *L = &(chromosome[c-1].locus[l-1]);
   L->nameOfcol1 = nm1.c_str();
   L->nameOfcol2 = nm2.c_str();
}
//RLF

void GeneticDist::genotypic_val(const unsigned c,const unsigned l,
                                const double* v)
{
   ///////////////////////////////////////////
   // requiements
   // v must have enough elements
   ///////////////////////////////////////////
   if (c < 1 || c>numchrom || l < 1 || l>chromosome[c-1].nloci()) throw exception("GeneticDist::genotypic_val(): bad args");
   unsigned i,j,k,na,cc,ll;
   cc = c-1; ll = l-1;
   na = chromosome[cc].locus[ll].nallele();
   double **me = chromosome[cc].locus[ll].genotypic_val_mat.begin();
   for (k=0,i=0; i<na; i++) for (j=0; j<=i; j++) {
      me[i][j] = v[k++];
      me[j][i] =  me[i][j];
   }
}

void GeneticDist::genotypic_val(const unsigned c,const unsigned l,
                                 const double v0,...)
{
   if (c < 1 || c>numchrom || l < 1 || l>chromosome[c-1].nloci()) throw exception("GeneticDist::genotypic_val(): bad args");
   unsigned i,j,na,cc,ll;
   cc = c-1; ll = l-1;
   na = chromosome[cc].locus[ll].nallele();
   double **me = chromosome[cc].locus[ll].genotypic_val_mat.begin();
   me[0][0] = v0;
   va_list ppt;
   va_start(ppt,v0);
   for (j=1; j<na; j++) me[j][0] = me[0][j] = va_arg(ppt,double);
   for (i=1;i<na;i++) for (j=i;j<na;j++) me[j][i]=me[i][j]= va_arg(ppt,double);
   va_end(ppt);
}

void GeneticDist::release(void)
{
   if (chromosome) {
      delete [] chromosome; chromosome = 0; numchrom = 0;
   }
}

double GeneticDist::recomb_rate(const unsigned c,const unsigned li,
                                 const unsigned lj)
{
   if (c==0 || li==0) throw exception("GeneticDist::recomb_rate(): bad args");
   double d=0.0;
   unsigned k, ii=li-1;
   for (k=li; k<lj; k++) d += chromosome[c-1].recomb_rate_mat[ii][k];
   return 0.5*tanh(2.0*d);
}

void GeneticDist::recomb_rate(const unsigned c,const unsigned li,
                              const unsigned lj, const double r)
{
   if (c <= 0) throw exception("GeneticDist.recomb_rate(): out of range");
   int n = chromosome[c-1].nloci();
   if (li<1 || li>n || lj<1 || lj>n) throw exception("GeneticDist::recomb_rate(): out of range");
   chromosome[c-1].recomb_rate_mat[li-1][lj-1] = r;
}

///////////////////  UnknownDist class /////////////////////////////

UnknownDist::UnknownDist(const Vector<double> mu, const doubleMatrix sigma)
{
   mu_vec = mu;
   dim = mu_vec.size();
   if (sigma.psd()) {
      if (dim == sigma.num_rows()) {
         var_mat = sigma;
      } else {
          warning(" UnknownDist(u,V): u and V incompatible");
      }
   } else {
      warning(" UnknownDist(u,V): V must be psd ");
   }
   strcpy(distname,"Unknown");
}

UnknownDist::UnknownDist(const double mu, const double sigma)
{
   resize(1);
   mu_vec[0] = mu;
   var_mat[0][0] = sigma;
   strcpy(distname,"Unknown");
}

void UnknownDist::copyfrom(const UnknownDist& A)
{
   if (this == &A) return;
   dim = A.dim;
   mu_vec = A.mu_vec;
   var_mat = A.var_mat;
   strcpy(distname,A.distname);
}

void UnknownDist::resize(const unsigned n)
{
   if (dim == n) return;
   dim = n;
   mu_vec.resize(dim);
   var_mat.resize(dim,dim);
}

//BRS
void GeneticDist::put_distance(const unsigned c,const unsigned l, double dist)
{
   if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::put_distance() : bad args");
   LocusStruct *L = &(chromosome[c-1].locus[l-1]);
   L->distance = dist;
}

double GeneticDist::get_distance(const unsigned c,const unsigned l)
{
   if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception(" GeneticDist.get_distance(): bad args");
   LocusStruct *L = &(chromosome[c-1].locus[l-1]);
   return (L->distance);
}

void GeneticDist::multi_loci(int num_loci)
{

   resize(1);
//   cout << "in GeneticDist::multi_loci: numtrait = " << numtrait<< endl;
   chromosome[0].resize(num_loci);
}

int GeneticDist::display(void) {
  if (chromosome) {
    std::cout << "GeneticDist Object: nloci = " << chromosome[0].nloci() << std::endl;
             std::cout << "numtrait = " << numtrait <<  std::endl;
  }
  else {
    std::cout << "GeneticDist Object: nloci not initialized yet " << std::endl;
    std::cout << "numtrait = " << numtrait <<  std::endl;
  }
  return 1;
}

//BRS

void GeneticDist::locus(const unsigned c, const unsigned l, string type, const unsigned na, ...)
{
	va_list param_pt;
	va_start(param_pt,na);

	if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::locus(): bad args");
	LocusStruct *L = &(chromosome[c-1].locus[l-1]);
	
	//RLF modified block
	L->gnodeType = "allelic";
	if (type == "qtl") {
		L->qtl_ml='q';
		L->genotypic_val_mat.resize(na,na);
	}
	// next else if added by LRT
	else if(type == "recessiveLocus") { 
		L->qtl_ml='r';
	}
	else if(type == "marker") { 
		L->qtl_ml = 'm';
		numMarkerLoci++;
	}
	else {
		throw("GeneticDist::locus(): locus type not recognized \n");
		
  }
	L->nallele(na);
	L->allele_freq.resize(na);
	//RLF modified block
	double sumf = 0.0;
	for (int i=0;i<na; i++) {
		double q = va_arg(param_pt,double);
		sumf += L->allele_freq[i] = q;
	}
	va_end(param_pt);
	if (fabs(sumf-1.0) > 1.0e-10) throw exception("GeneticDist::locus(): allele_freq does not sum to 1.0");
}



void GeneticDist::locus(const unsigned c,const unsigned l,string type,
						const unsigned na, Vector<double> init_freqs)
{

  if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::locus(): bad args");
  LocusStruct *L = &(chromosome[c-1].locus[l-1]);
  
  //RLF modified block
  L->gnodeType = "allelic";
  if (type=="qtl") {
    L->qtl_ml='q';
    L->genotypic_val_mat.resize(na,na);
  }
  // next else if added by LRT
  else if(type=="recessiveLocus") { 
    L->qtl_ml='r';
  }
  else if (type=="marker")  {
    L->qtl_ml = 'm';
	numMarkerLoci++;
  }
  else {
	throw("GeneticDist::locus(): locus type not recognized \n");
  }
  L->nallele(na);
  L->allele_freq.resize(na);
  //RLF modified block
  //   va_list param_pt;
  //   va_start(param_pt,na);
  double sumf = 0.0, temp=0.0;
  for (int i=0;i<na; i++) {
    temp=init_freqs[i];
    sumf += temp;
    L->allele_freq[i] = temp;
    //	   std::cout << i << " freq " << temp << " Sumf " << sumf << std::endl;
  }
  //   va_end(param_pt);
  if (fabs(sumf-1.0) > 1.0e-10) throw exception("GeneticDist::locus(): allele_freq does not sum to 1.0");
}

void GeneticDist::putGnodeType(const unsigned c,const unsigned l, string type)
{
   if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::putGnodeType() : bad args");
   LocusStruct *L = &(chromosome[c-1].locus[l-1]);
   L->gnodeType = type;
}

void GeneticDist::calcProbs(const unsigned c,const unsigned l, string type)
{
	if (c < 1 || c>numchrom || l < 1 || l> chromosome[c-1].nloci()) throw exception("GeneticDist::calcProbs() : bad args");
	LocusStruct *L = &(chromosome[c-1].locus[l-1]);
	if (type == "origin"){
		L->originProbs = true;
	}
	else if (type == "state"){
		L->stateProbs = true;
	}
	else {
		throw exception("GeneticDist::calcProbs() : bad arg for type");
	}	
}


} ////////////// end of namespace matvec
