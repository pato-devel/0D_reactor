/*
 * FiniteRateSurface.cpp
 * ---------------------
 *    file created:   January   16, 2013
 *    last modified:  August    24, 2013
 *    author:         Matthew MacLean
 *                    maclean@cubrc.org
 *    purpose:        implementation of FiniteRateSurface class functions
 *
 * notable revisions:
 * -----------------
 *    N/A
 *
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *
 *  This file is a port of the Finite Rate Surface Chemistry Model (FRM) code 
 *  for use as a Gas/Surface Interaction (GSI) module.
 *
 *  Original code developed by:
 *     Jochen Marschall (jochen.marschall@sri.com)
 *     Matthew MacLean  (maclean@cubrc.org)
 *
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *
 */
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <cmath>
#include "FiniteRateSurface.hxx"
#include "lumatrix.hxx"



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                  LOCAL PREPROCESSOR DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


#define INPUT_FILE_BUFFER_SIZE 1024

#define INPUT_STRTOK_TOKENS " \t',"

#define PI 3.1415926536

#define R 8.31451

#define AV 6.0221367e+23

#define Pref 100000.0

#define H 6.62618e-34


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                 LOCAL (PRIVATE) VARIABLE DECLARATIONS                  ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###               LOCAL (PRIVATE) FUNCTION DECLARATIONS                    ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/* 
 * name:      get_gibbs_data
 *
 * arguments: data     = Lewis coefficient data array (NLEWIS_COEFFS=10,NLEWIS_FITS=6,nspecies)
 *            t        = temperature
 *            irange   = range index (parameter 2)
 *            ispecies = species index (parameter 3)
 *            h        = computed enthalpy (return)
 *            s        = computed entropy (return)
 *            cp       = computed computed specific heat (return)
 *
 * return:    nothing
 *
 * purpose:   compute gibbs data from Lewis arrays (form 1)
 *
 */ 
void get_gibbs_data(Array3D<double> &data, double t, int irange, int ispecies, double *h, double *s, double *cp);



/* 
 * name:      get_gibbs_data_2
 *
 * arguments: data     = Lewis coefficient data array (NLEWIS_COEFFS=10,NLEWIS_FITS=6,nspecies)
 *            trange   = Lewis range data array (min/max=0 or 1,NLEWIS_FITS=6,nspecies)
 *            t        = temperature
 *            ispecies = species index (parameter 3)
 *            h        = computed enthalpy (return)
 *            s        = computed entropy (return)
 *            dh       = computed enthalpy derivative (return)
 *            ds       = computed entropy derivative (return)
 *
 * return:    nothing
 *
 * purpose:   compute gibbs data from Lewis arrays with derivative information (form 2)
 *
 */ 
void get_gibbs_data_2(Array3D<double> &data, Array3D<double> &trange, double t, int ispecies, double *h, double *s, double *dh, double *ds);



/* 
 * name:      matchelemname
 *
 * arguments: nelements   = reference to # elements counter (incremented if no match found)
 *            elem_names  = string array of existing names
 *            name        = C-string to search for
 *            name_length = length of C-string (may be a substring)
 *
 * return:    index of match
 *
 * purpose:   match element name to working list for EQ solver
 *
 */ 
int matchelemname(int &nelements, string *elem_names, const char *name, const int name_length);



/* 
 * name:      simplex
 *
 * arguments: ebal     = matrix of elemental constraints
 *            gibbs    = list of Gibbs energies
 *            ibasis   = list of basis vector for each row
 *
 * return:    0=no error
 *
 * purpose:   perform simplex method on elemental contrainsts to obtain a basis 
 *            set with minimum Gibbs energy choices
 *
 */ 
int simplex(Array2D<double> &ebal, Array1D<double> &gibbs, Array1D<int> &ibasis);




/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                     PUBLIC FUNCTION DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */




/* @@@@@@@@@@@@@@@@@@@@@@@@@ FiniteRateSurfaceCell @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


/* ------------------------- FiniteRateSurfaceCell ------------------------------- */
/* 
 * constructor
 *
 */
FiniteRateSurfaceCell::FiniteRateSurfaceCell()
{
	model = NULL;
	tmpr = 300.0;
	pres = 100000.0;
}
/* ------------------------------------------------------------------------------- */



/* ------------------------ ~FiniteRateSurfaceCell ------------------------------- */
/* 
 * destructor.
 *
 */
FiniteRateSurfaceCell::~FiniteRateSurfaceCell()
{

}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- initialize_cell ---------------------------------- */
/* 
 * initialize surface cell with a model
 *
 */
void FiniteRateSurfaceCell::initialize_cell(FiniteRateSurfaceModel *frsmodel)
{
	int nspecies, neqns;

	model = frsmodel;

	if (model != NULL)	{

		//size arrays for this model
		nspecies = model->get_num_species();
		delta_conc.resize(nspecies);
		neqns = nspecies + 1;
		concentration.resize(nspecies);
		srhs.resize(neqns);
		sjacobian.resize(neqns,neqns);

		//complete initialization
		srhs.initialize(0.0);
		sjacobian.initialize(0.0);
		delta_conc.initialize(0.0);
		model->initialize_cell(*this);
	}
}
/* ------------------------------------------------------------------------------- */



/* --------------------------- compute_frm_system -------------------------------- */
/* 
 * initialize surface cell with a model
 *
 */
void FiniteRateSurfaceCell::compute_frm_system(double tmpr)
{
	int i;
	
	//initialize here
	delta_conc.initialize(0.0);
	this->tmpr = tmpr;
	
	//populate if we've set a model
	if (model) model->compute_frm_system(*this, tmpr);
	else	{
		srhs.initialize(0.0);
		sjacobian.initialize(0.0);
		for (i=0; i<get_num_eqns(); i++)
			sjacobian(i,i) = 1.0;
	}
}
/* ------------------------------------------------------------------------------- */



/* ------------------------- compute_equilibrium_system -------------------------- */
/* 
 * compute equilibrium composition
 *
 */
double FiniteRateSurfaceCell::compute_equilibrium_system(double tmpr, int igastype, ostream *stream_ptr)
{
	this->tmpr = tmpr;
	this->pres = 0.0;
	for (int i=0; i<get_num_gas(); i++)
		this->pres += concentration(i);
	this->pres = (this->pres)*8.31451*tmpr;
	
	if (model)	return model->compute_equil_system(*this, igastype, stream_ptr);
	else return 0.0;
}
/* ------------------------------------------------------------------------------- */



/* --------------------------- compute_qss_system -------------------------------- */
/* 
 * compute QSS composition
 *
 */
double FiniteRateSurfaceCell::compute_qss_system(double tmpr)
{
	
	if (model)	return model->compute_qss_system(*this, tmpr);
	else return 0.0;
}
/* ------------------------------------------------------------------------------- */



/* @@@@@@@@@@@@@@@@@@@@@@@@ FiniteRateSurfaceModel @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


/* ------------------------ FiniteRateSurfaceModel ------------------------------- */
/* 
 * constructor
 *
 */
FiniteRateSurfaceModel::FiniteRateSurfaceModel()
{
	title[0] = '\0';
	cfrm_error_string[0] = '\0';
	nphases = 0;
	nspecies = 0;
	
	igstart = 0;
	igend   = 0;
	isstart = 0;
	isend   = 0;
	ibstart = 0;
	ibend   = 0;
	
	initialized = false;
}
/* ------------------------------------------------------------------------------- */



/* ----------------------- ~FiniteRateSurfaceModel ------------------------------- */
/* 
 * destructor.
 *
 */
FiniteRateSurfaceModel::~FiniteRateSurfaceModel()
{
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- initialize_cell ---------------------------------- */
/* 
 * initialize surface cell concentrations
 *
 */
void FiniteRateSurfaceModel::initialize_cell(FiniteRateSurfaceCell& cell)
{
	Array1D<double> &concentrations = cell.get_concentration();
	int nsurf, nsa, nbulk, j, k;
	
	//double check sizing is correct
	if (concentrations.getRows() != nspecies)	return;
	//initialize all to default
	concentrations.initialize(1.0e-15);
	//seek through to set all empty site species to site density
	k = ngps - 1;
	for (nsurf=0; nsurf<nsp; nsurf++)	{
		for (nsa=0; nsa<nspas(nsurf); nsa++)	{
			for (j=0; j<nspass(nsurf,nsa); j++)	{
				k++;
				if (j==0) concentrations(k) = sdenas(nsurf,nsa);
			}
		}
	}
	//set all bulk concentration to mole fraction
	for (nbulk=0; nbulk<nbp; nbulk++)	{
		for (j=0; j<nbps(nbulk); j++)	{
			k++;
			concentrations(k) = bmf(nbulk,j);
		}
	}
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- compute_frm_system ------------------------------- */
/* 
 * compute finite rate system with Jacobians
 *
 */
void FiniteRateSurfaceModel::compute_frm_system(FiniteRateSurfaceCell& cell, double tmpr)
{
	Array1D<double> &srhs = cell.get_rhs();
	Array2D<double> &sjacobian = cell.get_jacobian();
	Array1D<double> &concentration = cell.get_concentration();

	int ispec,ireact,iphase,ibulk,index,j,ier,k;
	double  tmp, prate, rflux, pf, pb, dFdtg, dFdtw, mdotc;
	Array1D<double> frates(isrt), brates(isrt), eqa(isrt), eqc(isrt),
			    df_dtg(isrt), df_dtw(isrt), deqa_dtg(isrt), deqa_dtw(isrt),
			    deqc_dtg(isrt), deqc_dtw(isrt), db_dtg(isrt), db_dtw(isrt);
	Array1D<double> gibbs(nspecies), dg_dtw(nspecies), dg_dtg(nspecies);
	
	/* --- zero cell arrays --- */
	srhs.initialize(0.0);
	sjacobian.initialize(0.0);
	
	/* --- initialize --- */
	frates.initialize(0.0);
	brates.initialize(0.0);
	eqa.initialize(0.0);
	eqc.initialize(0.0);
	df_dtg.initialize(0.0);
	df_dtw.initialize(0.0);
	deqa_dtg.initialize(0.0);
	deqa_dtw.initialize(0.0);
	deqc_dtg.initialize(0.0);
	deqc_dtw.initialize(0.0);
	db_dtg.initialize(0.0);
	db_dtw.initialize(0.0);
	gibbs.initialize(0.0);
	dg_dtw.initialize(0.0);
	dg_dtg.initialize(0.0);
	
	dFdtg = 0.0;
	dFdtw = 0.0;
	
	/* --- compute forward rates --- */
	for (ireact=0; ireact<isrt; ireact++)	{
		
		// identify gas phase participant
		k = RR(ireact,0);
		if (istype(ireact) == 4) k = RP(ireact,0);
		//compute forward rate
		frates(ireact) = srate(tmpr, tmpr, wmass(k), sdensity(ireact), sitep(ireact), istype(ireact),
						Cf(ireact), eta(ireact), Eact(ireact), &(df_dtg(ireact)), &(df_dtw(ireact)));
	}
	
	/* --- obtain Gibbs energies for all species --- */
	compute_gibbs_energies(gibbs, dg_dtg, dg_dtw, tmpr, tmpr, frates, df_dtg, df_dtw);

	/* --- calculate activity based equilibrium constants and derivatives --- */
	compute_eqcon(eqc, deqc_dtg, deqc_dtw, eqa, deqa_dtg, deqa_dtw, tmpr, tmpr, gibbs, dg_dtg, dg_dtw);

	/* --- calculate backward rates & derivatives --- */
	for (ireact=0; ireact<isrt; ireact++)	{
		brates(ireact) = frates(ireact) / eqc(ireact);
		db_dtg(ireact) = (df_dtg(ireact) - brates(ireact)*deqc_dtg(ireact))/eqc(ireact);
		db_dtw(ireact) = (df_dtw(ireact) - brates(ireact)*deqc_dtw(ireact))/eqc(ireact);
	}
	
	/* --- calculate reaction flux for each reaction & derivatives --- */
	for (ireact=0; ireact<isrt; ireact++)	{
	
		//sum forward and backward concentration product term
		pf = 1.0;
		pb = 1.0;
		for (ispec=0; ispec<nspecies; ispec++)	{
			pf = pf * pow(concentration(ispec), vr(ireact,ispec));
			pb = pb * pow(concentration(ispec), vp(ireact,ispec));
		}
		//account for rates turned on and compute reaction flux
		pf = pf * ((double)(isrfon(ireact))) * sfraction(ireact);
		pb = pb * ((double)(isrbon(ireact))) * sfraction(ireact);
		rflux = frates(ireact)*pf - brates(ireact)*pb;
		//compute the contribution to all species production and derivatives
		for (ispec=0; ispec<nspecies; ispec++)	{
			prate = rflux * v(ireact,ispec);

			//explicit contribution
			srhs(ispec) += prate;

			//derivatives w.r.t. concentration (concentration must not be exactly zero!)
			for (j=0; j<nspecies; j++)
				sjacobian(ispec,j) += v(ireact,ispec) * (frates(ireact)*pf*vr(ireact,j) - brates(ireact)*pb*vp(ireact,j)) / concentration(j);

			//derivative w.r.t. temperature
			dFdtg = v(ireact,ispec) * (pf*df_dtg(ireact)-pb*db_dtg(ireact));
			dFdtw = v(ireact,ispec) * (pf*df_dtw(ireact)-pb*db_dtw(ireact));
			sjacobian(ispec,nspecies) += dFdtg + dFdtw;
		
		}
	
	}

	/* --- add pyrolysis blowing term in --- */
	if (nblwflag == 1)	{

		if (nebc >= 1)	{
		
			//sum up bulk blowing rates to add to pyrolysis
			for (ibulk=ibstart; ibulk<ibend; ibulk++)	{
				for (ispec=0; ispec<ngps; ispec++)
					srhs(ispec) -= pyroadd(ispec)*wmass(ibulk)*srhs(ibulk);
				for (j=0; j<=nspecies; j++)
					sjacobian(ispec,j) -= pyroadd(ispec)*wmass(ibulk)*sjacobian(ibulk,j);
			}
		
		} else {
		
			//explicit mass flow addition
			for (ispec=0; ispec<ngps; ispec++)
				srhs(ispec) += pyroadd(ispec);
		
		}

	}
	
	/* --- diagonalize temperature row --- */
	sjacobian(nspecies,nspecies) = 1.0;

}
/* ------------------------------------------------------------------------------- */



/* --------------------------- compute_qss_system -------------------------------- */
/* 
 * compute QSS composition
 *
 */
double FiniteRateSurfaceModel::compute_qss_system(FiniteRateSurfaceCell& cell, double tmpr)
{
	Array1D<double> &delta_conc = cell.get_delta_conc();
	Array1D<double> &srhs = cell.get_rhs();
	Array2D<double> &sjacobian = cell.get_jacobian();
	Array1D<double> &concentration = cell.get_concentration();

	int j, k, irow, icol, nsurf, nsa, num_gas, num_surf, num_sites, neqns, istep;
	double resid;
	Array2D<double> qsslhs;
	Array1D<double> qssrhs;
	
	//obtain # of equations
	num_gas   = get_num_gas();
	num_surf  = get_num_surf();
	num_sites = get_num_sitesets();
	neqns = num_surf + num_sites;
	//resize matrices
	qssrhs.resize(neqns);
	qsslhs.resize(neqns,neqns);
	//reinitialize surface concentrations
	k = ngps - 1;
	for (nsurf=0; nsurf<nsp; nsurf++)	{
		for (nsa=0; nsa<nspas(nsurf); nsa++)	{
			for (j=0; j<nspass(nsurf,nsa); j++)	{
				k++;
				if (j==0) concentration(k) = sdenas(nsurf,nsa);
				else      concentration(k) = 1.0e-15;
			}
		}
	}
	//perform Newton iterations
	istep = 0;
	resid = 1.0;
	while ((istep<25)&&(resid>1.0e-30))	{
	
		//initialize matrices
		qssrhs.initialize(0.0);
		qsslhs.initialize(0.0);

		//compute rates
		compute_frm_system(cell, tmpr);

		//build matrices
		for (k=0; k<num_surf; k++)	{
			qssrhs(k) = -srhs(num_gas+k);
			for (j=0; j<num_surf; j++)
				qsslhs(k,j) = sjacobian(num_gas+k,num_gas+j);
		}
		irow = num_surf-1;
		icol = -1;
		for (nsurf=0; nsurf<nsp; nsurf++)	{
			for (nsa=0; nsa<nspas(nsurf); nsa++)	{
				irow++;
				qssrhs(irow) = sdenas(nsurf,nsa);
				for (k=0; k<nspass(nsurf,nsa); k++)	{
					icol++;
					qsslhs(irow,icol) = 1.0;
					qsslhs(icol,irow) = 1.0;
					qssrhs(irow) = qssrhs(irow) - concentration(num_gas+icol);
				}
			}
		}
		
		//solve matrices
		lusolve(qsslhs, qssrhs, qssrhs);

		//update variables
		resid = 0.0;
		for (k=0; k<num_surf; k++)	{
			concentration(num_gas+k) += qssrhs(k);
			if (concentration(num_gas+k) < 1.0e-20) concentration(num_gas+k) = 1.0e-20;
			resid = resid + qssrhs(k)*qssrhs(k);
		}

		istep++;
	}

	return resid;
}
/* ------------------------------------------------------------------------------- */



/* -------------------------- compute_equil_system ------------------------------- */
/* 
 * compute equilibrium distribution of a system
 *
 */
double FiniteRateSurfaceModel::compute_equil_system(FiniteRateSurfaceCell& cell, int igastype, ostream *stream_ptr)
{
	Array1D<double> &concentration = cell.get_concentration();
	Array1D<double> initconc, frates(isrt), df_dtg(isrt), df_dtw(isrt), 
	                gibbs(nspecies), dg_dtw(nspecies), dg_dtg(nspecies),
			    sref(nspecies);
	Array1D<double> gtot, rhs, uval, bi0_array;
	Array2D<double> ebal, rjacobian;
	Array1D<int>    inc_species(nspecies), ibasis;
	char buffer[INPUT_FILE_BUFFER_SIZE], unixnewline = 0x0A;  //LineFeed
	int i, j, k, ireact, ispec, nfrozen, nvol, idx, neqns, icount;
	double  tmpr, volume, tterm, rtmp, pres, residual, tmp;

	/* --- set up calculation --- */
	
	//save the initial concentrations
	initconc  = concentration;
	tmpr = cell.get_temperature();
	pres = cell.get_pressure();

	//compute Gibbs energies for all species
	frates.initialize(0.0);
	df_dtg.initialize(0.0);
	df_dtw.initialize(0.0);
	for (ireact=0; ireact<isrt; ireact++)	{
		
		// identify gas phase participant
		k = RR(ireact,0);
		if (istype(ireact) == 4) k = RP(ireact,0);
		//compute forward rate
		frates(ireact) = srate(tmpr, tmpr, wmass(k), sdensity(ireact), sitep(ireact), istype(ireact),
						Cf(ireact), eta(ireact), Eact(ireact), &(df_dtg(ireact)), &(df_dtw(ireact)));

	}
	
	/* --- obtain Gibbs energies for all species --- */
	
	//obtain Gibbs energies
	gibbs.initialize(0.0);
	dg_dtw.initialize(0.0);
	dg_dtg.initialize(0.0);
	compute_gibbs_energies(gibbs, dg_dtg, dg_dtw, tmpr, tmpr, frates, df_dtg, df_dtw);

	//Gibbs energies must be non-dimensional
	for (ispec=0; ispec<nspecies; ispec++)
		gibbs(ispec) = gibbs(ispec) / (R*tmpr);
	
	//get reference site densities for the surface species
	sref.initialize(1.0);  //make ln(sref)=0.0 by default
	for (k=isstart; k<isend; k++)
		sref(k) = sden(kphase(k)-1);

	/* --- set range variables --- */
	
	nfrozen = get_num_bulk();
	if (igastype==1) nfrozen += get_num_gas();
	
	nvol = 0;
	if (igastype==2) nvol++;
	volume = 1.0;
	
	inc_species.initialize(1);
	
	tterm = log(8.314*tmpr/1.0e5);

	/* --- eliminate redundancies --- */
	
	//remove frozen redundancies using simplex method
	if (nfrozen > 0)	{
	
		ebal.resize(nelements,nfrozen);
		ibasis.resize(nelements);
		
		gtot.resize(nfrozen);
		gtot.initialize(0.0);

		//catch frozen redundancies
		idx = 0;
		if (igastype==1)	{
			for (j=igstart; j<igend; j++)	{
				inc_species(j) = 0;
				rtmp = initconc(j);
				if (rtmp<1e-20) rtmp = 1e-20;
				gtot(idx) = gibbs(j) + log(rtmp) + tterm;
				for (i=0; i<nelements; i++)
					ebal(i,idx) = etable(j,i);
				idx++;
			}
		}
		for (j=ibstart; j<ibend; j++) {
			inc_species(j) = 0;
			gtot(idx) = gibbs(j);
			for (i=0; i<nelements; i++)
				ebal(i,idx) = etable(j,i);
			idx++;
		}
		
		//simplex method to select species
		simplex(ebal, gtot, ibasis);
		
		//set active species
		for (i=0; i<nelements; i++)	{
			j = ibasis(i);
			if (j >= 0)	{
				if (igastype == 1)	{
					if (j >= igend) j = j - igend + isend;
				} else
					j += isend;
				inc_species(j) = 1;
			}
		}
		
		//write output to stream if requested
		if (stream_ptr)	{
			//set active species
			(*stream_ptr) << "     Gibbs Data of Species:" << unixnewline;
			for (j=0; j<nspecies; j++)	{
				if (j<igend)	{
					rtmp = initconc(j);
					if (rtmp<1e-20) rtmp = 1e-20;
					rtmp = gibbs(j) + log(rtmp) + tterm;
				} else
					rtmp = gibbs(j);
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "       %-6s %14.4e %14.4e",species(j).c_str(),gibbs(j),rtmp);
				(*stream_ptr) << buffer << unixnewline;
			}
			//write active species out
			(*stream_ptr) << "     Species Considered for System:" << unixnewline << "       ";
			for (j=0; j<nspecies; j++)	{
				if (inc_species(j) > 0)	{
					snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-6s",species(j).c_str());
					(*stream_ptr) << buffer;
				}
			}
			(*stream_ptr) << unixnewline;
		}
	}


	/* --- allocate matrices --- */
	neqns = nspecies + nelements + nvol;
	bi0_array.resize(nelements);
	uval.resize(neqns);
	rhs.resize(neqns);
	rjacobian.resize(neqns,neqns);	


	/* --- initialize field --- */
	bi0_array.initialize(0.0);
	for (i=0; i<nelements; i++)	{
		for (j=0; j<nspecies; j++)
			bi0_array(i) += etable(j,i)*initconc(j);
	}
	
	uval.initialize(0.0);
	for (j=0; j<nspecies; j++)	{
		uval(j) = initconc(j);
		if (uval(j) < 1.0e-7) uval(j) = 1.0e-7;
	}

	/* --- MAIN LOOP --- */
	for (icount=0; icount<100; icount=icount+1)	{
		
		rhs.initialize(0.0);
		rjacobian.initialize(0.0);	
	
		// gas phase Gibbs equations
		for (j=igstart; j<igend; j++)	{
		
			if (igastype != 1) rjacobian(j,j) = 1.0;
			
			for (i=0; i<nelements; i++)
				rjacobian(j,nspecies+i) = etable(j,i);
			
			rhs(j) = -(gibbs(j) + log(uval(j)) + tterm);
		
		}
		
		// surface phase Gibbs equations
		for (j=isstart; j<isend; j++)	{
		
			rjacobian(j,j) = 1.0;
		
			for (i=0; i<nelements; i++)
				rjacobian(j,nspecies+i) = etable(j,i);
			
			rhs(j) = -(gibbs(j) + log(uval(j)) + log(sref(j)));
		
		}
		
		// bulk phase Gibbs equations
		for (j=ibstart; j<ibend; j++)	{
		
			for (i=0; i<nelements; i++)
				rjacobian(j,nspecies+i) = etable(j,i);
			
			rhs(j) = -(gibbs(j) + log(initconc(j)));
		
		}
		
		// elemental constraints
		for (i=0; i<nelements; i++)	{
		
			rhs(nspecies+i) = bi0_array(i);
			
			for (j=igstart; j<igend; j++)	{
				tmp = volume*etable(j,i)*uval(j);
				rjacobian(nspecies+i,j) = tmp;
				rhs(nspecies+i) -= tmp;
			}
			
			for (j=isstart; j<ibend; j++)	{
				tmp = etable(j,i)*uval(j);
				rjacobian(nspecies+i,j) = tmp;
				rhs(nspecies+i) -= tmp;
			}
			if (igastype==2)	{
				for (j=igstart; j<igend; j++)
					rjacobian(nspecies+i,neqns-1) += etable(j,i)*uval(j);
			}
		}
		
		// volume constraint (if constant pressure)
		if (igastype==2)	{
			rhs(neqns-1) = pres/(8.314*tmpr);
			for (j=igstart; j<igend; j++)	{
				rhs(neqns-1) -= uval(j);
				rjacobian(neqns-1,j) = uval(j);
			}
		}
		
		// force out non-participating species
		for (j=0; j<nspecies; j++)	{
			if (inc_species(j) == 0)	{
				for (i=0; i<neqns; i++)
					rjacobian(j,i) = 0.0;
				rjacobian(j,j) = 1.0;
				rhs(j) = 0.0;
			}
		}
/*for (i=0; i<neqns; i++)			{
	for (j=0; j<neqns; j++)
		printf(" %10.2e ",rjacobian(i,j));
	printf(" || %10.2e\n",rhs(i));
}
printf("---------------\n");
if (icount==3) return 0.0;*/
		// solve
		lusolve(rjacobian, rhs, rhs);
		
		// limit and update
		for (j=0; j<nspecies; j++)	{
			if (rhs(j) > 1.0) rhs(j) = 1.0;
			if (rhs(j) < -1.0) rhs(j) = -1.0;
			
			uval(j) = exp(log(uval(j)) + rhs(j));
		}
	
	}
	
	/* --- store values for return --- */
	residual = 0.0;
	for (j=0; j<nspecies; j++)	{
		if (uval(j) > 1.0e-15) residual += rhs(j)*rhs(j);
	}
	
	if (igastype != 1)	{
		for (j=igstart; j<igend; j++)
			concentration(j) = uval(j);
	}
	
	for (j=isstart; j<isend; j++)
		concentration(j) = uval(j);

	return residual;
}
/* ------------------------------------------------------------------------------- */



/* ------------------------ get_mean_thermal_speed ------------------------------- */
/* 
 * compute mean thermal speed of a gas species
 *
 */
double FiniteRateSurfaceModel::get_mean_thermal_speed(int i, double tmpr)
{	
	if ((i>=igstart)&&(i<igend)) 
		return sqrt( R * tmpr / (2.0 * PI * wmass(i)) );
	else
		return 0.0;
}
/* ------------------------------------------------------------------------------- */



/* -------------------------- get_species_name_pos ------------------------------- */
/* 
 * get position of species with a certain name for array ordering
 *
 */
int FiniteRateSurfaceModel::get_species_name_pos(const string &look_for, bool strip_braces, bool search_gas, bool search_surf, bool search_bulk)
{
	int i;
	size_t pos;
	
	// search gas species
	if (search_gas)	{
		for (i=igstart; i<igend; i++)	{
			if ( species(i).compare(look_for) == 0 ) return i;
		}
	}
	
	//search surface species
	if (search_surf)	{
		for (i=isstart; i<isend; i++)	{
			pos = std::string::npos;
			if (strip_braces) pos = species(i).find_first_of('(');
			if ( species(i).compare(0,pos,look_for) == 0 ) return i;
		}
	}
	
	//search bulk species
	if (search_bulk)	{
		for (i=ibstart; i<ibend; i++)	{
			pos = std::string::npos;
			if (strip_braces) pos = species(i).find_first_of('(');
			if ( species(i).compare(0,pos,look_for) == 0 ) return i;
		}
	}


	return -1;
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- readBlowingFile ---------------------------------- */
/* 
 * read in blowing file
 *
 */
int FiniteRateSurfaceModel::readBlowingFile(istream &stream)
{
	int status=0, n, k;
	double cval, sumkk;
	char buffer[INPUT_FILE_BUFFER_SIZE], *token=NULL;

	/* --- read header --- */
	// ignore two blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// read header
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off nsp token (integer)
	sscanf(token,"%d",&nebc);
	//
	if (nebc == 1)	{
		nblw = 1;
		ha   = 0.0;
	} else if (nebc == 2)	{
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off cyield token (fortran double)
		string_replace_char(token, 'd', 'e');
		string_replace_char(token, 'D', 'e');
		sscanf(token,"%lg",&cyield);
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off ha token (fortran double)
		string_replace_char(token, 'd', 'e');
		string_replace_char(token, 'D', 'e');
		sscanf(token,"%lg",&ha);
		nblw = 1;
		ha   = 0.0;
	} else {
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off nblw token (integer)
		sscanf(token,"%d",&nblw);
	}
	
	/* --- initialize memory --- */
	rmdotg.resize(nblw);
	blwmf.resize(nblw,ngps);
	// initialize
	rmdotg.initialize(0.0);
	blwmf.initialize(0.0);

	/* --- read body --- */
	if (nebc == 0)	{
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //read rmdotg values (fortran double)
		for (n=0; n<nblw; n++)	{
			string_replace_char(token, 'd', 'e');
			string_replace_char(token, 'D', 'e');
			sscanf(token,"%lg",&(rmdotg(n)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //read blwmf values (fortran double)
		}
	}
	for (k=0; k<ngps; k++)	{
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //read blwmf values (fortran double)
		for (n=0; n<nblw; n++)	{
			string_replace_char(token, 'd', 'e');
			string_replace_char(token, 'D', 'e');
			sscanf(token,"%lg",&(blwmf(n,k)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //read blwmf values (fortran double)
		}
	}
	
	/* --- compute effective mass addition terms --- */
	if (nebc > 0)	{
		cval = 0.0;
		if (cyield > 0.0) cval = 1.0/cyield - 1.0;
		sumkk = 0.0;
		for (k=0; k<ngps; k++)
			sumkk += blwmf(0,k)*wmass(k);
		if (sumkk == 0.0) sumkk = 1.0;
		for (k=0; k<ngps; k++)
			pyroadd(k) = blwmf(0,k)*cval/sumkk;

	} else {
	
		for (n=0; n<nblw; n++)	{
			sumkk = 0.0;
			for (k=0; k<ngps; k++)
				sumkk += blwmf(0,k)*wmass(k);
			if (sumkk == 0.0) sumkk = 1.0;
			for (k=0; k<ngps; k++)
				pyroadd(k) = pyroadd(k) + blwmf(n,k)*rmdotg(n)/sumkk;
		}
	
	}

	/* --- turn blowing on --- */
	isblwon = 1;

	return status;
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- readLewisGasFile --------------------------------- */
/* 
 * read Lewis gas data
 *
 */
int FiniteRateSurfaceModel::readLewisGasFile(istream &stream)
{
	int count=0, nrange, i, j, n;
	char buffer[INPUT_FILE_BUFFER_SIZE], name[INPUT_FILE_BUFFER_SIZE], *token=NULL;
	double h, s, cp;
	Array2D<double> cof(NLEWIS_FITS,NLEWIS_COEFFS);
	Array1D<double> tmprs(NLEWIS_FITS+1);

	//must have initialized the model with data
	if (initialized == false) return 1;
	
	//temperatures ranges in the gas file are fixed
	tmprs.initialize(0.0);
	tmprs(0) = 0.0;
	tmprs(1) = 200.0;
	tmprs(2) = 1000.0;
	tmprs(3) = 6000.0;
	tmprs(4) = 20000.0;
	tmprs(5) = 1.0e+20;

	// ignore two blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	
	//read Lewis data one species at a time and look for species name matches
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //species name
	strcpy(name, token);
	token = strtok(NULL, INPUT_STRTOK_TOKENS);   //# of ranges
	sscanf(token,"%d",&nrange);
	while	( (count<10000) && (stream.good()) && (string_is_end(name)) )	{
		
		// check range limits
		if (nrange > (NLEWIS_FITS-2)) return 2;
		//zero out coefficients
		cof.initialize(0.0);
		//read temperature lines
		for (i=0; i<nrange; i++)	{
			// line #1
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
			string_replace_char(buffer, 'd', 'e');
			string_replace_char(buffer, 'D', 'e');
			token = strtok(buffer, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,0)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,1)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,2)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,3)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,4)));
			// line #2
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
			string_replace_char(buffer, 'd', 'e');
			string_replace_char(buffer, 'D', 'e');
			token = strtok(buffer, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,5)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,6)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,7)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,8)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,9)));
		}
		
		//check to see if this matches any species in our model
		for (n=igstart; n<igend; n++)	{

			if (strcmp(name,species(n).c_str()) == 0)	{

				//set # of ranges
				rlewis_nrange(n) = nrange+2;

				//copy coefficients
				for (i=0; i<nrange; i++)	{
					for (j=0;j<NLEWIS_COEFFS; j++)
						rlewis_coeffs(j,i+1,n) = cof(i,j);
				}

				//extrapolate min range with constant Cp
				get_gibbs_data(rlewis_coeffs, tmprs(1), 1, n, &h, &s, &cp);
				rlewis_coeffs(2,0,n) = cp;
				rlewis_coeffs(8,0,n) = (h-cp)*tmprs(1);
				rlewis_coeffs(9,0,n) = s - cp*log(tmprs(1));
				//extrapolate max range with constant Cp
				get_gibbs_data(rlewis_coeffs, tmprs(nrange+1), nrange, n, &h, &s, &cp);
				rlewis_coeffs(2,nrange+1,n) = cp;
				rlewis_coeffs(8,nrange+1,n) = (h-cp)*tmprs(nrange+1);
				rlewis_coeffs(9,nrange+1,n) = s - cp*log(tmprs(nrange+1));

				// set tmpr ranges
				for (i=0; i<rlewis_nrange(n); i++)	{
					rlewis_trange(0,i,n) = tmprs(i);
					rlewis_trange(1,i,n) = tmprs(i+1);
				}

			}

		}
		
		//skip blank line
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		//read new name
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //species name
		strcpy(name, token);
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //# of ranges
		sscanf(token,"%d",&nrange);
		//increment count
		count++;
	}
      
      return 0;
}
/* ------------------------------------------------------------------------------- */



/* --------------------------- readLewisBulkFile --------------------------------- */
/* 
 * read Lewis bulk data
 *
 */
int FiniteRateSurfaceModel::readLewisBulkFile(istream &stream)
{
	int count=0, nrange, i, j, n;
	char buffer[INPUT_FILE_BUFFER_SIZE], name[INPUT_FILE_BUFFER_SIZE], 
	     modelname[INPUT_FILE_BUFFER_SIZE], *token=NULL;
	double h, s, cp;
	Array2D<double> cof(NLEWIS_FITS,NLEWIS_COEFFS);
	Array2D<double> tmprs(2,NLEWIS_FITS);

	//must have initialized the model with data
	if (initialized == false) return 1;

	//read Lewis data one species at a time and look for species name matches
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //species name
	strcpy(name, token);
	token = strtok(NULL, INPUT_STRTOK_TOKENS);   //# of ranges
	sscanf(token,"%d",&nrange);

	while	( (count<10000) && (stream.good()) && (string_is_end(name)) )	{
		// check range limits
		if (nrange > (NLEWIS_FITS-2)) return 2;
		//zero out coefficients
		cof.initialize(0.0);
		tmprs.initialize(0.0);
		//read temperature lines
		for (i=0; i<nrange; i++)	{
			//temperature ranges
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
			string_replace_char(buffer, 'd', 'e');
			string_replace_char(buffer, 'D', 'e');
			token = strtok(buffer, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&( tmprs(0,i) ));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&( tmprs(1,i) ));
			// line #1
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
			string_replace_char(buffer, 'd', 'e');
			string_replace_char(buffer, 'D', 'e');
			token = strtok(buffer, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,0)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,1)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,2)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,3)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,4)));
			// line #2
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
			string_replace_char(buffer, 'd', 'e');
			string_replace_char(buffer, 'D', 'e');
			token = strtok(buffer, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,5)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,6)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,7)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,8)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS);
			sscanf(token,"%lg",&(cof(i,9)));
		}

		//check to see if this matches any species in our model
		for (n=ibstart; n<nspecies; n++)	{

			//in the model, species names have something like (b1) at end, so we strip that off
			strcpy(modelname, species(n).c_str());
			token = strchr(modelname,'(');
			if (token) token[0] = '\0';  //just null this character out so the name ends here
			if (strcmp(name,modelname) == 0)	{

				//set # of ranges
				rlewis_nrange(n) = nrange+2;

				//copy coefficients
				for (i=0; i<nrange; i++)	{
					for (j=0;j<NLEWIS_COEFFS; j++)
						rlewis_coeffs(j,i+1,n) = cof(i,j);
					rlewis_trange(0,i+1,n)  = tmprs(0,i);
					rlewis_trange(1,i+1,n)  = tmprs(1,i);
				}

				//extrapolate min range with constant Cp
				get_gibbs_data(rlewis_coeffs, tmprs(0,0), 1, n, &h, &s, &cp);
				rlewis_coeffs(2,0,n) = cp;
				rlewis_coeffs(8,0,n) = (h-cp)*tmprs(0,0);
				rlewis_coeffs(9,0,n) = s - cp*log(tmprs(0,0));
				rlewis_trange(0,0,n) = 0.0;
				rlewis_trange(1,0,n) = tmprs(0,0);
				//extrapolate max range with constant Cp
				get_gibbs_data(rlewis_coeffs, tmprs(1,nrange-1), nrange, n, &h, &s, &cp);
				rlewis_coeffs(2,nrange+1,n) = cp;
				rlewis_coeffs(8,nrange+1,n) = (h-cp)*tmprs(1,nrange-1);
				rlewis_coeffs(9,nrange+1,n) = s - cp*log(tmprs(1,nrange-1));
				rlewis_trange(0,nrange+1,n) = tmprs(1,nrange-1);
				rlewis_trange(1,nrange+1,n) = 1.0e+20;
			}

		}
		
		//skip blank line
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		//read new name
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //species name
		strcpy(name, token);
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //# of ranges
		sscanf(token,"%d",&nrange);
		//increment count
		count++;
	}
      
      return 0;
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- sanityCheck -------------------------------------- */
/* 
 * perform sanity checks on the model data
 *
 */
int FiniteRateSurfaceModel::sanityCheck()
{
	int i,j,n,m,k,iph,ier,irsum[3],ipsum[3],
	    iphar[ISMODEL_MAX_REACT], iphap[ISMODEL_MAX_PROD], sanity_check=0;
	double rpsum, rrsum;

	//---- WARNING #1: Make sure we have desorption rates for all surface
	//                 species if calculating thermo data is required
	//
	//     This is a WARNING, but the following occurs:
	//             - ALL backward rates for reactions involving this  
	//               species are implicitly turned off
	//             - Gibbs energy is set to zero explicitly (kthermo = 2)
	k = ngps-1;
	for (n=0; n<nsp; n++)	{
		for (j=0; j<nsps(n); j++)	{
			k++;
			if ( (kthermo(k)==0) && (kadsr(k)<0) )	{
				if(sanity_check <= 1) sanity_check = 1;
				snprintf(cfrm_error_string, 1024, "missing desorption data for species %5d: some backward rates set to zero",k+1);
				kthermo(k) = 2;
				for (i=0; i<isrt; i++)	{
					for (m=0; m<ISMODEL_MAX_REACT; m++)
						if (RR(i,m) == k) isrbon(i) = 0;
					for (m=0; m<ISMODEL_MAX_PROD; m++)
						if (RP(i,m) == k) isrbon(i) = 0;
				}
			}
		}
	}

	//---- WARNING #2: Pyrolysis cannot be active with no bulk phase
	if ((isblwon > 0) && (nbp < 1))	{
		if(sanity_check <= 2) sanity_check = 2;
		snprintf(cfrm_error_string, 1024, "pyrolysis blowing cannot be active without at least one bulk phase");
		isblwon  = 0;
		nblwflag = 0;
	}

	//---- WARNING #3: reactions must have proper participant types
	for (i=0; i<isrt; i++)	{
		
		// count species in each phase
		irsum[0] = 0;	irsum[1] = 0;	irsum[2] = 0;
		ipsum[0] = 0;	ipsum[1] = 0;	ipsum[2] = 0;
		for (m=0; m<ISMODEL_MAX_REACT; m++)	{
			iph = -1;
			if (RR(i,m) >= igstart) iph=0;
			if (RR(i,m) >= igend)   iph=1;
			if (RR(i,m) >= isend)   iph=2;
			if (iph>=0) irsum[iph] += RVR(i,m);
		}
		for (m=0; m<ISMODEL_MAX_PROD; m++)	{
			iph = -1;
			if (RP(i,m) >= igstart) iph=0;
			if (RP(i,m) >= igend)   iph=1;
			if (RP(i,m) >= isend)   iph=2;
			if (iph>=0) ipsum[iph] += RVP(i,m);
		}
		
		//issue warnings
		if (istype(i) == 2)	{   //E-R
			if ((irsum[0] < 1) || (irsum[1] < 1) || (ipsum[0] < 1) || (ipsum[1] < 1))	{
				if(sanity_check <= 3) sanity_check = 3;
				snprintf(cfrm_error_string, 1024, "E-R reaction not of form A+[B]<-->AB+[] for reaction #%4d",i+1);
			}
		}
		else if (istype(i) == 3)	{   //L-H
			if ((irsum[0] < 0) || (irsum[1] < 2) || (ipsum[0] < 1) || (ipsum[1] < 2))	{
				if(sanity_check <= 3) sanity_check = 3;
				snprintf(cfrm_error_string, 1024, "L-H reaction not of form [A]+[B]<-->AB+2[] for reaction #%4d",i+1);
			}
		}
		else if (istype(i) == 4)	{   //sublimation
			if ((irsum[0] < 0) || (irsum[1] < 1) || (ipsum[0] < 1) || (ipsum[1] < 1) || ((irsum[2]+ipsum[2]) < 1))	{
				if(sanity_check <= 3) sanity_check = 3;
				snprintf(cfrm_error_string, 1024, "sublimation reaction not of form [A]+<B><-->AB+[] for reaction #%4d",i+1);
			}
		}
		else { // adsorption (type 1 or 5) or Arrhenius (type 0)
			if ((irsum[0] < 1) || (irsum[1] < 1) || (ipsum[0] < 0) || (ipsum[1] < 1))	{
				if(sanity_check <= 3) sanity_check = 3;
				snprintf(cfrm_error_string, 1024, "Adsorption reaction not of form A+[]<-->[A] for reaction #%4d",i+1);
			}
		}
	
	}

	//---- WARNING #4: mass balance for reactions
	for (i=0; i<isrt; i++)	{
		
		//sum mass for each reaction
		rrsum = 0.0;
		rpsum = 0.0;
		for (m=0; m<ISMODEL_MAX_REACT; m++)
			if (RR(i,m) >= 0)	rrsum += (double)(RVR(i,m)) * wmass(RR(i,m));
		for (m=0; m<ISMODEL_MAX_PROD; m++)
			if (RP(i,m) >= 0)	rpsum += (double)(RVP(i,m)) * wmass(RP(i,m));
		if (rrsum < 1.0e-30) rrsum = 1.0e-30;
		if (rpsum < 1.0e-30) rpsum = 1.0e-30;
		if ( abs( (rrsum-rpsum)/(rrsum+rpsum) ) >= 0.01 )	{
			if(sanity_check <= 4) sanity_check = 4;
			snprintf(cfrm_error_string, 1024, "mass does not balance for reaction #%4d",i+1);
		}
	}

	//---- WARNING #5: reactions can have only one gas or bulk
	for (i=0; i<isrt; i++)	{
		
		// count species in each phase
		irsum[0] = 0;	irsum[1] = 0;	irsum[2] = 0;
		ipsum[0] = 0;	ipsum[1] = 0;	ipsum[2] = 0;
		for (m=0; m<ISMODEL_MAX_REACT; m++)	{
			iph = -1;
			if (RR(i,m) >= igstart) iph=0;
			if (RR(i,m) >= igend)   iph=1;
			if (RR(i,m) >= isend)   iph=2;
			if (iph>=0) irsum[iph] += 1;
		}
		for (m=0; m<ISMODEL_MAX_PROD; m++)	{
			iph = -1;
			if (RP(i,m) >= igstart) iph=0;
			if (RP(i,m) >= igend)   iph=1;
			if (RP(i,m) >= isend)   iph=2;
			if (iph>=0) ipsum[iph] += 1;
		}
		if ( (irsum[0] > 1) || (irsum[2] > 1) || (ipsum[0] > 1) || (ipsum[2] > 1) )	{
			if(sanity_check <= 5) sanity_check = 5;
			snprintf(cfrm_error_string, 1024, "reaction #%4d is not an elementary reaction",i+1);
		}
	}
	
	//---- ERROR #101: Make sure we do not have an empty surface site
	//                 as the first reactant or product
	for (i=0; i<isrt; i++)	{
		// identify gas phase reactant participant
		k = RR(i,0);
		if ((istype(i)>0)&&(istype(i)<4)&&(wmass(k)==0.0))	{
			sanity_check = 101;
			snprintf(cfrm_error_string, 1024, "cannot set empty surface site as first reactant for types 1,2,3");
			return sanity_check;
		}
		// identify gas phase product participant
		k = RP(i,0);
		if ((istype(i)==4)&&(wmass(k)==0.0))	{
			sanity_check = 101;
			snprintf(cfrm_error_string, 1024, "cannot set empty surface site as first product for sublimation");
			return sanity_check;
		}
	}

	//---- ERROR #102: Insufficient array size for Lewis data
	for (n=0; n<nspecies; n++)	{
		if (rlewis_nrange(n) > NLEWIS_FITS)	{
			sanity_check = 102;
			snprintf(cfrm_error_string, 1024, "too many Lewis ranges; increase the size of the NLEWIS_FITS parameter");
			return sanity_check;
		}
	}
	
	//---- ERROR #103: Make sure participants are in order: gas/surface/bulk
	for (i=0; i<isrt; i++)	{
		// get phase of each participant
		for (m=0; m<ISMODEL_MAX_REACT; m++)	{
			iph = -1;
			if (RR(i,m) == -1)        iph=3;  //zero entry
			if (RR(i,m) >= igstart)   iph=0;  //gas phase
			if (RR(i,m) >= isstart)   iph=1;  //surface phase
			if ( (RR(i,m) >= ibstart) 
			  && (RR(i,m) < ibend) )  iph=2;  //bulk phase
			iphar[m] = iph;
		}
		for (m=0; m<ISMODEL_MAX_PROD; m++)	{
			iph = -1;
			if (RP(i,m) == -1)        iph=3;  //zero entry
			if (RP(i,m) >= igstart)   iph=0;  //gas phase
			if (RP(i,m) >= isstart)   iph=1;  //surface phase
			if ( (RP(i,m) >= ibstart) 
			  && (RP(i,m) < ibend) )  iph=2;  //bulk phase
			iphap[m] = iph;
		}
		for (m=1; m<ISMODEL_MAX_REACT; m++)	{
			if (iphar[m] < iphar[m-1])	{
				sanity_check = 103;
				snprintf(cfrm_error_string, 1024, "reaction #%4d participants must be ordered: gas, surface, bulk, then zero entries",i+1);
				return sanity_check;
			}
		}
		for (m=1; m<ISMODEL_MAX_PROD; m++)	{
			if (iphap[m] < iphap[m-1])	{
				sanity_check = 103;
				snprintf(cfrm_error_string, 1024, "reaction #%4d products must be ordered: gas, surface, bulk, then zero entries",i+1);
				return sanity_check;
			}
		}
	}
	
	//---- ERROR #104: Look for negative stoichiometric coefficients
	for (i=0; i<isrt; i++)	{
		for (m=0; m<ISMODEL_MAX_REACT; m++)	
			if (RVR(i,m) < 0) sanity_check = 104;
		for (m=0; m<ISMODEL_MAX_PROD; m++)
			if (RVP(i,m) < 0) sanity_check = 104;
		if (sanity_check > 100)	{
			snprintf(cfrm_error_string, 1024, "reaction #%4d cannot have negative stoichiometric coefficients",i+1);
			return sanity_check;
		}
	}
	
	//---- ERROR #105: Look for bad participants
	for (i=0; i<isrt; i++)	{
		for (m=0; m<ISMODEL_MAX_REACT; m++)	
			if ((RR(i,m) < -1)||(RR(i,m) >= nspecies)) sanity_check = 105;
		for (m=0; m<ISMODEL_MAX_PROD; m++)
			if ((RP(i,m) < -1)||(RP(i,m) >= nspecies)) sanity_check = 105;
		if (sanity_check > 100)	{
			snprintf(cfrm_error_string, 1024, "reaction #%4d lists a bad participant",i+1);
			return sanity_check;
		}
	}
	
	//---- ERROR #106: bad reaction types
	for (i=0; i<isrt; i++)	{
		if ((istype(i) < 0)||(istype(i) > 5))	{
			sanity_check = 106;
			snprintf(cfrm_error_string, 1024, "reaction #%4d is not a valid type",i+1);
			return sanity_check;
		}
	}
	
	//---- ERROR #107: bad molecular weight
	for (k=0; k<nspecies; k++)	{
		if (wmass(k) < 0.0)	{
			sanity_check = 107;
			snprintf(cfrm_error_string, 1024, "species #%4d cannot have negative molecular weight",k+1);
			return sanity_check;
		}
	}

	return sanity_check;
}
/* ------------------------------------------------------------------------------- */



/* ------------------------------ operator>> ------------------------------------- */
/* 
 * input operator
 *
 */
istream &operator>>(istream &stream, FiniteRateSurfaceModel &frm)
{
	char buffer[INPUT_FILE_BUFFER_SIZE], *token=NULL;
	int i, j, k, nphase, nsurf, nsa, nbulk, etable_error=0;
	double Edsum;


	/* --- Header --- */
	//read title line
	get_line_of_text(stream, frm.title, 200);
	//ignore 3 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off nsp token (integer)
	sscanf(token,"%d",&frm.nsp);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off nbp token (integer)
	sscanf(token,"%d",&frm.nbp);
	//ignore 3 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off nblwflag token (integer)
	sscanf(token,"%d",&frm.nblwflag);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off initsurf token (integer)
	sscanf(token,"%d",&frm.initsurf);
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// compute number of phases
	frm.nphases = 1 + frm.nsp + frm.nbp;


	/* --- gas phase declaration --- */
	// allocate memory
	frm.phases.resize(frm.nphases);
	// read gas declaration line
	nphase = 0;
	// skip 2 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); // split off name token (char *)
	frm.phases(nphase).assign(token);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off ngps token (integer)
	sscanf(token,"%d",&frm.ngps);
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	frm.nspecies = frm.ngps;


	/* --- surface phase declarations --- */
	// allocate memory
	frm.nspas.resize(frm.nsp);
	frm.sfrc.resize(frm.nsp);
	frm.ithermo.resize(frm.nsp);
	// skip 3 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	frm.nspasmax = 0;
	frm.nssites  = 0;
	for (nsurf=0; nsurf<frm.nsp; nsurf++)	{
		nphase++;
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);

		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off name token (char *)
		frm.phases(nphase).assign(token);
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off sfrc token (fortran double)
		string_replace_char(token, 'd', 'e');
		sscanf(token,"%lg",&(frm.sfrc(nsurf)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off nspas (integer)
		sscanf(token,"%d",&(frm.nspas(nsurf)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off ithermo (integer)
		sscanf(token,"%d",&(frm.ithermo(nsurf)));
		// find maximum size of any phase
		if ( frm.nspas(nsurf) > frm.nspasmax ) frm.nspasmax = frm.nspas(nsurf);
		// count number of site sets
		frm.nssites += frm.nspas(nsurf);

	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);


	/* --- surface active site declarations --- */
	// allocate
	frm.nspass.resize(frm.nsp,frm.nspasmax);
	frm.sdenas.resize(frm.nsp,frm.nspasmax);
	frm.nsps.resize(frm.nsp);
	frm.sden.resize(frm.nsp);
	// initialize
	frm.nspass.initialize(0);
	frm.sdenas.initialize(0.0);
	frm.sden.initialize(0.0);
	frm.nsps.initialize(0);
	frm.sdentot = 0.0;
	frm.sfrcI   = 0.0;
	//  skip 3 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// read surface phase active site declarations
	for (nsurf=0; nsurf<frm.nsp; nsurf++)	{
		for (nsa=0; nsa<frm.nspas(nsurf); nsa++)	{
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
			token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off sdenas token (fortran double)
			string_replace_char(token, 'd', 'e');
			sscanf(token,"%lg",&(frm.sdenas(nsurf,nsa)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off nspass (integer)
			sscanf(token,"%d",&(frm.nspass(nsurf,nsa)));
			frm.nsps(nsurf) = frm.nsps(nsurf) + frm.nspass(nsurf,nsa);
			frm.sden(nsurf) = frm.sden(nsurf) + frm.sdenas(nsurf,nsa);
		}
		frm.nspecies = frm.nspecies + frm.nsps(nsurf);
		frm.sdentot  = frm.sdentot  + frm.sfrc(nsurf) * frm.sden(nsurf);
		if (frm.nspas(nsurf) == 0) frm.sfrcI = frm.sfrcI + frm.sfrc(nsurf);
	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);


	/* --- bulk phase declarations --- */
	// allocate memory
	frm.den.resize(frm.nbp);
	frm.nbps.resize(frm.nbp);
	frm.por.resize(frm.nbp);
	frm.vfb.resize(frm.nbp);
	// initialize
	frm.vfbI    = 0.0;
	frm.bdentot = 0.0;
	frm.bportot = 0.0;
	if (frm.nbp > 0)	{
		// skip 3 blank lines
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		// read bulk delcarations
		for (nbulk=0; nbulk<frm.nbp; nbulk++)	{
			nphase++;
			get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);

			token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off name token (char *)
			frm.phases(nphase).assign(token);
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off den token (fortran double)
			string_replace_char(token, 'd', 'e');
			sscanf(token,"%lg",&(frm.den(nbulk)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off por token (fortran double)
			string_replace_char(token, 'd', 'e');
			sscanf(token,"%lg",&(frm.por(nbulk)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off vfb token (fortran double)
			string_replace_char(token, 'd', 'e');
			sscanf(token,"%lg",&(frm.vfb(nbulk)));
			token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off nbps (integer)
			sscanf(token,"%d",&(frm.nbps(nbulk)));

			frm.nspecies = frm.nspecies + frm.nbps(nbulk);
			if (frm.nbps(nbulk) == 0) frm.vfbI = frm.vfbI + frm.vfb(nbulk);
			frm.bdentot = frm.bdentot + frm.vfb(nbulk) * frm.den(nbulk);
			frm.bportot = frm.bportot + frm.vfb(nbulk) * frm.por(nbulk);
		}
		// skip 1 blank line
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	}


	/* --- compute start/end indices --- */
	// note that C++ indexes from 0 whereas Fortran indexes from 1
	frm.igstart = 0;
	frm.igend   = frm.ngps;
	frm.isstart = frm.igend;
	frm.isend = frm.isstart;
	for (nsurf=0; nsurf<frm.nsp; nsurf++) frm.isend = frm.isend + frm.nsps(nsurf);
	frm.ibstart = frm.isend;
	frm.ibend   = frm.nspecies;


	/* --- read species lists --- */
	// allocate memory
	frm.species.resize(frm.nspecies);
	frm.wmass.resize(frm.nspecies);
	frm.Ediss.resize(frm.nspecies);
	frm.Ed.resize(frm.nspecies);
	frm.kphase.resize(frm.nspecies);
	frm.kthermo.resize(frm.nspecies);
	frm.ksites.resize(frm.nspecies);
	frm.bmf.resize(frm.nbp,frm.nspecies);
	// initialize
	frm.Ediss.initialize(0.0);
	frm.Ed.initialize(0.0);
	frm.kphase.initialize(0);
	frm.kthermo.initialize(1);
	frm.ksites.initialize(0);
	frm.bmf.initialize(0.0);
	// skip 2 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// read gas phase species list
	for (i=frm.igstart; i<frm.igend; i++)	{
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);

		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off name token (char *)
		frm.species(i).assign(token);
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off wmass token (fortran double)
		string_replace_char(token, 'd', 'e');
		sscanf(token,"%lg",&(frm.wmass(i)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off Ediss token (fortran double)
		string_replace_char(token, 'd', 'e');
		sscanf(token,"%lg",&(frm.Ediss(i)));

		frm.kphase(i) = 0;
	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//read surface phase species list
	i = frm.ngps-1;
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	for (nsurf=0; nsurf<frm.nsp; nsurf++)	{
		for (nsa=0; nsa<frm.nspas(nsurf); nsa++)	{
			for (j=0; j<frm.nspass(nsurf,nsa); j++)	{
				i++;
				get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);

				token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off name token (char *)
				frm.species(i).assign(token);
				token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off wmass token (fortran double)
				string_replace_char(token, 'd', 'e');
				sscanf(token,"%lg",&(frm.wmass(i)));
				token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off Ed token (fortran double)
				string_replace_char(token, 'd', 'e');
				sscanf(token,"%lg",&(frm.Ed(i)));

				frm.kphase(i)  = 1+nsurf;
				frm.ksites(i)  = nsa;
				frm.kthermo(i) = frm.ithermo(nsurf);

			}
		}
	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	if (frm.nbp > 0)	{
		// skip 2 blank lines
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		// read bulk phase species list
		for (nbulk=0; nbulk<frm.nbp; nbulk++)	{
			for (j=0; j<frm.nbps(nbulk); j++)	{
				i++;
				get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);

				token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off name token (char *)
				frm.species(i).assign(token);
				token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off wmass token (fortran double)
				string_replace_char(token, 'd', 'e');
				sscanf(token,"%lg",&(frm.wmass(i)));
				token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off bmf token (fortran double)
				string_replace_char(token, 'd', 'e');
				sscanf(token,"%lg",&(frm.bmf(nbulk,j)));

				frm.kphase(i)  = 1+frm.nsp+nbulk;
			}
		}
		// skip 1 blank line
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	}


	/* --- read reaction information --- */
	// skip 2 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// read header
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off isrt token (integer)
	sscanf(token,"%d",&frm.isrt);

	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// allocate memory
	frm.RR.resize(frm.isrt,ISMODEL_MAX_REACT);
	frm.RP.resize(frm.isrt,ISMODEL_MAX_PROD);
	frm.RVR.resize(frm.isrt,ISMODEL_MAX_REACT);
	frm.RVP.resize(frm.isrt,ISMODEL_MAX_PROD);
	frm.vr.resize(frm.isrt,frm.nspecies);
	frm.vp.resize(frm.isrt,frm.nspecies);
	frm.v.resize(frm.isrt,frm.nspecies);
	// initialize
	frm.RR.initialize(0);
	frm.RP.initialize(0);
	frm.RVR.initialize(0);
	frm.RVP.initialize(0);
	frm.vr.initialize(0.0);
	frm.vp.initialize(0.0);
	frm.v.initialize(0.0);
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// read reactant/product species
	for (i=0; i<frm.isrt; i++)	{
		// read from stream
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off RR1 (integer)
		sscanf(token,"%d",&(frm.RR(i,0)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RR2 (integer)
		sscanf(token,"%d",&(frm.RR(i,1)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RR3 (integer)
		sscanf(token,"%d",&(frm.RR(i,2)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RP1 (integer)
		sscanf(token,"%d",&(frm.RP(i,0)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RP2 (integer)
		sscanf(token,"%d",&(frm.RP(i,1)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RP3 (integer)
		sscanf(token,"%d",&(frm.RP(i,2)));
		// here, we must be careful.  Species1 is called "1", but C++ indexing starts with 0
		frm.RR(i,0)--;
		if ( (frm.RR(i,0)<0) || (frm.RR(i,0)>=frm.nspecies) ) frm.RR(i,0) = -1;
		frm.RR(i,1)--;
		if ( (frm.RR(i,1)<0) || (frm.RR(i,1)>=frm.nspecies) ) frm.RR(i,1) = -1;
		frm.RR(i,2)--;
		if ( (frm.RR(i,2)<0) || (frm.RR(i,2)>=frm.nspecies) ) frm.RR(i,2) = -1;
		frm.RP(i,0)--;
		if ( (frm.RP(i,0)<0) || (frm.RP(i,0)>=frm.nspecies) ) frm.RP(i,0) = -1;
		frm.RP(i,1)--;
		if ( (frm.RP(i,1)<0) || (frm.RP(i,1)>=frm.nspecies) ) frm.RP(i,1) = -1;
		frm.RP(i,2)--;
		if ( (frm.RP(i,2)<0) || (frm.RP(i,2)>=frm.nspecies) ) frm.RP(i,2) = -1;
	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// read reactant/product stoichiometric coefficients
	for (i=0; i<frm.isrt; i++)	{
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off RVR1 (integer)
		sscanf(token,"%d",&(frm.RVR(i,0)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RVR2 (integer)
		sscanf(token,"%d",&(frm.RVR(i,1)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RVR3 (integer)
		sscanf(token,"%d",&(frm.RVR(i,2)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RVP1 (integer)
		sscanf(token,"%d",&(frm.RVP(i,0)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RVP2 (integer)
		sscanf(token,"%d",&(frm.RVP(i,1)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS); //split off RVP3 (integer)
		sscanf(token,"%d",&(frm.RVP(i,2)));
		for (j=0; j<ISMODEL_MAX_REACT; j++)	{
			if (frm.RR(i,j) >= 0)
				frm.vr(i,frm.RR(i,j)) += (double)(frm.RVR(i,j));
		}
		for (j=0; j<ISMODEL_MAX_PROD; j++)	{
			if (frm.RP(i,j) >= 0)
				frm.vp(i,frm.RP(i,j)) += (double)(frm.RVP(i,j));
		}
		for (j=0; j<frm.nspecies; j++)
			frm.v(i,j) = frm.vp(i,j) - frm.vr(i,j);
	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);

	
	/* --- read rate constant information --- */
	// skip 8 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// allocate memory
	frm.istype.resize(frm.isrt);
	frm.isrfon.resize(frm.isrt);
	frm.isrbon.resize(frm.isrt);
	frm.Cf.resize(frm.isrt);
	frm.eta.resize(frm.isrt);
	frm.Ea.resize(frm.isrt);
	// initialize
	frm.istype.initialize(0);
	frm.isrfon.initialize(1);
	frm.isrbon.initialize(1);
	frm.Cf.initialize(0.0);
	frm.eta.initialize(0.0);
	frm.Ea.initialize(0.0);
	//
	for (i=0; i<frm.isrt; i++)	{
		get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
		
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off istype token (integer)
		sscanf(token,"%d",&(frm.istype(i)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off Cf token (fortran double)
		string_replace_char(token, 'd', 'e');
		sscanf(token,"%lg",&(frm.Cf(i)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off eta token (fortran double)
		string_replace_char(token, 'd', 'e');
		sscanf(token,"%lg",&(frm.eta(i)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off Ea token (fortran double)
		string_replace_char(token, 'd', 'e');
		sscanf(token,"%lg",&(frm.Ea(i)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off isrfon token (integer)
		sscanf(token,"%d",&(frm.isrfon(i)));
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off isrbon token (integer)
		sscanf(token,"%d",&(frm.isrbon(i)));
	}
	// skip 1 blank line
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);


	/* --- read rate constant information --- */
	// skip 13 blank lines
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE); 
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE); 
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
	// allocate memory
	frm.iadtype.resize(frm.isrt);
	frm.iadform.resize(frm.isrt);
	frm.Cfad.resize(frm.isrt);
	frm.etaad.resize(frm.isrt);
	frm.vad.resize(frm.isrt);
	frm.Edes.resize(frm.isrt);
	frm.kadsr.resize(frm.nspecies);
	// initialize
	frm.iadtype.initialize(0);
	frm.iadform.initialize(0);
	frm.kadsr.initialize(-1);
	frm.Cfad.initialize(0.0);
	frm.etaad.initialize(0.0);
	frm.Edes.initialize(0.0);
	frm.vad.initialize(0.0);
	//
	for (i=0; i<frm.isrt; i++)	{
		if ( (frm.istype(i) == 1) || (frm.istype(i) == 5) )	{
			//identify species
			k = frm.RP(i,0);
			if (k >= 0) {
				frm.kadsr(k) = i;
				//if no thermodynamic data for this species, read
				if (frm.kthermo(k) == 0)	{
					get_line_of_text(stream, buffer, INPUT_FILE_BUFFER_SIZE);
				
					token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off iadtype token (integer)
					sscanf(token,"%d",&(frm.iadtype(i)));
					token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off iadform token (integer)
					sscanf(token,"%d",&(frm.iadform(i)));
					token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off Cfad token (fortran double)
					string_replace_char(token, 'd', 'e');
					sscanf(token,"%lg",&(frm.Cfad(i)));
					token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off etaad token (fortran double)
					string_replace_char(token, 'd', 'e');
					sscanf(token,"%lg",&(frm.etaad(i)));
					token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off vad token (fortran double)
					string_replace_char(token, 'd', 'e');
					sscanf(token,"%lg",&(frm.vad(i)));
					token = strtok(NULL, INPUT_STRTOK_TOKENS);   //split off Edes token (fortran double)
					string_replace_char(token, 'd', 'e');
					sscanf(token,"%lg",&(frm.Edes(i)));

					frm.Ed(k) = frm.Edes(i);
				}
			}
			
		}
	}

	/* --- compute derived quantities for model --- */
	// allocate memory
	frm.Eact.resize(frm.isrt);
	frm.sitep.resize(frm.isrt);
	frm.sdensity.resize(frm.isrt);
	frm.sfraction.resize(frm.isrt);
	frm.pyroadd.resize(frm.ngps);
	frm.sumv.resize(frm.isrt,frm.nphases);
	// initialize
	frm.Eact.initialize(0.0);
	frm.sitep.initialize(0.0);
	frm.sdensity.initialize(0.0);
	frm.sfraction.initialize(0.0);
	frm.pyroadd.initialize(0.0);
	frm.sumv.initialize(0.0);
	// compute activation energy and total site density for all reactions
	for (i=0; i<frm.isrt; i++)	{
		frm.Eact(i) = frm.Ea(i);
		if (frm.istype(i) == 3) {
			Edsum = 0.0;
			for (j=0; j<ISMODEL_MAX_REACT; j++)	{
				k = frm.RR(i,j);
				if (k >= 0)
					Edsum = Edsum + frm.RVR(i,j) * frm.Ed(k);
			}
			k = frm.RP(i,0);
			if (k >= 0) Edsum = Edsum - frm.Ediss(k);
			if (Edsum > frm.Eact(i)) frm.Eact(i) = Edsum;
		}
		
		for (j=0; j<ISMODEL_MAX_REACT; j++)	{
			k = frm.RR(i,j);
			if (k >= 0)	{
				nphase = frm.kphase(k);
				if ( (nphase>0) && (nphase <= frm.nsp) )	{  //index shift from FORTRAN
					frm.sitep(i) += (double)(frm.RVR(i,j));
					frm.sdensity(i)  = frm.sden(nphase-1);
					frm.sfraction(i) = frm.sfrc(nphase-1);
				}
			}
		}
		
		for (j=frm.igstart; j<frm.isend; j++)	{
			nphase = frm.kphase(j);
			frm.sumv(i,nphase) += frm.v(i,j);
		}
	}

	/* --- set sizing of Lewis matrices --- */
	frm.rlewis_coeffs.resize(NLEWIS_COEFFS,NLEWIS_FITS,frm.nspecies);
	frm.rlewis_coeffs.initialize(0.0);
	frm.rlewis_trange.resize(2,NLEWIS_FITS,frm.nspecies);
	frm.rlewis_trange.initialize(0.0);
	frm.rlewis_nrange.resize(frm.nspecies);
	frm.rlewis_nrange.initialize(0);

	/* --- most species don't use Gibbs data but empty sites MUST use Gibbs data so it zeros out --- */
	i = frm.ngps-1;
	for (nsurf=0; nsurf<frm.nsp; nsurf++)	{
		for (nsa=0; nsa<frm.nspas(nsurf); nsa++)	{
			for (j=0; j<frm.nspass(nsurf,nsa); j++)	{
				i++;
				if (j == 0)	{
					frm.kthermo(i) = 1;
					frm.rlewis_trange(0,0,i) = 0.0;
					frm.rlewis_trange(1,0,i) = 1.0e+20;
					frm.rlewis_nrange(i)     = 1;
				}
			}
		}
	}

	/* --- set miscellaneous fields --- */
	//turn off debugging
	frm.ifrm_debug = 0;
	//build element table
	etable_error = frm.build_etable();
	//turn off blowing for now (file gets read separately)
	frm.isblwon = 0;
	frm.nblw    = 0;
	frm.cyield  = 1.0;
	frm.ha      = 0.0;
	//initialized
	frm.initialized = true;

	return stream; //must return the stream 
}
/* ------------------------------------------------------------------------------- */




/* ------------------------------ operator<< ------------------------------------- */
/* 
 * output operator
 *
 */
ostream &operator<<(ostream &stream, FiniteRateSurfaceModel &frm)
{
	char buffer[INPUT_FILE_BUFFER_SIZE], unixnewline = 0x0A;  //LineFeed
	int i, n, k, m, j;
	double Edes;

	/* --- header --- */
	stream << frm.title << unixnewline;
	stream << "------------------------------------------------------------------------" << unixnewline;
	stream << "Number of Surface and bulk phases (gas phase=1 by definition)" << unixnewline;
	stream << "nsp, nbp" << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-4d %-4d",frm.nsp, frm.nbp);
	stream << buffer << unixnewline;
	stream << unixnewline;

	stream << "Blowing/pyrolyzing gas flows? (0=NO, 1=Yes) and Surface Initialization (0=empty sites, 1=QSS)" << unixnewline;
	stream << "nblwflag     initsurf" << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-1d            %-1d",frm.nblwflag, frm.initsurf);
	stream << buffer << unixnewline;
	stream << unixnewline;

	/* --- gas phase declaration --- */
	stream << "Number of gas phase species participating in surface reactions" << unixnewline;
	stream << "Name        ngps   phase#" << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-10s  %-3d      %1d",frm.phases(0).c_str(), frm.ngps, 1);
	stream << buffer << unixnewline;
	stream << unixnewline;

	/* --- surface phase declaration --- */
	stream << "For each surface phase: list name, surface fraction and" << unixnewline;
	stream << "number of active site sets, and thermo availability (0=No,1=Yes)" << unixnewline;
	stream << "Name       sfrc     nspas  iThermo  phase#" << unixnewline;
	for (n=0; n<frm.nsp; n++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-10s %-5.3lf    %-3d    %-1d        %3d",frm.phases(n+1).c_str(), frm.sfrc(n), frm.nspas(n),frm.ithermo(n),n+2);
		stream << buffer << unixnewline;
	}
	stream << unixnewline;

	/* --- surface active site declaration --- */
	stream << "For each surface phase with 1 or more sets of active sites, list" << unixnewline;
	stream << "the site density and number of species for each active site set" << unixnewline;
	stream << "sdenas (mol/m2)      nspass  phase#/site#" << unixnewline;
	for (n=0; n<frm.nsp; n++)	{
		for (k=0; k<frm.nspas(n); k++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-12.5le         %-3d        %3d/%-3d",frm.sdenas(n,k),frm.nspass(n,k),n+2,k+1);
			stream << buffer << unixnewline;
		}
	}
	stream << unixnewline;

	/* --- bulk phase declaration --- */
	if (frm.nbp > 0)	{
		stream << "For each bulk phase: list name, density, porosity, " << unixnewline;
		stream << "volume fraction and number of bulk species " << unixnewline;
		stream << "Name       sden         porosity  vol. fract.  nbps   phase#" << unixnewline;
		for (n=0; n<frm.nbp; n++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-10s %-10.4lf   %-7.4lf   %-7.4lf      %-3d    %3d",frm.phases(frm.nsp+1+n).c_str(), frm.den(n), frm.por(n),frm.vfb(n),frm.nbps(n),frm.nsp+2+n);
			stream << buffer << unixnewline;
		}
		stream << unixnewline;
	}

	/* --- gas species list --- */
	stream << "Order of gas species" << unixnewline;
	stream << "Name        Molar mass     Ediss          #Species" << unixnewline;
	for (n=frm.igstart; n<frm.igend; n++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-10s %-11.8lf    %-11.4le    %-3d",frm.species(n).c_str(), frm.wmass(n), frm.Ediss(n),n+1);
		stream << buffer << unixnewline;
	}
	stream << unixnewline;

	/* --- surface species list --- */
	stream << "Order of surface species for (number each consecutively)" << unixnewline;
	stream << "Name        Molar mass     Edes           Species#   Phase#/Site#" << unixnewline;
	k = frm.igend;
	for (j=0; j<frm.nsp; j++)	{
		for (m=0; m<frm.nspas(j); m++)	{
			for (n=0; n<frm.nspass(j,m); n++)	{
				Edes = 0.0;
				if (frm.kadsr(k) >= 0) Edes = frm.Edes( frm.kadsr(k) );
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-10s %-11.8lf    %-11.4le    %-3d           %3d/%-3d",frm.species(k).c_str(), frm.wmass(k), Edes, k+1, j+2, m+1);
				stream << buffer << unixnewline;
				k++;
			}
		}
	}
	stream << unixnewline;

	/* --- bulk species list --- */
	stream << "Order of bulk species (number each consecutively)" << unixnewline;
	stream << "Name        Molar mass     Mole fraction  Species#   Phase#" << unixnewline;
	for (j=0; j<frm.nbp; j++)	{
		for (n=0; n<frm.nbps(j); n++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-10s %-11.8lf    %-11.4le    %-3d           %-3d",frm.species(k).c_str(), frm.wmass(k), frm.bmf(j,n), k+1, j+2+frm.nsp);
			stream << buffer << unixnewline;
			k++;
		}
	}
	stream << unixnewline;

	/* --- reaction header --- */
	stream << "Total number of surface reactions" << unixnewline;
	stream << "nsrt" << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-4d",frm.isrt);
	stream << buffer << unixnewline;
	stream << unixnewline;

	/* --- reaction participants --- */
	stream << "Reactant/product species for each forward surface reaction" << unixnewline;
	for (i=0; i<frm.isrt; i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%2d,%2d,%2d,%2d,%2d,%2d      #%-3d",1+frm.RR(i,0),1+frm.RR(i,1),1+frm.RR(i,2),1+frm.RP(i,0),1+frm.RP(i,1),1+frm.RP(i,2),i+1);
		stream << buffer << unixnewline;
	}
	stream << unixnewline;

	/* --- reaction stoichiometry --- */
	stream << "Stoichiometric coefficients for each surface reaction" << unixnewline;
	for (i=0; i<frm.isrt; i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%2d,%2d,%2d,%2d,%2d,%2d      #%-3d",frm.RVR(i,0),frm.RVR(i,1),frm.RVR(i,2),frm.RVP(i,0),frm.RVP(i,1),frm.RVP(i,2),i+1);
		stream << buffer << unixnewline;
	}
	stream << unixnewline;

	/* --- forward reaction rate data --- */
	stream << "Reaction parameters for each type of reaction:" << unixnewline;
	stream << "Arrhenius:             0,   Cf, beta,   Ea, isrfon, isrbon" << unixnewline;
	stream << "Adsorption:            1,   S0, beta, Eads, isrfon, isrbon" << unixnewline;
	stream << "Eley-Rideal:           2,  Ger, beta,  Eer, isrfon, isrbon" << unixnewline;
	stream << "Langmuir-Hinschelwood: 3,  Clh, beta,   Em, isrfon, isrbon" << unixnewline;
	stream << "Sublimation:           4, a0Pv, beta, Esub, isrfon, isrbon" << unixnewline;
	stream << "Arrhenius Adsorption:  5,   Cf, beta,   Ea, isrfon, isrbon" << unixnewline;
	stream << "Type, (Param), beta,   Ea, isrfon, isrbon" << unixnewline;
	for (i=0; i<frm.isrt; i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-1d, %-10.4le, %-10.4le, %-10.4le, %-1d, %-1d     #%-3d",frm.istype(i), frm.Cf(i), frm.eta(i), frm.Ea(i), frm.isrfon(i), frm.isrbon(i), i+1);
		stream << buffer << unixnewline;
	}
	stream << unixnewline;
	
	/* --- desorption rate data --- */
	stream << "Desorption reaction or equilibrium constant parameters:" << unixnewline;
	stream << "Type 1: Desorption:" << unixnewline;
	stream << "        Form 0: Arrhenius" << unixnewline;
	stream << "        Form 1: Constant attempt frequency" << unixnewline;
	stream << "        Form 2: Simple transition state theory" << unixnewline;
	stream << "        Form 3: Complex transition state theory" << unixnewline;
	stream << "Type 2: Equilibrium:" << unixnewline;
	stream << "        Form 0: Arrhenius" << unixnewline;
	stream << "        Form 1: Immobile adsorption - simple transition state theory" << unixnewline;
	stream << "        Form 2: Immobile adsorption - complex transition state theory" << unixnewline;
	stream << "        Form 3: Mobile adsorption - simple transition state theory" << unixnewline;
	stream << "        Form 4: Mobile adsorption - complex transition state theory" << unixnewline;
	stream << "Type, Form, Cf, eta, vdes, Edes" << unixnewline;
	for (i=0; i<frm.isrt; i++)	{

		if ( (frm.istype(i) == 1)||(frm.istype(i) == 5) )	{

			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-1d, %-1d, %-10.4le, %-10.4le, %-10.4le, %-10.4le     #%-3d",frm.iadtype(i), frm.iadform(i), frm.Cfad(i), frm.etaad(i), frm.vad(i), frm.Edes(i), i+1);
			stream << buffer << unixnewline;

		}
	}
	stream << unixnewline;
	
	/* --- 
	   Extra data not read by input file
	   --- */

	stream << "     --------------------     " << unixnewline;

	/* --- version information --- */
	stream << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Finite Rate Model Version %1d-%-2d", FINITE_RATE_MODEL_MAJOR_VERSION, FINITE_RATE_MODEL_MINOR_VERSION);
	stream << buffer << unixnewline;
	
	/* --- derived surface parameters --- */
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Inert Surface Fraction   %7.4lf", frm.sfrcI);
	stream << buffer << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Surface Density           %10.4le", frm.sdentot);
	stream << buffer << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Inert Bulk Fraction      %7.4lf", frm.vfbI);
	stream << buffer << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Bulk density              %10.4le", frm.bdentot);
	stream << buffer << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Bulk Porosity            %7.4lf", frm.bportot);
	stream << buffer << unixnewline;
	stream << unixnewline;
	
	/* --- blowing data --- */
	if (frm.isblwon)	{
		stream << "     --------------------     " << unixnewline;
		stream << unixnewline;
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Blowing gas data for %s", frm.title);
		stream << buffer << unixnewline;
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "----------------------------------------------");
		stream << buffer << unixnewline;
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %1d      | EBC type (0: mdotg (kg/s)) (1: cyield) (2: cyield, ha (J/kg))",frm.nebc);
		stream << buffer << unixnewline;
		if (frm.nebc > 0)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %7.4lf  %12.5le      | yield/enthalpy",frm.cyield,frm.ha);
			stream << buffer << unixnewline;
		} else {
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %1d      | Number of blowing/pyrolyzing flows",frm.nblw);
			stream << buffer << unixnewline;
			for (n=0; n<frm.nblw; n++)	{
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %11.4le ",frm.rmdotg(n));
				stream << buffer;
			}
			stream << "           | each mdotg (kg/s)" << unixnewline;
		}
		for (k=frm.igstart; k<frm.igend; k++)	{
			for (n=0; n<frm.nblw; n++)	{
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %7.4lf     ",frm.blwmf(n,k));
				stream << buffer;
			}
			if (k==0)
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-10s | <= begin Mole fractions of blowing gas",frm.species(k).c_str());
			else
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%-10s",frm.species(k).c_str());
			stream << buffer << unixnewline;
		}
		stream << unixnewline;
	}
	
	stream << "     --------------------     " << unixnewline;
	stream << unixnewline;
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "Gibbs Free-energy data from Lewis database");
	stream << buffer << unixnewline;
	for (k=0; k<frm.nspecies; k++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %s %3d %3d %15.6le %15.6le %15.6le",frm.species(k).c_str(), frm.rlewis_nrange(k), 0, frm.wmass(k)*1000.0, 0.0, 0.0);
		stream << buffer << unixnewline;
		if (frm.kthermo(k) > 0) {
			for (j=0; j<frm.rlewis_nrange(k); j++)	{
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "    %14.3le %14.3le",frm.rlewis_trange(0,j,k), frm.rlewis_trange(1,j,k));
				stream << buffer << unixnewline;
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "    %14.6le %14.6le %14.6le %14.6le %14.6le", frm.rlewis_coeffs(0,j,k), frm.rlewis_coeffs(1,j,k), frm.rlewis_coeffs(2,j,k), frm.rlewis_coeffs(3,j,k), frm.rlewis_coeffs(4,j,k) );
				stream << buffer << unixnewline;
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "    %14.6le %14.6le %14.6le %14.6le %14.6le", frm.rlewis_coeffs(5,j,k), frm.rlewis_coeffs(6,j,k), frm.rlewis_coeffs(7,j,k), frm.rlewis_coeffs(8,j,k), frm.rlewis_coeffs(9,j,k) );
				stream << buffer << unixnewline;
			}
		} else {
			stream << "    Not using Gibbs data!" << unixnewline;
		}
	}

	return stream; //must return the stream 
}
/* ------------------------------------------------------------------------------- */



/* ----------------------------- print_etable ------------------------------------ */
/* 
 * print element table to stream
 *
 */
void FiniteRateSurfaceModel::compute_recession_rates(FiniteRateSurfaceCell& cell, double &sblow, double &pblow, double &recess)
{
	int k;
	double bb;
	Array1D<double> &srhs = cell.get_rhs();

	sblow  = 0.0;
	pblow  = 0.0;
	recess = 0.0;

	//surface blowing
	for (k=ibstart; k<ibend; k++)
		sblow -= wmass(k) * srhs(k);
	//pyrolysis blowing
	if (nblwflag == 1)	{
		if (nebc > 0)
			pblow = sblow * (1.0/cyield - 1.0);
		else {
			for (k=ibstart; k<ibend; k++)
				pblow += pyroadd(k);
		}
	}
	//recession
	bb = bdentot;
	if (bb < 1.0e-30) bb = 1.0e-30;
	recess = sblow / bb;

}
/* ------------------------------------------------------------------------------- */



/* ----------------------------- print_etable ------------------------------------ */
/* 
 * print element table to stream
 *
 */
void FiniteRateSurfaceModel::print_etable(ostream &stream)
{
	char buffer[INPUT_FILE_BUFFER_SIZE], unixnewline = 0x0A;
	int i, j;

	stream << "      Element List" << unixnewline;
	stream << "       ";
	for (i=0; i<nelements; i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-4s",elemnames(i).c_str());
		stream << buffer;
	}
	stream << unixnewline;
	
	stream << "      Stoichiometry" << unixnewline;
	for (j=0; j< nspecies; j++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "        %-10s:  ",species(j).c_str());
		stream << buffer;
		for (i=0; i<nelements; i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %-3d",(int)(etable(j,i)));
			stream << buffer;
		}
		stream << unixnewline;
	}

}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- print_diagnostics -------------------------------- */
/* 
 * print out additional diagnostics itemizing rates
 *
 */
void FiniteRateSurfaceModel::print_diagnostics(ostream &stream, FiniteRateSurfaceCell& cell, double tmpr, int iscan)
{
	Array1D<double> &srhs = cell.get_rhs();
	Array2D<double> &sjacobian = cell.get_jacobian();
	Array1D<double> &concentration = cell.get_concentration();

	int ispec,ireact,iphase,ibulk,index,j,ier,k;
	double  tmp, prate, rflux, pf, pb, dFdtg, dFdtw, eqact;
	Array1D<double> frates(isrt), brates(isrt), eqa(isrt), eqc(isrt),
			    df_dtg(isrt), df_dtw(isrt), deqa_dtg(isrt), deqa_dtw(isrt),
			    deqc_dtg(isrt), deqc_dtw(isrt), db_dtg(isrt), db_dtw(isrt);
	Array1D<double> gibbs(nspecies), dg_dtw(nspecies), dg_dtg(nspecies);
	
	char buffer[INPUT_FILE_BUFFER_SIZE], unixnewline = 0x0A;  //LineFeed
	
	/* --- zero cell arrays --- */
	srhs.initialize(0.0);
	sjacobian.initialize(0.0);
	
	/* --- initialize --- */
	frates.initialize(0.0);
	brates.initialize(0.0);
	eqa.initialize(0.0);
	eqc.initialize(0.0);
	df_dtg.initialize(0.0);
	df_dtw.initialize(0.0);
	deqa_dtg.initialize(0.0);
	deqa_dtw.initialize(0.0);
	deqc_dtg.initialize(0.0);
	deqc_dtw.initialize(0.0);
	db_dtg.initialize(0.0);
	db_dtw.initialize(0.0);
	gibbs.initialize(0.0);
	dg_dtw.initialize(0.0);
	dg_dtg.initialize(0.0);

	dFdtg = 0.0;
	dFdtw = 0.0;

	/* --- compute forward rates --- */
	for (ireact=0; ireact<isrt; ireact++)	{
		
		// identify gas phase participant
		k = RR(ireact,0);
		if (istype(ireact) == 4) k = RP(ireact,0);
		//compute forward rate
		frates(ireact) = srate(tmpr, tmpr, wmass(k), sdensity(ireact), sitep(ireact), istype(ireact),
						Cf(ireact), eta(ireact), Eact(ireact), &(df_dtg(ireact)), &(df_dtw(ireact)));
	}

	/* --- obtain Gibbs energies for all species --- */
	compute_gibbs_energies(gibbs, dg_dtg, dg_dtw, tmpr, tmpr, frates, df_dtg, df_dtw);

	/* --- calculate activity based equilibrium constants and derivatives --- */
	compute_eqcon(eqc, deqc_dtg, deqc_dtw, eqa, deqa_dtg, deqa_dtw, tmpr, tmpr, gibbs, dg_dtg, dg_dtw);

	/* --- calculate backward rates & derivatives --- */
	for (ireact=0; ireact<isrt; ireact++)	{
		brates(ireact) = frates(ireact) / eqc(ireact);
		db_dtg(ireact) = (df_dtg(ireact) - brates(ireact)*deqc_dtg(ireact))/eqc(ireact);
		db_dtw(ireact) = (df_dtw(ireact) - brates(ireact)*deqc_dtw(ireact))/eqc(ireact);
	}


	/* --- write reaction information to screen --- */
	if (iscan==0) stream << "      <<Reaction Constants>>" << unixnewline;
	if (iscan==0) stream << "      react      kf          kb         KeqC     KeqC-actual   wf-actual   wb-actual" << unixnewline;
	for (ireact=0; ireact<isrt; ireact++)	{
		//sum forward and backward concentration product term
		pf = frates(ireact) * ((double)(isrfon(ireact))) * sfraction(ireact);
		pb = brates(ireact) * ((double)(isrbon(ireact))) * sfraction(ireact);
		eqact = 1.0;
		for (ispec=0; ispec<nspecies; ispec++)	{
			pf    = pf * pow(concentration(ispec), vr(ireact,ispec));
			pb    = pb * pow(concentration(ispec), vp(ireact,ispec));
			eqact = eqact * pow(concentration(ispec), v(ireact,ispec));
		}
		//keep track of species RHS
		for (ispec=0; ispec<nspecies; ispec++)
			srhs(ispec) += v(ireact,ispec) * (pf - pb);
		//write reaction diagnostic statement
		if (iscan==0) {
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "       %3d  %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e",ireact+1,frates(ireact),brates(ireact),eqc(ireact),eqact,pf,pb);
			stream << buffer << unixnewline;
		}
	}


	/* --- write species global production information to screen --- */
	if (iscan==0) {
		stream << "      <<Species Production>>" << unixnewline;
		stream << "      species      concentration     molar prod.      mass prod." << unixnewline;
		for (ispec=0; ispec<nspecies; ispec++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "      %-10s %15.4e %15.4e %15.4e",species(ispec).c_str(), concentration(ispec), srhs(ispec),  wmass(ispec)*srhs(ispec));
			stream << buffer << unixnewline;
		}
	}


	/* --- write reaction specific production information to screen --- */
	if (iscan==0) {
		stream << "      <<Species Production Per Reaction>>" << unixnewline;
		stream << "      species           w-dot";
		for (ireact=0; ireact<isrt; ireact++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "           w%03d",ireact+1);
			stream << buffer;
		}
		stream << unixnewline;
	}
		
	for (ispec=0; ispec<nspecies; ispec++)	{
		if (iscan==0) snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "      %-10s %12.4e",species(ispec).c_str(), srhs(ispec));
		if (iscan==0) stream << buffer;
		for (ireact=0; ireact<isrt; ireact++)	{
			pf = frates(ireact) * ((double)(isrfon(ireact))) * sfraction(ireact);
			pb = brates(ireact) * ((double)(isrbon(ireact))) * sfraction(ireact);
			for (k=0; k<nspecies; k++)	{
				pf    = pf * pow(concentration(k), vr(ireact,k));
				pb    = pb * pow(concentration(k), vp(ireact,k));
			}
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "   %12.4e",v(ireact,ispec)*(pf-pb));
			stream << buffer;
		}
		if (iscan==0) stream << unixnewline;
	}

}
/* ------------------------------------------------------------------------------- */



/* --------------------------- compute_gibbs_energies ---------------------------- */
/* 
 * compute Gibbs energies of all species
 *
 */
void FiniteRateSurfaceModel::compute_gibbs_energies(Array1D<double> &gibbs, Array1D<double> &dgdTg, Array1D<double> &dgdTw,
					 double Tw, double Tg, Array1D<double> const &frate, Array1D<double> const &dfdTg, Array1D<double> const &dfdTw)
{
	int n, ir, kg;
	double h, s, dh, ds, eqKa, dKag, dKaw;

	/* --- compute all gas phase species using Tg --- */
	for (n=igstart; n<igend; n++)	{
		
		//get Gibbs data for this species
		get_gibbs_data_2(rlewis_coeffs, rlewis_trange, Tg, n, &h, &s, &dh, &ds);
		
		//Gibbs energies made dimensional to avoid temperature vagueness
		gibbs(n) = (h - s)*R*Tg;
		dgdTg(n) = (dh - ds)*R*Tg + (h - s)*R;
		dgdTw(n) = 0.0;
		
	}
	
	/* --- surface/bulk species calculated using Tw or interpolated --- */
	for (n=igend; n<nspecies; n++)	{
		
		if (kthermo(n) == 1)	{

			//get Gibbs data for this species
			get_gibbs_data_2(rlewis_coeffs, rlewis_trange, Tw, n, &h, &s, &dh, &ds);
		
			//Gibbs energies made dimensional to avoid temperature vagueness
			gibbs(n) = (h - s)*R*Tw;
			dgdTg(n) = 0.0;
			dgdTw(n) = (dh - ds)*R*Tw + (h - s)*R;
		
		} else if (kthermo(n) == 2)	{
		
			gibbs(n) = 0.0;
			dgdTg(n) = 0.0;
			dgdTw(n) = 0.0;
		
		} else {
		
			ir = kadsr(n);  //get adsorption reaction #
			kg = RR(ir,0);  //get gas species that anchors this reaction
			eqKa = adeq(Tw, Tg, wmass(kg), iadtype(ir), iadform(ir), Cfad(ir), etaad(ir), vad(ir), Edes(ir), Ea(ir), frate(ir), dfdTg(ir), dfdTw(ir), &dKag, &dKaw);
			
			gibbs(n) = gibbs(kg) - R*Tw*log(eqKa);
			dgdTg(n) = dgdTg(kg) - R*Tw/eqKa*dKag;
			dgdTw(n) = dgdTw(kg) - R*(log(eqKa) + Tw/eqKa*dKaw);
		
		}
	
	}

}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- compute_eqcon ------------------------------------ */
/* 
 * compute equilibrium constants for all reactions
 *
 */
void FiniteRateSurfaceModel::compute_eqcon(Array1D<double> &EqC, Array1D<double> &dEqCdTg, Array1D<double> &dEqCdTw,
					 Array1D<double> &EqA, Array1D<double> &dEqAdTg, Array1D<double> &dEqAdTw,
					 double Tw, double Tg, Array1D<double> const &gibbs, Array1D<double> const &dgdTg, Array1D<double> const &dgdTw)
{
	int ireact, ispec, iphase;
	double tmp;
	
	/* --- initialize --- */
	EqA.initialize(0.0);
	dEqAdTg.initialize(0.0);
	dEqAdTw.initialize(0.0);
	EqC.initialize(0.0);
	dEqCdTg.initialize(0.0);
	dEqCdTw.initialize(0.0);
	
	/* --- calculate activity based equilibrium constants --- */
	for (ireact=0; ireact<isrt; ireact++)	{
		for (ispec=0; ispec<nspecies; ispec++)	{
			tmp = v(ireact,ispec)/(R*Tw);
			EqA(ireact) = EqA(ireact) - tmp*gibbs(ispec);
			dEqAdTg(ireact) = dEqAdTg(ireact) - tmp * dgdTg(ispec);
			dEqAdTw(ireact) = dEqAdTw(ireact) - tmp * (dgdTw(ispec) - gibbs(ispec)/Tw);
		}
		EqA(ireact) = exp(EqA(ireact));
		dEqAdTg(ireact) = EqA(ireact) * dEqAdTg(ireact);
		dEqAdTw(ireact) = EqA(ireact) * dEqAdTw(ireact);
	}

	/* --- calculate concentration based equilibrium constants --- */
	for (ireact=0; ireact<isrt; ireact++)	{
		tmp = pow(Pref/(R*Tw),sumv(ireact,0));
		for (iphase=1; iphase<nsp+1; iphase++)
			tmp = tmp * pow(sden(iphase-1),sumv(ireact,iphase));
		EqC(ireact)     = tmp * EqA(ireact);
		dEqCdTg(ireact) = tmp * dEqAdTg(ireact);
		dEqCdTw(ireact) = tmp * (dEqAdTw(ireact) - sumv(ireact,0) * EqA(ireact)/Tw);
	}

}
/* ------------------------------------------------------------------------------- */



/* ------------------------------- srate ----------------------------------------- */
/* 
 * compute forward reaction rate
 *
 * itype=0   Arrhenius 		
 * type=1   Adsoprtion
 * itype=2   Eley-Rideal
 * itype=3   Langmuir-Hinschelwood
 * itype=4   Sublimation based on vapor pressure
 * itype=5   Arrhenius adsorption
 *
 */
double FiniteRateSurfaceModel::srate(double Tw, double Tg, double wm, double stot, double p, int itype, double A, double B, double E, double *dkdTg, double *dkdTw)
{

	double rate,v,S0c,S0,gERc,gER,v2D,clh,a0vp0;

	*dkdTg = 0.0;
	*dkdTw = 0.0;
	rate   = 0.0;

	switch(itype)	{

	case 0:	//Arrhenius reaction
	case 5:

		rate   = A * pow(Tw, B) * exp(-E/(R*Tw));
		*dkdTw = rate * (B + E/(R*Tw))/Tw;

		break;
	case 1:	//Adsorption Reaction: A=sticking coefficient, E=Adsorption Energy barrier
	
		v=sqrt(8.0*R*Tg/(PI*wm));
		S0c=A*pow(Tw,B);
		S0=min(1.0,S0c);
		rate=(v/(4.0*pow(stot,p)))*S0*exp(-E/(Tw*R));
		if (S0c > 1.0)
			*dkdTw = rate * E/(Tw*R*Tw);
		else
			*dkdTw = rate * (B+E/(Tw*R))/Tw;
		*dkdTg = rate * 0.5/Tg;
	
		break;
	case 2:	//Eley-Rideal Reaction: A=reaction efficiency, E=Reaction Energy barrier

		v=sqrt(8.0*R*Tg/(PI*wm));
		gERc=A*pow(Tw,B);
		gER=min(1.0,gERc);
		rate=(v/(4.0*pow(stot,p)))*gER*exp(-E/(Tw*R));
		if (gERc > 1.0)
			*dkdTw = rate * E/(Tw*Tw*R);
		else
			*dkdTw = rate * (B+E/(Tw*R))/Tw;
		*dkdTg = rate * 0.5/Tg;

		break;
	case 3:	//Langmuir-Hinschelwood Reaction: A=reaction efficiency, E=Reaction barrier

		v2D=sqrt(PI*R*Tw/(2.0*wm));
		clh = A*pow(Tw,B);
		rate=v2D*sqrt(AV*pow(stot,(3.0-2.0*p)))*clh*exp(-E/(Tw*R));
		*dkdTw = rate * (0.5+B+E/(Tw*R))/Tw;

		break;
	case 4:	//Sublimation Reaction: based on vapor pressure

		v=sqrt(8.0*R*Tg/(PI*wm));
		a0vp0=(A*pow(Tw,B))/(R*Tg);
		rate=(v/(4.0*pow(stot,p)))*a0vp0*exp(-E/(Tw*R));
		*dkdTw = rate * (B + E/(Tw*R))/Tw;
		*dkdTg = -0.5*rate/Tg;

		break;
	}

	return rate;
}
/* ------------------------------------------------------------------------------- */



/* ------------------------------- adeq ------------------------------------------ */
/* 
 * perform sanity checks on the model data
 *
 * Type 1  Desorption
 *  Form 0:  Arrhenius
 *  Form 1:  constant frequency
 *  Form 2:  Simple TST frequency
 *  Form 3:  Complex TST frequency
 * 
 * Type 2  Equilibrium constant
 *  Form 0:  Arrhenius
 *  Form 1:  Immobile simple TST
 *  Form 2:  Immobile complex TST
 *  Form 3:  Mobile simple TST
 *  Form 4:  Mobile complex TST
 *
 */
double FiniteRateSurfaceModel::adeq(double Tw, double Tg, double wm, int itype, int iform, double A, double B, double v,
						double Ed, double Ea, double frate, double dfdTg, double dfdTw, double *dEQdTg, double *dEQdTw)
{
	double brate, vfac, eqKc, qtr, qvib, qtr1, qtr2, dkbdTw, dvf, eterm, eqKa, dKCdTw, dKCdTg;
	
	brate   = 1.0;
	dkbdTw  = 0.0;
	
	if (itype == 1)	{
	
		switch (iform)	{
		case 0:	//Arrhenius

			brate = A * pow(Tw,B) * exp(-Ed/(R*Tw));
			dkbdTw = brate * (B + Ed/(R*Tw))/Tw;

			break;
		case 1:	//constant attempt frequency

			brate = A * pow(Tw,B) * v * exp(-Ed/(R*Tw));
			dkbdTw = brate * (B + Ed/(R*Tw))/Tw;

			break;
		case 2:	//simple TST frequency

			brate = A * pow(Tw,B) * (R*Tw/(AV*H)) * exp(-Ed/(R*Tw));
			dkbdTw = brate * (B + 1.0 + Ed/(R*Tw))/Tw;

			break;
		case 3:	//complex TST frequency
			
			eterm  = exp(R*Tw/(AV*H*v));
			vfac   = sqrt(eterm)/(.0-eterm);
			brate  = A * pow(Tw,B) * (R*Tw/(AV*H)) / vfac * exp(-Ed/(R*Tw));
			dvf    = 0.5 * R / (AV*H*v) * (1.0+eterm)/(1.0-eterm);
			dkbdTw = brate * ((B + 1.0 + Ed/(R*Tw))/Tw - dvf);
			
			break;
		}
		
		eqKc   = frate / brate;
		dKCdTw = dfdTw/brate - frate*dkbdTw/(brate*brate);
		dKCdTg = dfdTg/brate;
	
	} else if (itype == 2)	{

		switch (iform)	{
		case 0:	//Arrhenius

			eqKc   = A * pow(Tw,B) * exp(Ed/(R*Tw));
			dKCdTw = eqKc * (B - Ed/(R*Tw))/Tw;

			break;
		case 1:	//Immobile simple TST

			qtr    = pow( 2.0*PI*wm*R*Tg/(AV*AV*H*H) , -1.5 );
			eqKc   = A * pow(Tw,B) * qtr * exp(-Ea/(R*Tw)) * exp(Ed/(R*Tw));
			dKCdTw = eqKc * (B + (Ea-Ed)/(R*Tw))/Tw;
			dKCdTg = -1.5 * eqKc / Tg;

			break;
		case 2:	//Immobile complex TST

			eterm  = exp( -H*AV*v/(R*Tw));
			qtr    = pow( 2.0*PI*wm*R*Tg/(AV*AV*H*H) , -1.5 );
			qvib   = pow( sqrt(eterm)/(1.0-eterm), 3.0);
			eqKc   = A * pow(Tw,B) * qtr * qvib * exp(-Ea/(R*Tw)) * exp(Ed/(R*Tw));
			dKCdTw = eqKc * (B + (Ea-Ed)/(R*Tw) + 1.5*(H*AV*v/(R*Tw))*(1.0+eterm)/(1.0-eterm))/Tw;
			dKCdTg = -1.5 * eqKc / Tg;

			break;
		case 3:	//Mobile simple TST
			
			qtr1   = pow( 2.0*PI*wm*R*Tg/(AV*AV*H*H) , -1.5 );
			qtr2   = 2.0*PI*wm*R*Tw/(AV*AV*H*H);
			eqKc   = A * pow(Tw,B) * qtr1 * qtr2 * exp(-Ea/(R*Tw)) * exp(Ed/(R*Tw));
			dKCdTw = eqKc * (B + 1.0 + (Ea-Ed)/(R*Tw))/Tw;
			dKCdTg = -1.5 * eqKc / Tg;
			
			break;
			
		case 4:	//Mobile complex TST
			
			qtr1   = pow( 2.0*PI*wm*R*Tg/(AV*AV*H*H) , -1.5 );
			qtr2   = 2.0*PI*wm*R*Tw/(AV*AV*H*H);
			eterm  = exp( -H*AV*v/(R*Tw));
			qvib   = sqrt(eterm)/(1.0-eterm);
			eqKc   = A * pow(Tw,B) * qtr1 * qtr2 * qvib * exp(-Ea/(R*Tw)) * exp(Ed/(R*Tw));
			dKCdTw = eqKc * (B + 1.0 + (Ea-Ed)/(R*Tw) + 0.5*(H*AV*v/(R*Tw))*(1.0+eterm)/(1.0-eterm))/Tw;
			dKCdTg = -1.5 * eqKc / Tg;
			
			break;
		}

	}

	// Calculate activity-based equilibrium constant
      eqKa   = eqKc*Pref/(R*Tg);
      *dEQdTw = dKCdTw*Pref/(R*Tg);
      *dEQdTg = dKCdTg*Pref/(R*Tg) - eqKc*Pref/(R*Tg*Tg);

	return eqKa;
}
/* ------------------------------------------------------------------------------- */



/* ----------------------------- build_etable ------------------------------------ */
/* 
 * build elemental stoichiometry table
 *
 */
int FiniteRateSurfaceModel::build_etable()
{
	const char *name;
	string *tmp_elem_names=NULL;
	Array2D<int> tmp_elem_matrix;
	int ntmpelem, i, j, iclose, ipopen, iqstart, iqend, isurf, istart, iend, ierror,
	    is, ic, ncurrent, ncharge, itype;

//	----- initialize error -----
	ierror = 0;

//     ----- precompute the number of possible species that may occur -----
//           n = 1 to start to account for electrons
//           anything with 's' might be a surface site
	ntmpelem = 1;
	for (i=0; i<nspecies; i++)	{
		name = species(i).c_str();
		for (j=0; j<strlen(name); j++)	{
			if (name[j] == 's')	ntmpelem++;
			if (isupper(name[j]))	ntmpelem++;
		}
	}

//	----- resize temporary arrays -----
	tmp_elem_matrix.resize(nspecies, ntmpelem);
	tmp_elem_matrix.initialize(0);
	tmp_elem_names = new string[ntmpelem];

//	----- initialization -----
	nelements = 0;
	
	iclose  = 0;
      ipopen  = -1;
      istart  = -1;
      iend    = -1;
      iqstart = -1;
      iqend   = -1;
      isurf   = 0;

//	----- cycle through all characters looking for unique elements -----
	for (is=0; is<nspecies; is++)	{

		// easier to access C-strings
		name = species(is).c_str();

		//go through all characters in the string
		for (ic=0; ic< species(is).length(); ic++)	{
		
				
			//character type (copied from Fortran code)
			itype = 0;
			if (isupper(name[ic])) itype = 1;
			if (islower(name[ic])) itype = 2;
			if (isdigit(name[ic])) itype = 3;
			if ( name[ic] == '+' ) itype = 4;
			if ( name[ic] == '-' ) itype = 5;
			if ( name[ic] == '(' ) itype = 6;
			if ( name[ic] == ')' ) itype = 7;
			if ( name[ic] == ' ' ) itype = 8;

			// special consideration for open parenthesis
			if (ipopen >= 0)	{
			
				if (itype == 7)	{
					if (ic > (ipopen+1))	{
						if (name[ipopen+1] != 'b')	{
							ncurrent = matchelemname(nelements, tmp_elem_names, name+ipopen+1, ic-ipopen-1);
							tmp_elem_matrix(is,ncurrent) = 1;
						}
						ipopen = -1;
					} else {
						ipopen = -1;  //empty parenthesis
					}
				}
			
			} else {
				switch (itype)	{
				
				case (1):  //upper case - start a new element
					if (istart >= 0)	{
						ncurrent = matchelemname(nelements, tmp_elem_names, name+istart, iend-istart+1);
						tmp_elem_matrix(is,ncurrent) = 1;
					}
					
					if (name[ic] != 'E')	{
						istart = ic;
						iend   = ic;
					} else
						istart = -1;
					
					break;
				case (2):  //lower case - continue
				
					if (istart >= 0) iend++;
					else	{
						if (name[ic] == 'e')	{  //matched free-electron
							ncharge = matchelemname(nelements, tmp_elem_names, "CHARGE", 6);
							tmp_elem_matrix(is,ncharge) = tmp_elem_matrix(is,ncharge) - 1;
						} else  //error: no other element starts lowercase
							ierror = 2;
					}
				
					break;
				case (3):  //number - quantity info
					if (iqstart < 0)	{
						iqstart = ic;
						iqend   = ic;
					} else
						iqend++;
					break;
				case (4):  //plus or minus - charge
				case (5):
					if (istart >= 0)	{
						ncurrent = matchelemname(nelements, tmp_elem_names, name+istart, iend-istart+1);
						tmp_elem_matrix(is,ncurrent) = 1;
						istart = -1;
					}
					ncharge = matchelemname(nelements, tmp_elem_names, "CHARGE", 6);
					if (itype==4)
						tmp_elem_matrix(is,ncharge) = tmp_elem_matrix(is,ncharge) + 1;
					else
						tmp_elem_matrix(is,ncharge) = tmp_elem_matrix(is,ncharge) - 1;
					
					break;
				case (6):  //left parenthesis (open)
					if (ipopen>=0)	ierror = 6;
					else			ipopen = ic;
					break;
				case (7):  //right parenthesis (close)
					if (ipopen < 0)	ierror = 7;
					else			ipopen = -1;
					break;
				default:
					break;
				
				}
			}
			
					
			if ( (istart >= 0)&&(itype>2) ) iclose = 1;
			if ( (istart >= 0)&&(ic==(species(is).length()-1)) ) iclose = 1;
			
			if (iclose)	{
				ncurrent = matchelemname(nelements, tmp_elem_names, name+istart, iend-istart+1);
				tmp_elem_matrix(is,ncurrent) = 1;
				iclose = 0;
				istart = -1;
			}
			
			if ((iqstart >= 0)&&(itype != 3)) iclose = 1;
			if ((iqstart >= 0)&&(ic==(species(is).length()-1))) iclose = 1;
			
			if (iclose)	{
			
				tmp_elem_matrix(is,ncurrent) = atoi(name+iqstart);
				iclose = 0;
				iqstart = -1;
			}
			
		}  //end ic loop
	
	} //end is loop

	etable.resize(nspecies,nelements);
	etable.initialize(0.0);
	elemnames.resize(nelements);
	for (j=0; j<nelements; j++)	{
		elemnames(j) = tmp_elem_names[j];
		for (i=0; i<nspecies; i++)
			etable(i,j) = (double)(tmp_elem_matrix(i,j));
	}
	
	if (tmp_elem_names) delete [] tmp_elem_names;
	
	return ierror;
}
/* ------------------------------------------------------------------------------- */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                    PRIVATE FUNCTION DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


/* --------------------------- get_line_of_text ---------------------------------- */
/* 
 * get a line of text searching for redundant newline characters
 *
 */
void get_line_of_text(istream &stream, char *buffer, int size)
{
	int ipos=0;
	char c='\0', c2='\0';
	bool flag = false;

	while ( (flag==false) && (ipos<(size-1)) )	{
		c='\0';
		stream.peek();
		if (stream.rdstate() == 0)	{
			stream.get(c);
			if ( (c == 0x0A) || (c == 0x0D) )	{
				flag = true;
				c2 = stream.peek();
				if ( (stream.rdstate()==0) && ((c2 == 0x0A)||(c2 == 0x0D)) && (c2 != c) ) stream.get(c2);
			} else	{
				buffer[ipos] = c;
				ipos++;
			}
		} else
			flag = true;
	}
	
	buffer[ipos] = '\0';  //NULL terminator
}
/* ------------------------------------------------------------------------------- */



/* ------------------------- string_replace_char --------------------------------- */
/* 
 * replace all occurances of one character with another
 *
 */
void string_replace_char(char *str, char search_for, char replace_with)
{
	int ipos=0;
	
	if (str == NULL) return;
	if ((replace_with == '\0')||(search_for == '\0')) return;
	
	while (str[ipos] != '\0')	{
		
		if (str[ipos] == search_for) str[ipos] = replace_with;
		ipos++;

	}
}
/* ------------------------------------------------------------------------------- */



/* ----------------------------- string_is_end ----------------------------------- */
/* 
 * determine if a string is "END" with case insensitivity
 *
 */
int string_is_end(char *str)
{
	int length = strlen(str);
	
	if (length != 3) return 1;
	if ((str[0] != 'e') && (str[0] != 'E')) return 1;
	if ((str[1] != 'n') && (str[1] != 'N')) return 1;
	if ((str[2] != 'd') && (str[2] != 'D')) return 1;
	return 0;
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- get_gibbs_data ----------------------------------- */
/* 
 * get gibbs data from Lewis arrays (form 1)
 *
 */
void get_gibbs_data(Array3D<double> &data, double t, int irange, int ispecies, double *h, double *s, double *cp)
{
	double t1, t2, t3, t4, ti, ts, tl;

	t1 = t;
	t2 = t1*t1;
	t3 = t1*t2;
	t4 = t2*t2;
	ti = 1.0 / t1;
	ts = ti*ti;
	tl = log(t1);
	
	*h = data(2,irange,ispecies)
	   + data(1,irange,ispecies) * ti * tl
	   - data(0,irange,ispecies) * ts
	   + 0.5*data(3,irange,ispecies) * t1
	   + 1.0/3.0*data(4,irange,ispecies) * t2
	   + 0.25*data(5,irange,ispecies) * t3
	   + 0.2*data(6,irange,ispecies) * t4
	   + data(8,irange,ispecies) * ti ;

	*s = data(2,irange,ispecies) * tl
	   - data(1,irange,ispecies) * ti
	   - 0.5*data(0,irange,ispecies) * ts
	   + data(3,irange,ispecies) * t1
	   + 0.5*data(4,irange,ispecies) * t2
	   + 1.0/3.0*data(5,irange,ispecies) * t3
	   + 0.25*data(6,irange,ispecies) * t4
	   + data(9,irange,ispecies) ;

	*cp = data(0,irange,ispecies) * ts
	    + data(1,irange,ispecies) * ti
	    + data(2,irange,ispecies)
	    + data(3,irange,ispecies) * t1
	    + data(4,irange,ispecies) * t2
	    + data(5,irange,ispecies) * t3
	    + data(6,irange,ispecies) * t4;

}
/* ------------------------------------------------------------------------------- */



/* --------------------------- get_gibbs_data_2 ---------------------------------- */
/* 
 * get gibbs data from Lewis arrays with derivative (form 2)
 *
 */
void get_gibbs_data_2(Array3D<double> &data, Array3D<double> &trange, double t, int ispecies, double *h, double *s, double *dh, double *ds)
{
	double t1, t2, t3, t4, ti, ts, tl;
	int irange, i;
	
	// determine which range we fall into; we don't need nrange for this since trange is initialized to zero
	irange = 0;
	for (i=0; i<NLEWIS_FITS; i++)	{
		if ( (t >= trange(0,i,ispecies)) && (t <= trange(1,i,ispecies)) ) irange = i;
	}
	
	//compute temperature parameters
	t1 = t;
	t2 = t1*t1;
	t3 = t1*t2;
	t4 = t2*t2;
	ti = 1.0 / t1;
	ts = ti*ti;
	tl = log(t1);

	//compute enthalpy
	*h = data(2,irange,ispecies)
	   + data(1,irange,ispecies) * ti * tl
	   - data(0,irange,ispecies) * ts
	   + 0.5*data(3,irange,ispecies) * t1
	   + 1.0/3.0*data(4,irange,ispecies) * t2
	   + 0.25*data(5,irange,ispecies) * t3
	   + 0.2*data(6,irange,ispecies) * t4
	   + data(8,irange,ispecies) * ti ;

	*dh = data(1,irange,ispecies) * ts * (1.0-tl)
	    + 2.0*data(0,irange,ispecies) * ts * ti
	    + 0.5*data(3,irange,ispecies)
	    + 2.0/3.0*data(4,irange,ispecies) * t1
	    + 0.75*data(5,irange,ispecies) * t2
	    + 0.8*data(6,irange,ispecies) * t3
	    - data(8,irange,ispecies) * ts ;

	//compute entropy
	*s = data(2,irange,ispecies) * tl
	   - data(1,irange,ispecies) * ti
	   - 0.5*data(0,irange,ispecies) * ts
	   + data(3,irange,ispecies) * t1
	   + 0.5*data(4,irange,ispecies) * t2
	   + 1.0/3.0*data(5,irange,ispecies) * t3
	   + 0.25*data(6,irange,ispecies) * t4
	   + data(9,irange,ispecies) ;

	*ds = data(0,irange,ispecies) * ts * ti
	    + data(1,irange,ispecies) * ts
	    + data(2,irange,ispecies) * ti
	    + data(3,irange,ispecies)
	    + data(4,irange,ispecies) * t1
	    + data(5,irange,ispecies) * t2
	    + data(6,irange,ispecies) * t3 ;
}
/* ------------------------------------------------------------------------------- */



/* ------------------------------ matchelemname ---------------------------------- */
/* 
 * search through elemental names looking for a match
 *
 */
int matchelemname(int &nelements, string *elem_names, const char *name, const int name_length)
{
	int match, i;
	char *local_name;
	
	//must copy string because the source is a part of a larger one
	local_name = new char [name_length+1];
	strncpy(local_name,name,name_length);
	local_name[name_length] = '\0';
	//default to a new element
	match = nelements;
	for (i=0; i<nelements; i++)
		if ( strcmp(local_name,elem_names[i].c_str()) == 0 ) match = i;
	
	if (match == nelements)	{
		elem_names[nelements].assign(local_name);
		nelements++;
	}
	
	//reclaim memory
	delete [] local_name;

	return match;
}
/* ------------------------------------------------------------------------------- */



/* --------------------------------- simplex ------------------------------------- */
/* 
 * search through elemental names looking for a match
 *
 */
int simplex(Array2D<double> &ebal, Array1D<double> &gibbs, Array1D<int> &ibasis)
{
	int nelements, nspecies, nrows, ncols, i, j, k, iprow, ipcol, nactive, iflag;
	double tmp;
	Array2D<double> rsys;
	Array1D<int>    isactive;
	
	//--- sizing ---
	nelements = ebal.getRows();
	nspecies  = ebal.getCols();
//printf("SIMPLEX nelements=%d nspecies=%d\n",nelements,nspecies);

	if (gibbs.getRows() != nspecies) return 1;
	if (ibasis.getRows() != nelements) return 2;
	
	//--- initialize & allocate ---
	ibasis.initialize(-1);
	
	isactive.resize(nelements);
	isactive.initialize(0);
	
	nactive = 0;
	for (i=0; i<nelements; i++)	{
		for (j=0; j<nspecies; j++)
			if (ebal(i,j)>0) isactive(i) = 1;
		if (isactive(i)>0) nactive++;
	}
	
	nrows = nactive + 1 + nspecies + 2;
	ncols = nspecies + 2*nactive + 1 + nspecies + 1;
	
	rsys.resize(nrows,ncols);
	rsys.initialize(0.0);
/*printf("nelements=%d nactive=%d nspecies=%d\n",nelements,nactive,nspecies);	
for (i=0; i<nelements; i++)	{
		for (j=0; j<nspecies; j++)
			printf(" %d ",(int)(ebal(i,j)));
		printf("\n");
	}
	for (j=0; j<nspecies; j++)
		printf(" gibbs = %lf\n",gibbs(j));

	for (i=0; i<nelements; i++)
		printf(" isactive = %d\n",isactive(i));
*/
	//--- build system with constraints ---

	//constraints: elemental composition
	k = 0;
	for (i=0; i<nelements; i++)	{
		if (isactive(i) > 0)	{
			for (j=0; j<nspecies; j++)
				if (ebal(i,j) > 0) rsys(k,j) = 1.0;
			rsys(k,nspecies+k)           = -1.0;
			rsys(k,2*nspecies+nactive+k) = 1.0;
			rsys(k,ncols-1)              = 1.0;
			k++;
		}
	}
	
	//constraint: sum[xj] = Ne
	for (j=0; j<nspecies; j++)
		rsys(nactive,j) = 1.0;
	rsys(nactive,ncols-2) = 1.0;
	rsys(nactive,ncols-1) = (double)(nactive);
	
	//constraints: xj <= 1
	for (j=0; j<nspecies; j++)	{
		rsys(nactive+1+j,j) = 1.0;
		rsys(nactive+1+j,nspecies+nactive+j) = 1.0;
		rsys(nactive+1+j,ncols-1) = 1.0;
	}
	
	//cost function: sum[gibbs(j)*xj]
	tmp = 0.0;
	for (j=0; j<nspecies; j++)
		if (tmp < gibbs(j)) tmp = gibbs(j);
	for (j=0; j<nspecies; j++)
		rsys(nrows-2,j) = gibbs(j) - tmp - 1.0;
	
	//artificial cost function: sum[zk]
	for (i=0; i<=nactive; i++)	{
		for (j=0; j<(nspecies+nactive); j++)
			rsys(nrows-1,j) -= rsys(i,j);
		rsys(nrows-1,ncols-1) -= rsys(i,ncols-1);
	}

/*      Concentrationsprintf(" ---------------------------------------------------\n");
for (i=0; i<nrows; i++)	{
for (j=0; j<ncols; j++)
	printf(" %6.2lf",rsys(i,j));
printf("\n");
}*/
	
	//--- perform phase I simplex ---
	for (k=0; k<10*ncols; k++)	{
		ipcol = -1;
		iprow = -1;
		tmp = 0.0;
		
		//pivot column is most negative artificial cost function
		for (j=0; j<(ncols-1); j++)	{
			if (rsys(nrows-1,j) < tmp)	{
				tmp = rsys(nrows-1,j);
				ipcol = j;
			}
		}
		if (ipcol < 0) break;
		
		//pivot row has the smallest positive ratio
		tmp = 0.0;
		for (j=0; j<(nrows-2); j++)	{
			if ( (rsys(j,ipcol)>0.0) && (rsys(j,ncols-1)>=0.0) )	{
				if ( (tmp > rsys(j,ncols-1)/rsys(j,ipcol)) || (iprow<0) )	{
					tmp = rsys(j,ncols-1) / rsys(j,ipcol);
					iprow = j;
				}
			}
		}
		if (iprow < 0) break;
		
		//normalize pivot row
		tmp = rsys(iprow,ipcol);
		for (j=0; j<ncols; j++)
			rsys(iprow,j) = rsys(iprow,j) / tmp;
		
		//remove from all other rows
		for (i=0; i<nrows; i++)	{
			if (i != iprow)	{
				tmp = rsys(i,ipcol) / rsys(iprow,ipcol);
				for (j=0; j<ncols; j++)
					rsys(i,j) -= tmp*rsys(iprow,j);
			}
		}
	
/*printf("\niprow=%d ipcol=%d\n",iprow,ipcol);
printf("---------------------------------------------------\n");
for (i=0; i<nrows; i++)	{
for (j=0; j<ncols; j++)
	printf(" %6.2lf",rsys(i,j));
printf("\n");
}*/
		
		
	}
//printf("\n\n******** END PHASE I ********\n");	
	//--- perform phase II simplex ---
	for (k=0; k<10*ncols; k++)	{
		ipcol = -1;
		iprow = -1;
		tmp = 0.0;
	
		//pivot column is most negative cost function
		for (j=0; j<(2*nspecies+nactive); j++)	{
			if (rsys(nrows-2,j) < tmp)	{
				tmp = rsys(nrows-2,j);
				ipcol = j;
			}
		}
		if (ipcol < 0) break;
		
		//pivot row has the smallest positive ratio
		tmp = 0.0;
		for (j=0; j<(nrows-2); j++)	{
			if ( (rsys(j,ipcol)>0.0) && (rsys(j,ncols-1)>=0.0) )	{
				if ( (tmp > rsys(j,ncols-1)/rsys(j,ipcol)) || (iprow<0) )	{
					tmp = rsys(j,ncols-1) / rsys(j,ipcol);
					iprow = j;
				}
			}
		}
		if (iprow < 0) break;
		
		//normalize pivot row
		tmp = rsys(iprow,ipcol);
		for (j=0; j<ncols; j++)
			rsys(iprow,j) = rsys(iprow,j) / tmp;
		
		//remove from all other rows
		for (i=0; i<nrows; i++)	{
			if (i != iprow)	{
				tmp = rsys(i,ipcol) / rsys(iprow,ipcol);
				for (j=0; j<ncols; j++)
					rsys(i,j) -= tmp*rsys(iprow,j);
			}
		}

/*printf("\niprow=%d ipcol=%d\n",iprow,ipcol);
printf("---------------------------------------------------\n");
for (i=0; i<nrows; i++)	{
for (j=0; j<ncols; j++)
	printf(" %6.2lf",rsys(i,j));
printf("\n");
}*/
	}
	
	//--- save basis species ---
	for (j=0; j<nspecies; j++)	{
		iflag = 0;
		if (rsys(nrows-2,j)==0.0)	{
			k=0;
			for (i=0; i<nelements; i++)	{
				if ((isactive(i)>0)&&(iflag==0))	{
					if (rsys(k,j) > 0.0)	{
						ibasis(i) = j;
						iflag = 1;
					}
					k++;
				}
			}
		}
	}
	
/*printf(" ibasis = ");
for (i=0; i<nelements; i++)
	printf("  %d ",ibasis(i));
printf("\n");*/
	
	return 0;
}
/* ------------------------------------------------------------------------------- */
