/*   
 * FiniteRateSurface.hxx
 * ---------------------
 *    file created:   January   16, 2013
 *    last modified:  August    24, 2013
 *    author:         Matthew MacLean
 *                    maclean@cubrc.org
 *    purpose:        define classes for Finite Rate Surface Chemistry
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

#if !defined(_FINITE_RATE_SURFACE_DOT_H)
#define _FINITE_RATE_SURFACE_DOT_H


/* global includes */
#include <string>
#include <iostream>
#include <iomanip>

/* local includes */
#include "FRMArray.hxx"


/* namespace */
using namespace std;



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                     PREPROCESSOR DEFINITIONS                           ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


/* This is the major revision of the original FRM
   Incremented based on compatibility with Joe Marschall's original codes */
#define FINITE_RATE_MODEL_MAJOR_VERSION  3

/* This is the minor revision of the original FRM 
   Incremented when a bug is found in the finite rate model formulation */
#define FINITE_RATE_MODEL_MINOR_VERSION  8

/* This is a local release to indicate a change in this C++ source code */
#define FINITE_RATE_MODEL_LOCAL_RELEASE  1

/* Hardwired Array sizing parameters */
#define ISMODEL_MAX_REACT 3   //maximum number of reactant species
#define ISMODEL_MAX_PROD  3   //maximum number of product species
#define NLEWIS_FITS 6         //maximum number of Lewis fits
#define NLEWIS_COEFFS 10      //number of Lewis coefficients in each range


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###              PUBLIC STRUCTURE (CLASS) DECLARATIONS                     ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


/*class declarations */
class FiniteRateSurfaceModel;
class FiniteRateSurfaceCell;



/*  ----------------- <CLASS> FiniteRateSurfaceModel <CLASS> ---------------------
 *
 *  description:  encapsulation of all FRM chemistry database for a given model
 *
 *  ------------------------------------------------------------------------------
 */
class FiniteRateSurfaceModel
{
/* PUBLIC Fields */
public:

	/* ----- PUBLIC CONSTRUCTOR ----- */
	FiniteRateSurfaceModel();
	
	/* ----- PUBLIC DESTRUCTOR ----- */
	
	/* arguments: none
       *
       * purpose:   default destructor
       */
	~FiniteRateSurfaceModel();
	
	/* ----- PUBLIC ACCESSOR METHODS ----- */
	
	/* arguments: none
       *
       * purpose:   trivial file accessors
       */
	int get_num_species()	{ return nspecies;	}
	int get_num_gas()		{ return igend-igstart;	}
	int get_num_surf()	{ return isend-isstart;	}
	int get_num_bulk()	{ return ibend-ibstart;	}
	int get_num_phases()	{ return nphases;	}
	int get_num_sitesets()	{ return nssites; }
	int get_num_reactions()	{ return isrt; }
	int get_nblwflag()	{ return nblwflag; }
	int get_initsurf()	{ return initsurf;	}
	char *get_error_string(){ return cfrm_error_string;	}

	/* arguments: i = species #
       *
       * purpose:   retrieve species name string
       */
	const char *get_species_name(int i)
	{	if ((i>=0)&&(i<nspecies)) return species(i).c_str(); else return "";	}
	
	/* arguments: i = species #
       *
       * purpose:   retrieve species molecular weight (returns in kg/gmol)
       */
	double get_species_molewgt(int i)
	{	if ((i>=0)&&(i<nspecies)) return wmass(i); else return 0.0;	}

	/* arguments: look_for     = string name to look for
	 *            strip_braces = make the comparison by removing the parenthesis designation from surface and bulk
	 *            search_gas   = include gas species in search
	 *            search_surf  = include surface species in search
	 *            search_bulk  = include bulk species in search
       *
       * purpose:   find the position of a species with a given name
	 *
	 * return:    position of species (0 - nspecies) or -1 if not present
       */
	int get_species_name_pos(const string &look_for, bool strip_braces, bool search_gas, bool search_surf, bool search_bulk);

	/* arguments: i    = species #
	 *            tmpr = temperature
       *
       * purpose:   compute mean thermal speed
       */
	double get_mean_thermal_speed(int i, double tmpr);

	/* arguments: cell = cell reference to modify
	 *
	 * return:    nothing
	 *
	 * purpose:   initialize cell contents and associate the model
	 */ 
	void initialize_cell(FiniteRateSurfaceCell& cell);

	/* arguments: cell = cell reference to modify
	 *            tmpr = temperature
	 *
	 * return:    nothing
	 *
	 * purpose:   compute finite rate system at tmpr, pres and set into cell 
	 */ 
	void compute_frm_system(FiniteRateSurfaceCell& cell, double tmpr);

	/* arguments: cell = cell reference to modify
	 *            tmpr = temperature
	 *
	 * return:    algorithm residual
	 *
	 * purpose:   compute Quasi-Steady state surface distribution at tmpr, pres and set into cell 
	 */ 
	double compute_qss_system(FiniteRateSurfaceCell& cell, double tmpr);

	/* arguments: cell       = cell reference to modify
	 *            igastype   = type of equilibrium condition ()
	 *            stream_ptr = pointer to an output stream to write diagnostics to (NULL to ignore)
	 *
	 * return:    algorithm residual
	 *
	 * purpose:   compute equilibrium state surface distribution at tmpr, pres and set into cell 
	 */ 
	double compute_equil_system(FiniteRateSurfaceCell& cell, int igastype, ostream *stream_ptr = NULL);
	
	/* arguments: stream = stream to write to
	 *            cell = cell reference to modify
	 *            tmpr = temperature
	 *            iscan = write diagnostics in "iscan" mode
	 *
	 * return:    nothing
	 *
	 * purpose:   print diagnostic rates to stream 
	 */ 
	void print_diagnostics(ostream &stream, FiniteRateSurfaceCell& cell, double tmpr, int iscan=0);
	
	/* arguments: stream = stream to write to
	 *
	 * return:    nothing
	 *
	 * purpose:   print the element table out
	 */ 
	void print_etable(ostream &stream);
	
	/* arguments: cell   = cell reference to modify
	 *            sblow  = surface blowing
	 *            pblow  = pyrolysis blowing
	 *            recess = recession rate
	 *
	 * return:    nothing
	 *
	 * purpose:   set surface blowing, pyrolysis blowing, and recession rate
	 */ 
	void compute_recession_rates(FiniteRateSurfaceCell& cell, double &sblow, double &pblow, double &recess);


	/* ----- PUBLIC MUTATOR METHODS ----- */
	
	/* arguments: stream = input stream of blowing file
       *
	 * return:    0 = no error (success), else error
	 *
       * purpose:   read in blowing file data
       */
	int readBlowingFile(istream &stream);
	
	/* arguments: stream = open istream for the lewis.thermo file
	 *
	 * return:    0 = no error (success), else error
	 *
	 * purpose:   read Lewis gas data
       */
	int readLewisGasFile(istream &stream);
	
	/* arguments: stream = open istream for the lewis.bulk file
	 *
	 * return:    0 = no error (success), else error
	 *
	 * purpose:   read Lewis bulk data
       */
	int readLewisBulkFile(istream &stream);

	/* arguments: none
       *
	 * return:    0 = no error (success), else error
	 *
       * purpose:   perform checks on the model data
       */
	int sanityCheck();
	
	
	/* ----- PUBLIC OPERATORS ----- */

	/* arguments: stream = stream to read from
	 *            frm    = reference to the object to put the data to
	 *
	 * return:    the same stream reference
	 *
	 * purpose:   read a complete object from a filestream
	 *
	 */ 
	friend istream &operator>>(istream &stream, FiniteRateSurfaceModel &frm);

	/* arguments: stream = stream to write to
	 *            frm    = reference to the object to write out
	 *
	 * return:    the same stream reference
	 *
	 * purpose:   output a complete object to a filestream
	 *
	 */ 
	friend ostream &operator<<(ostream &stream, FiniteRateSurfaceModel &frm);


	/* ----- PUBLIC CLASS DATA ----- */



/* Protected Fields */
protected:
	
	/* ----- PROTECTED CLASS METHODS ----- */

	/* arguments: gibbs   = array to store Gibbs energies (returned)
	 *            dgdTg   = array to store Gibbs gas tmpr derivatives (returned)
	 *            dgdTw   = array to store Gibbs wall tmpr derivatives (returned)
	 *            Tw      = surface temperature
	 *            Tg      = gas temperature
	 *            frate   = existing forward rate
	 *            dfdTg   = existing forward rate derivative w.r.t. gas temperature
	 *            dfdTw   = existing forward rate derivative w.r.t. wall temperature
	 *
	 * return:    nothing
	 *
       * purpose:   compute Gibbs energies of all species
       */
	void compute_gibbs_energies(Array1D<double> &gibbs, Array1D<double> &dgdTg, Array1D<double> &dgdTw,
					    double Tw, double Tg, Array1D<double> const &frate, Array1D<double> const &dfdTg, Array1D<double> const &dfdTw);

	/* arguments: EqC     = array to store EqC energies (returned)
	 *            dEqCdTg = array to store EqC gas tmpr derivatives (returned)
	 *            dEqCdTw = array to store EqC wall tmpr derivatives (returned)
	 *            EqA     = array to store EqA energies (returned)
	 *            dEqAdTg = array to store EqA gas tmpr derivatives (returned)
	 *            dEqAdTw = array to store EqA wall tmpr derivatives (returned)
	 *            Tw      = surface temperature
	 *            Tg      = gas temperature
	 *            gibbs   = array of existing Gibbs energies
	 *            dgdTg   = array of existing Gibbs gas tmpr derivatives
	 *            dgdTw   = array of existing Gibbs wall tmpr derivatives
	 *
	 * return:    nothing
	 *
       * purpose:   compute concentration-based equilibrium constants
       */
	void compute_eqcon(Array1D<double> &EqC, Array1D<double> &dEqCdTg, Array1D<double> &dEqCdTw,
				 Array1D<double> &EqA, Array1D<double> &dEqAdTg, Array1D<double> &dEqAdTw,
				 double Tw, double Tg, Array1D<double> const &gibbs, Array1D<double> const &dgdTg, Array1D<double> const &dgdTw);
	
	/* arguments: Tw      = surface temperature
	 *            Tg      = gas temperature
	 *            wm      = molecular weight
	 *            stot    = total surface density
	 *            p       = site density exponent
	 *            itype   = type of surface reaction
	 *            A, B, E = reaction specific parameters
	 *            dkdTg   = derivative w.r.t. gas temperature (return)
	 *            dkdTw   = derivative w.r.t. wall temperature (return)
	 *
	 * return:    reaction rate
	 *
       * purpose:   compute forward reaction rate for various forms
       */
	double srate(double Tw, double Tg, double wm, double stot, double p, int itype, double A, double B, double E, double *dkdTg, double *dkdTw);
	
	/* arguments: Tw          = surface temperature
	 *            Tg          = gas temperature
	 *            wm          = molecular weight
	 *            stot        = total surface density
	 *            itype       = type of equilibrium constant
	 *            iform       = form of equilibrium constant
	 *            A,B,v,Ed,Ea = reaction specific parameters
	 *            frate       = existing forward rate
	 *            dfdTg,dfdTw = existing forward rate derivatives
	 *            dEQdTg      = derivative w.r.t. gas temperature (return)
	 *            dEQdTw      = derivative w.r.t. wall temperature (return)
	 *
	 * return:    activity equilibrium constant
	 *
       * purpose:   compute equilibrium constant for various forms
       */
	double adeq(double Tw, double Tg, double wm, int itype, int iform, double A, double B, double v,
			double Ed, double Ea, double frate, double dfdTg, double dfdTw, double *dEQdTg, double *dEQdTw);

	/* arguments: none
	 *
	 * return:    error codes (0=success)
	 *
       * purpose:   build the elemental stoichiometry table from species names
       */
	int build_etable();


	/* ----- PROTECTED CLASS DATA ----- */

	//      ---- title
        char title[200];

	//      ---- initialization flag(s)
        int initsurf;

	//      ---- phase and species information
        int          ngp,nsp,nbp,nptot,kstot,ngps,nphases,nspecies,nspasmax,nssites,
                     igstart,igend,isstart,isend,ibstart,ibend;
        Array1D<int> nsps, nbps, ithermo, nspas,
                     kphase, kthermo, kadsr, ksites; //1D arrays
        Array2D<int> nspass;                         //2D arrays

	//      ---- Thermodynamic Lewis data ----
        Array3D<double> rlewis_coeffs,rlewis_trange;  //3D arrays
        Array1D<int>    rlewis_nrange;                //1D arrays

	//       ---- Reaction information
        int isrt;
        Array2D<int>    RR,RP,RVR,RVP;      //2D arrays
        Array2D<double> vr,vp,v;            //2D arrays
        Array1D<int> isrfon,isrbon,istype,
                     iadtype,iadform;       //1D arrays

	//       ---- Phase and species information
        double          sdentot,sfrcI,vfbI,bdentot,bportot;
        Array1D<double> sden,sfrc,den,por,vfb;  //1D arrays
        Array2D<double> sumv,bmf,sdenas;        //2D arrays
        Array1D<double> wmass,Ediss,Ed,Eact;    //1D arrays

	//      ---- Equation information
        Array1D<double> Cf,eta,Ea,sitep;      //1D arrays
        Array1D<double> Cfad,etaad,vad,Edes;  //1D arrays
        Array1D<double> sdensity,sfraction;   //1D arrays

	//      ---- Pyrolysis blowing information
        int             nblwflag,nebc,nblw,isblwon;
        double          cyield,ha;
        Array1D<double> rmdotg,pyroadd;  //1D arrays
        Array2D<double> blwmf;           //2D arrays

	//      ---- Species and phase names
	  Array1D<string> species, phases; //1D string arrays

	//      ---- debugging flag ---- 
        int ifrm_debug;

	//      ---- error string ----
        char cfrm_error_string[1024];

	//      ---- element table -----
	  Array2D<double> etable;
	  int nelements;
	  Array1D<string> elemnames;

	//      ---- initialization flag -----
	  bool initialized;

};
/* --------------------<END> FiniteRateSurfaceModel <END>------------------------- */




/*  ------------------ <CLASS> FiniteRateSurfaceCell <CLASS> ---------------------
 *
 *  description:  store data for state of finite rate surface location
 *
 *  ------------------------------------------------------------------------------
 */
class FiniteRateSurfaceCell
{
/* PUBLIC Fields */
public:

	/* ----- PUBLIC CONSTRUCTOR ----- */
	
	/* arguments: none
       *
       * purpose:   default constructor
       */
	FiniteRateSurfaceCell();


	/* ----- PUBLIC DESTRUCTOR ----- */
	
	/* arguments: none
       *
       * purpose:   default destructor
       */
	~FiniteRateSurfaceCell();


	/* ----- PUBLIC ACCESSOR METHODS ----- */
	
	/* arguments: none
       *
       * purpose:   series of trivial data accessors
       */
	int get_num_species()	{ if (model) return model->get_num_species(); else return 0;	}
	int get_num_gas()		{ if (model) return model->get_num_gas(); else return 0;	}
	int get_num_surf()	{ if (model) return model->get_num_surf(); else return 0;	}
	int get_num_bulk()	{ if (model) return model->get_num_bulk(); else return 0;	}
	int get_num_nongas()	{ if (model) return (model->get_num_surf()+model->get_num_bulk()); else return 0;	}
	int get_num_nonbulk()	{ if (model) return (model->get_num_gas()+model->get_num_surf()); else return 0;	}
	int get_num_eqns()	{ if (model) return (model->get_num_species()+1); else return 0;	}
	double get_pressure()	{ return pres;	}
	double get_temperature(){ return tmpr;	}
	
	/* arguments: none
       *
       * purpose:   reference accessors for Array objects
	 *            (use carefully)
       */
	Array1D<double>& get_concentration(){	return concentration;	}
	Array1D<double>& get_rhs()		{	return srhs;	}
	Array1D<double>& get_delta_conc()	{	return delta_conc;	}
	Array2D<double>& get_jacobian()	{	return sjacobian;	}
	
	Array1D<double> *get_concentration_ptr()	{	return &concentration;	}
	Array1D<double> *get_rhs_ptr()		{	return &srhs;	}
	Array1D<double> *get_delta_conc_ptr()	{	return &delta_conc;	}
	Array2D<double> *get_jacobian_ptr()		{	return &sjacobian;	}
	
	const Array1D<double>& get_concentration_const() const {	return concentration;	}
	const Array1D<double>& get_rhs_const() const		 {	return srhs;	}
	const Array1D<double>& get_delta_conc_const() const	 {	return delta_conc;	}
	const Array2D<double>& get_jacobian_const() const	 {	return sjacobian;	}

	
	/* ----- PUBLIC MUTATOR METHODS ----- */
	
	/* arguments: frsmodel = an initialized FiniteRateSurfaceModel object to associate to
	 *
	 * return:    nothing
	 *
	 * purpose:   initialize cell contents and associate the model
	 */ 
	void initialize_cell(FiniteRateSurfaceModel *frsmodel);
	
	/* arguments: tmpr = temperature
	 *
	 * return:    nothing
	 *
	 * purpose:   compute finite rate system at tmpr, pres and set in 
	 */ 
	void compute_frm_system(double tmpr);
	
	/* arguments: tmpr       = temperature
	 *            igastype   = type of equilibrium solution
	 *            stream_ptr = pointer to output stream for diagnostic write (NULL to ignore)
	 *
	 * return:    residual of the solution
	 *
	 * purpose:   compute equilibrium composition and set concentrations
	 */ 
	double compute_equilibrium_system(double tmpr, int igastype, ostream *stream_ptr);
	
	/* arguments: tmpr = temperature
	 *
	 * return:    residual of the solution
	 *
	 * purpose:   compute QSS distribution and set surface concentrations
	 */ 
	double compute_qss_system(double tmpr);


	/* ----- PUBLIC OPERATORS ----- */



	/* ----- PUBLIC CLASS DATA ----- */	



/* Protected Fields */
protected:

	FiniteRateSurfaceModel *model;

	double tmpr, pres;
	
	Array1D<double> concentration, srhs, delta_conc;
	Array2D<double> sjacobian;

};
/* ---------------------<END> FiniteRateSurfaceCell <END>------------------------- */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC VARIABLE DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */




/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC FUNCTION DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/* 
 * name:      get_line_of_text
 *
 * arguments: stream = stream to read from
 *            buffer = C-style character array
 *            size   = maximum length of buffer
 *
 * return:    nothing
 *
 * purpose:   read a line from a stream into a buffer including any combo of newline characters (thrown away)
 *
 */ 
void get_line_of_text(istream &stream, char *buffer, int size);



/* 
 * name:      string_replace_char
 *
 * arguments: str          = string to search through
 *            search_for   = char to search for
 *            replace_with = char to replace with
 *
 * return:    nothing
 *
 * purpose:   search a NULL terminated C-string and replace one character with another
 *
 */ 
void string_replace_char(char *str, char search_for, char replace_with);



/* 
 * name:      string_is_end
 *
 * arguments: str = string to compare
 *
 * return:    0 = END, != 0 otherwise
 *
 * purpose:   see if a string is equal to END (case insensitive)
 *
 */ 
int string_is_end(char *str);



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC OPERATOR DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */




#endif   //  _FINITE_RATE_SURFACE_DOT_H
