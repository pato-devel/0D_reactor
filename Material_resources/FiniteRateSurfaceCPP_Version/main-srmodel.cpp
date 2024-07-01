#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "FRMArray.hxx"
#include "FiniteRateSurface.hxx"
#include "lumatrix.hxx"

using namespace std;


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                  LOCAL PREPROCESSOR DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


/* buffer sizing and tokens */
#define INPUT_FILE_BUFFER_SIZE 1024

#define INPUT_STRTOK_TOKENS " \t',"

/* constants used by the solver */
#define PI 3.1415926536

#define R 8.31451

#define AV 6.0221367e+23

#define Pref 100000.0

#define H 6.62618e-34

/* C++ version main release number */
#define CPPRELEASE 1


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




double run_time_integration(FiniteRateSurfaceCell  &cell, double &volume, double &dt, int &nsteps, 
				    double &temperature, double &pressure, double &final_pres,
				    int &iconv, int &iout, int &iscan, int &itorder, int &igastype, ofstream &ofile);



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                       MAIN FUNCTION DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */

int main(int argc, char **argv)
{
	char                   buffer[INPUT_FILE_BUFFER_SIZE], *token;
	string                 fname, pname, bname, *snames=NULL;
	size_t                 spos;
	double                 tw, p, dt, volume, new_volume, smoltot, tmp, temperature, resid, resL1, resL2,
	                       alpha, beta, gamma, time, dti, final_pres, v, scantw, scanp, sblow, pblow, recess;
	Array1D<double>        *concentration, *srhs, *delta_conc, *qssconc, *equilconc;
	Array2D<double>        *sjacobian;
	int                    i, k, iconv, iout, igastype, iscan, nsteps, itorder, idiagnostic, ieqtype, icycle, 
	                       ierr, istep, ngas, iscanextra, iscanflag;
	ifstream               inputfile, ifile;
	ofstream               ofile;
	Array1D<double>        smolefr, sscanmole, dxold;
	FiniteRateSurfaceModel frm;
	FiniteRateSurfaceCell  cell, qss, equil;


	/*
	 * --- Write title and version info to screen ---
	 */
	cout << "SRModel\n";
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"   version:     %1d.%1d.%1d",FINITE_RATE_MODEL_MAJOR_VERSION,FINITE_RATE_MODEL_MINOR_VERSION,CPPRELEASE);
	cout << buffer << "\n";
	cout << "   written by:   Joe Marschall   (jochen.marschall@sri.com)\n";
	cout << "                 Matthew MacLean (maclean@cubrc.org)\n";
	cout << "\n";
	
	
	/*
	 * --- Read in input file ---
	 */
	inputfile.open("input.inp");
	if ( inputfile.bad() || inputfile.fail() )	{
		cout << "'ERROR: File input.inp could not be read.\n";
		return 1;
	}
	
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	// fname
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off file name token (char *)
	fname.assign(token);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	// iconv, iout, pname
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //iconv
	sscanf(token,"%d",&iconv);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //iout
	sscanf(token,"%d",&iout);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //pname
	pname.assign(token);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	// tw, p, igastype, iscan
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	string_replace_char(buffer, 'd', 'e');
	string_replace_char(buffer, 'D', 'e');
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //tw
	sscanf(token,"%lg",&tw);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //p
	sscanf(token,"%lg",&p);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //igastype
	sscanf(token,"%d",&igastype);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //iscan
	sscanf(token,"%d",&iscan);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	// dt, nsteps, itorder, idiagnostic, ieqtype
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	string_replace_char(buffer, 'd', 'e');
	string_replace_char(buffer, 'D', 'e');
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //dt
	sscanf(token,"%lg",&dt);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //nsteps
	sscanf(token,"%d",&nsteps);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //itorder
	sscanf(token,"%d",&itorder);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //idiagnostic
	sscanf(token,"%d",&idiagnostic);
	token = strtok(NULL, INPUT_STRTOK_TOKENS); //ieqtype
	sscanf(token,"%d",&ieqtype);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	// ngas
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	token = strtok(buffer, INPUT_STRTOK_TOKENS); //ngas
	sscanf(token,"%d",&ngas);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
	//   name, mole fraction (mole fraction lines)
	snames = new string [ngas];
	smolefr.resize(ngas);
	sscanmole.resize(ngas);
	for (i=0; i<ngas; i++)	{
		get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
		token = strtok(buffer, INPUT_STRTOK_TOKENS); //split off file name token (char *)
		snames[i].assign(token);
		token = strtok(NULL, INPUT_STRTOK_TOKENS);   //mole fraction
		string_replace_char(token, 'd', 'e');
		string_replace_char(token, 'D', 'e');
		sscanf(token,"%lg",&(smolefr(i)));
	}
	//
	get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);


	/*
	 * --- Limit flags ---
	 */
	
	// --- limit input flags -----
      if ((itorder != 2)&&(itorder != 0)) itorder = 1;

      if ((idiagnostic > 3)||(idiagnostic < 1)) idiagnostic=0;

	if (igastype < 0) igastype = 0;
	if (igastype > 2) igastype = 2;

      if ((ieqtype >= 3) || (ieqtype < 0)) ieqtype = igastype;

	// --- set scanning mode (turns off many flags) ---
	if (iscan > 0) {
		iscan = 1;        //start with scan #1
		iscanextra = 0;   //turn off extra reaction output
		if (idiagnostic > 0) iscanextra = 1;
		idiagnostic = 0;  //turn off all diagnostic output
		iout = -1;        //turn off all file output
		iconv = -1;       //turn off convergence output
		iscanflag = 1;
	}


	/*
	 * --- Read in input files and collect data ---
	 */
	
	ifile.open(fname.c_str());
	if ( ifile.bad() || ifile.fail() )	{
		cerr << "ERROR: " << fname << " could not be opened.\n";
		return 1;
	}
	ifile >> frm;
	if ( ifile.bad() || ifile.fail() )	{
		cerr << "ERROR: Call to read_finite_rate_surf_file failed with error code " << 1 << ". Execution halted.\n";
		return 1;
	}
	ifile.close();
	if ( ngas != frm.get_num_gas() )	{
		cerr << "ERROR: Number of species in input file must match chemistry database!\n";
		return 1;
	}


	/*
	 * --- Read pyrolysis blowing data if necessary ---
	 */
	if (frm.get_nblwflag())	{
		//look for blowing file locally first
		bname = fname;
		spos = bname.find_last_of("/");
		if (spos < bname.length()) bname.erase(0,spos+1);
		spos = bname.find_last_of(".");
		if (spos < bname.length()) bname.erase(spos);
		bname.append(".blw");
		ifile.open(bname.c_str());
		//try looking in model file directory
		if ( ifile.bad() || ifile.fail() )	{
			ifile.clear();
			bname = fname;
			spos = bname.find_last_of(".");
			if (spos < bname.length()) bname.erase(spos);
			bname.append(".blw");
		}
		if ( ifile.good() )	ierr = frm.readBlowingFile(ifile);
		if ( ifile.bad() || ifile.fail() || ierr )
			cout << "'WARNING: Blowing data file could not be read. Pyrolysis will be de-activated.\n";
		ifile.close();
	}


	/*
	 * --- Read Lewis data files ---
	 */
	
	//provide gas phase Lewis data
	ifile.open("lewis.thermo");
	if ( ifile.good() ) ierr = frm.readLewisGasFile(ifile);
	if ( ifile.bad() || ifile.fail() || ierr )	{
		cerr << "ERROR: Call to readLewisGasFile failed.\n";
		return 1;
	}
	ifile.close();

	//provide bulk phase Lewis data
	ifile.open("lewis.bulk");
	if ( ifile.good() ) ierr = frm.readLewisBulkFile(ifile);
	if ( ifile.bad() || ifile.fail() || ierr )	{
		cerr << "ERROR: Call to readLewisBulkFile failed.\n";
		return 1;
	}
	ifile.close();
	
	
	/*
	 * --- Write diagnostic file back out ---
	 */
	if (idiagnostic >= 3)	{
		ofile.open("model.out");
		ofile << frm;
		ofile.close();
	}


	/*
	 * --- Perform sanity check on the data ---
	 */
	ierr = frm.sanityCheck();
	if ((ierr>0)&&(ierr<100)) cout << "WARNING #" << ierr << ": " << frm.get_error_string() << "\n";
	if (ierr>=100)	{
		cout << "ERROR #" << ierr << ": " << frm.get_error_string() << "\n";
		return 1;
	}

	/*
	 * --- Open output stream ---
	 */
	if (iout  < 0) iout = 0;
      if (iconv < 0) iout = 0;
      if (iout  != 0) ofile.open(pname.c_str());


	/*
	 * --- Initialize single surface cell ---
	 */
	
	/* --- initialize cell --- */
	cell.initialize_cell(&frm);
	
	/* --- initialize volume & temperature --- */
	volume      = 1.0;
	new_volume  = 1.0;
	temperature = tw;

	/* --- normalize mole fraction values --- */
	smoltot = 0.0;
	for (i=0; i<ngas; i++)	{
		if (smolefr(i) < 1.0e-15) smolefr(i) = 1.0e-15;
		smoltot += smolefr(i);
	}
	for (i=0; i<ngas; i++)
		smolefr(i) = smolefr(i) / smoltot;
	
	/* --- set gas phase concentrations --- */
	concentration = cell.get_concentration_ptr();
	for (i=0; i<ngas; i++)
		(*concentration)(i) = p * smolefr(i) / (R*tw);
		
	/* --- initialize to QSS if needed --- */
	if (frm.get_initsurf()) cell.compute_qss_system(tw);

	/* --- output initial state --- */
	cout << "\n";
	cout << "---Initial Conditions-----------------------------" << "\n";
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       Wall Tmpr   = %13.4le K",tw);
	cout << buffer << "\n";
	snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       Pressure    = %13.4le Pa",p);
	cout << buffer << "\n";
	for (i=0; i<cell.get_num_gas(); i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le mol/m^3",frm.get_species_name(i),(*concentration)(i));
		cout << buffer << "\n";
	}
	for (i=cell.get_num_gas(); i<cell.get_num_nonbulk(); i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le mol/m^2",frm.get_species_name(i),(*concentration)(i));
		cout << buffer << "\n";
	}
	for (i=cell.get_num_nonbulk(); i<cell.get_num_species(); i++)	{
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le []",frm.get_species_name(i),(*concentration)(i));
		cout << buffer << "\n";
	}
	
	/* --- initialize equilibrium cell & save initial concentrations --- */
	equil.initialize_cell(&frm);
	equilconc = equil.get_concentration_ptr();
	for (i=0; i<cell.get_num_gas(); i++)
		(*equilconc)(i) = (*concentration)(i);


	
	
	/*
	 * --- Write header info and track scan mode ---
	 */
	cout << "\n";
	if (iscan <= 0)	{
		cout << "--------------------------------------------------" << "\n";
		cout << "iteration         time       L1-norm       L2-norm" << "\n";
	}
	
	if (iout >= 0)	{
		ofile << "iteration          time";
		for (i=0; i<cell.get_num_species(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE," %15s",frm.get_species_name(i));
			ofile << buffer; 
		}
		ofile << "\n";
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE, "%9d %13.4lg", 0, 0.0);
		ofile << buffer; 
		for (i=0; i<cell.get_num_species(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE, " %15.6lg", (*concentration)(i));
			ofile << buffer; 
		}
		ofile << "\n";
	}
	
	
	/*
	 * --- Main routine enclosed for scanning ---
	 */
	if (iscan <= 0) /* single mode */	{
	
		/* time integration */
		run_time_integration(cell, volume, dt, nsteps, temperature, p, final_pres, iconv, iout, iscan, itorder, igastype, ofile);
	
		/* final condition output */
		cout << "\n";
		cout << "---Final Conditions-------------------------------" << "\n";
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       Wall Tmpr   = %13.4le K",temperature);
		cout << buffer << "\n";
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       Pressure    = %13.4le Pa",final_pres);
		cout << buffer << "\n";
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       Rel. volume = %13.4le",volume);
		cout << buffer << "\n";
		for (i=0; i<cell.get_num_gas(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le mol/m^3",frm.get_species_name(i),(*concentration)(i));
			cout << buffer << "\n";
		}
		for (i=cell.get_num_gas(); i<cell.get_num_nonbulk(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le mol/m^2",frm.get_species_name(i),(*concentration)(i));
			cout << buffer << "\n";
		}
		for (i=cell.get_num_nonbulk(); i<cell.get_num_species(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le []",frm.get_species_name(i),(*concentration)(i));
			cout << buffer << "\n";
		}
	
	
	} else /* scanning mode on */ {

		while (iscan < 4)	{

			// - reset time integration terms
			dxold = 0.0;
			volume = 1.0;
			new_volume = volume;
			
			// - read a new line of settings
			if (iscanflag == 1)	{
			
				get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
				get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
				get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
				iscan++;
				iscanflag = 0;
				if (iscan == 2)	{
					cout << "---Pressure/Temperature Scan----------------------\n";
					cout << "            Pi             Tw         Volume             Pf          resL2";
				}
				if (iscan == 3)	{
					cout << "---Mole Fraction Scan----------------------\n";
					for (k=0; k<cell.get_num_gas(); k++)	{
						snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"         mfi_%s",frm.get_species_name(k));
						cout << buffer; 
					}
					cout << "        Volume            Pf         resL2";	
				}
				if (iscan <= 3)	{
					for (k=0; k<cell.get_num_species(); k++)	{
						snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"            x%s",frm.get_species_name(k));
						cout << buffer; 
					}
					for (k=0; k<cell.get_num_gas(); k++)	{
						snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"            g%s",frm.get_species_name(k));
						cout << buffer; 
					}
					if (iscanextra > 0)	{
						for (k=0; k<cell.get_num_species(); k++)	{
							for (i=0; i<frm.get_num_reactions(); i++)	{
								snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"            x%s_%d",frm.get_species_name(k),i+1);
								cout << buffer; 
							}
						}
					}
					cout << "\n";
				}		
			}
			
			// - read a new row of data
			if (iscan == 2)	{
				get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
				string_replace_char(buffer, 'd', 'e');
				string_replace_char(buffer, 'D', 'e');
				token = strtok(buffer, INPUT_STRTOK_TOKENS); //p
				if(sscanf(token,"%lg",&scanp) < 1) iscanflag = 1;
				token = strtok(NULL, INPUT_STRTOK_TOKENS); //tw
				if(sscanf(token,"%lg",&scantw) < 1) iscanflag = 1;
				
				sscanmole = smolefr;
				
				if ((scantw<0.0)||(scanp<0.0)) iscanflag = 1;
			
			} else if (iscan == 3)	{
			
				get_line_of_text(inputfile, buffer, INPUT_FILE_BUFFER_SIZE);
				string_replace_char(buffer, 'd', 'e');
				string_replace_char(buffer, 'D', 'e');
				token = strtok(buffer, INPUT_STRTOK_TOKENS); //xi(0)
				if(sscanf(token,"%lg",&(sscanmole(0))) < 1) iscanflag = 1;
				if(sscanmole(0) < 0.0) iscanflag = 1;
				for (i=1; i<cell.get_num_gas(); i++)	{
					token = strtok(NULL, INPUT_STRTOK_TOKENS); //xi
					if(sscanf(token,"%lg",&(sscanmole(i))) < 1) iscanflag = 1;
					if(sscanmole(i) < 0.0) iscanflag = 1;
				}
				scantw = tw;
				scanp  = p;
			}

			// - run analysis unless we're scanning for the next round
			if ( (iscanflag == 0) && (iscan<4) )	{
				//reset concentrations
				smoltot = 0.0;
				for (i=0; i<ngas; i++)	{
					if (sscanmole(i) < 1.0e-15) sscanmole(i) = 1.0e-15;
					smoltot += sscanmole(i);
				}
				for (i=0; i<ngas; i++)
					sscanmole(i) = sscanmole(i) / smoltot;
	
				// set gas phase concentrations
				concentration = cell.get_concentration_ptr();
				for (i=0; i<ngas; i++)
					(*concentration)(i) = scanp * sscanmole(i) / (R*scantw);
		
				// reset surface concentrations (steal from equilconc)
				for (i=cell.get_num_gas(); i<cell.get_num_nonbulk(); i++)
					(*concentration)(i) = (*equilconc)(i);
				
				// initialize to QSS if needed
				if (frm.get_initsurf()) cell.compute_qss_system(scantw);

				// run time integration
				resL2 = run_time_integration(cell, volume, dt, nsteps, scantw, scanp, final_pres, iconv, iout, iscan, itorder, igastype, ofile);

				//write output to screen
				if (iscan == 2)	{
					snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"%14.6le %14.6le %14.6le %14.6le %14.6le", scanp, scantw, volume, final_pres, resL2);
					cout << buffer;
				} else if (iscan == 3)	{
					for (i=0; i<cell.get_num_gas(); i++)	{
						snprintf(buffer, INPUT_FILE_BUFFER_SIZE," %14.6le", sscanmole(i));
						cout << buffer;
					}
					snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"%14.6le %14.6le %14.6le", volume, final_pres, resL2);
					cout << buffer;
				}
				cell.compute_frm_system(scantw);
				concentration = cell.get_concentration_ptr();
				srhs = cell.get_rhs_ptr();
				for (i=0; i<cell.get_num_species(); i++)	{
					snprintf(buffer, INPUT_FILE_BUFFER_SIZE," %14.6le",(*concentration)(i));
					cout << buffer;
				}
				for (i=0; i<cell.get_num_gas(); i++)	{
					v = frm.get_mean_thermal_speed(i, scantw);
					snprintf(buffer, INPUT_FILE_BUFFER_SIZE," %14.6le",-(*srhs)(i)/(*concentration)(i)/v);
					cout << buffer;
				}
				if (iscanextra > 0) frm.print_diagnostics(cout, cell, scantw, 1);
				cout << "\n";
			}
		}
	}

	/* --- close output file --- */
	if (iout != 0) ofile.close();
	
	
	/*
	 * --- Post-processed quantities ---
	 */
	
	if (iscan <= 0)	{
	
		/* --- calculate loss efficiencies (recombination probability) --- */
		cout << "\n" << "---Loss Efficiencies------------------------------\n";
		cell.compute_frm_system(temperature);
		concentration = cell.get_concentration_ptr();
		srhs = cell.get_rhs_ptr();
		for (i=0; i<cell.get_num_gas(); i++)	{
			v = frm.get_mean_thermal_speed(i, temperature);
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       g%-10s = %13.4le",frm.get_species_name(i),-((*srhs)(i))/((*concentration)(i))/v);
			cout << buffer << "\n";
		}
		
		/* --- calculate blowing and recession rates --- */
		frm.compute_recession_rates(cell, sblow, pblow, recess);
		cout << "\n" << "---Blowing & Recession Rate-----------------------\n";
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"blowing rate   %13.4e kg/m^2-s",sblow+pblow);
		cout << buffer << "\n";
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"recession rate %13.4e m/s",recess);
		cout << buffer << "\n";
	}


	/*
	 * --- Perform QSS (Quasi-Steady State) Analysis ---
	 */
	if (idiagnostic >= 1)	{
	
		/* --- write header --- */
		cout << "\n";
		cout << "---QSS Surface Coverage Distribution--------------\n";

		/* --- initialize --- */
		qss.initialize_cell(&frm);

		/* --- set gas phase concentrations --- */
		qssconc = qss.get_concentration_ptr();
		for (i=0; i<cell.get_num_gas(); i++)
			(*qssconc)(i) = (*concentration)(i);
		
		/* --- QSS Solver --- */
		resid = qss.compute_qss_system(temperature);
	
		/* --- Write to screen --- */
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"      Residual = %13.4le",sqrt(resid));
		cout << buffer << "\n";
		cout << "      Concentrations\n";
		for (i=cell.get_num_gas(); i<cell.get_num_nonbulk(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%s    = %13.4le mol/m^2",frm.get_species_name(i),(*qssconc)(i));
			cout << buffer << "\n";
		}
	}


	/*
	 * --- Perform Equilibrium Computation ---
	 */
	if (idiagnostic >= 1)	{

		/* --- write header --- */
		cout << "\n";
		cout << "---Gibbs Free-energy Minimization-----------------\n";
		frm.print_etable(cout);
		
		/* --- Equilibrium Solver --- */
		equil.compute_equilibrium_system(temperature, igastype, &cout);

		/* --- write output --- */
		equilconc = equil.get_concentration_ptr();
		cout << "      Concentrations\n";
		for (i=0; i<cell.get_num_gas(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le mol/m^3",frm.get_species_name(i),(*equilconc)(i));
			cout << buffer << "\n";
		}
		for (i=cell.get_num_gas(); i<cell.get_num_nonbulk(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le mol/m^2",frm.get_species_name(i),(*equilconc)(i));
			cout << buffer << "\n";
		}
		for (i=cell.get_num_nonbulk(); i<cell.get_num_species(); i++)	{
			snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"       x%-10s = %13.4le []",frm.get_species_name(i),(*equilconc)(i));
			cout << buffer << "\n";
		}

	}
	
	/*
	 * --- Write extra diagnostics ---
	 */
	if (idiagnostic >= 2)	{
	
		/* --- write header --- */
		cout << "\n";
		cout << "---Additional Diagnostics-----------------------------------------------------------\n";
		
		/* --- write diagnostics --- */
		frm.print_diagnostics(cout, cell, temperature);
		
	
	}


	/*
	 * --- Clean-up Memory & Close Input ---
	 */
	inputfile.close();
	
	
	if (snames) delete [] snames;
	
	return 0;
}





double run_time_integration(FiniteRateSurfaceCell  &cell, double &volume, double &dt, int &nsteps, 
				    double &temperature, double &pressure, double &final_pres,
				    int &iconv, int &iout, int &iscan, int &itorder, int &igastype, ofstream &ofile)
{
	int             istep, i, k;
	double          resL1, resL2, time, tmp, new_volume, alpha, beta, gamma, dti;
	Array1D<double> *concentration, *srhs, *delta_conc, dxold;
	Array2D<double> *sjacobian;
	char            buffer[INPUT_FILE_BUFFER_SIZE];

	/*
	 * --- Set time integration flags ---
	 */
	
	if (itorder == 2)	{
		alpha = 1.5;
		beta  = 0.5;
		gamma = 1.0;
	} else if (itorder == 1)	{
		alpha = 1.0;
		beta  = 0.0;
		gamma = 1.0;
	} else {
		alpha = 1.0;
		beta  = 0.0;
		gamma = 0.0;
	}
	
	dti  = alpha / dt;
	beta = beta / dt;
	
	dxold.resize(cell.get_num_eqns());
	dxold.initialize(0.0);


	/* --- Set pointers --- */
	concentration = cell.get_concentration_ptr();
	srhs          = cell.get_rhs_ptr();
	delta_conc    = cell.get_delta_conc_ptr();
	sjacobian     = cell.get_jacobian_ptr();
	
	/* --- Integrate in time --- */
	time = 0.0;
	
	for (istep=0; istep < nsteps; istep++)	{
	
		// - compute implicit system -
		cell.compute_frm_system(temperature);
		
		// - freeze gas phase at user request -
		if (igastype == 1)	{
			for (k=0; k<cell.get_num_gas(); k++)	{
				for (i=0; i<cell.get_num_eqns(); i++)
					(*sjacobian)(k,i) = 0.0;
				(*sjacobian)(k,k) = 1.0;
				(*srhs)(k) = 0.0;
			}
		}
		
		// - freeze bulk phase (always) -
		for (k=cell.get_num_nonbulk(); k<cell.get_num_species(); k++)	{
			for (i=0; i<cell.get_num_eqns(); i++)
				(*sjacobian)(k,i) = 0.0;
			(*sjacobian)(k,k) = 1.0;
			(*srhs)(k) = 0.0;
		}
		
		// - add time step term and negate -
		for (k=0; k<cell.get_num_eqns(); k++)	{
			for (i=0; i<cell.get_num_eqns(); i++)	{
				if (i == k) (*sjacobian)(i,k) = dti - gamma*((*sjacobian)(i,k));
				else        (*sjacobian)(i,k) =     - gamma*((*sjacobian)(i,k));
			}
		}

		// - add 2nd order time explicit term -
		for (k=0; k<cell.get_num_eqns(); k++)
			(*srhs)(k) += beta * dxold(k);

		// - solve system -
		lusolve(cell.get_jacobian(), cell.get_rhs(), cell.get_rhs());

		// - compute new volume for constant pressure -
		if (igastype == 2)	{
			new_volume = 0.0;
			for (k=0; k<cell.get_num_gas(); k++)	{
				tmp = (*srhs)(k) + volume * (*concentration)(k);
				if (tmp < 1.0e-20) tmp = 1.0e-20;
				new_volume += tmp;
			}
			new_volume = new_volume * (R*temperature/pressure);
		} else
			new_volume = volume;
		
		// - update variables (accounting for volume change) -
		for (k=0; k<cell.get_num_gas(); k++)	{
			tmp = volume * (*concentration)(k);
			(*concentration)(k) = ( tmp + (*srhs)(k) ) / new_volume;
			if ((*concentration)(k) < 1.0e-20) (*concentration)(k) = 1.0e-20;
		}
		for (k=cell.get_num_gas(); k<cell.get_num_species(); k++)	{
			(*concentration)(k) += (*srhs)(k);
			if ((*concentration)(k) < 1.0e-20) (*concentration)(k) = 1.0e-20;
		}
		temperature = temperature + (*srhs)(cell.get_num_species());
		volume = new_volume;
		
		// - increment time index -
		time += dt;
		
		// - save solution for next iteration (2nd order time)
		for (k=0; k<cell.get_num_eqns(); k++)
			dxold(k) = (*srhs)(k);
		
		// - compute residual -
		resL1 = 0.0;
		resL2 = 0.0;
		for (k=0; k<cell.get_num_species(); k++)	{
			tmp = abs((*srhs)(k));
			if (tmp > resL1) resL1 = tmp;
			resL2 += tmp * tmp;
		}
		resL2 = sqrt(resL2);
		
		// - write output to screen & file -
		snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"%8d %13.4le %13.4le %13.4le",istep,time,resL1,resL2);
		if (iconv > 0)
			if ( ((istep+1)%iconv == 0) || (istep==(nsteps-1)) )	cout << buffer << "\n";
		else
			if ( (istep==(nsteps-1)) && (iscan <= 0) )	cout << buffer << "\n";
		
		if (iout > 0)	{
			if ((istep+1)%iout == 0)	{
				snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"%8d %13.4le",istep,time);
				ofile << buffer;
				for (k=0; k<cell.get_num_species(); k++)	{
					snprintf(buffer, INPUT_FILE_BUFFER_SIZE,"%15.6le",(*concentration)(k));
					ofile << buffer;
				}
				ofile << "\n";
			}
		}
	}

	/* --- compute final gas pressure --- */
	final_pres = 0.0;
	for (k=0; k<cell.get_num_gas(); k++)
		final_pres += (*concentration)(k) * R * temperature;

	return resL2;
}
