#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#endif

#if defined(_WIN32)
#include <windows.h> // for the CreateDirectory function
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CoolPropTools.h"
#include "CoolProp.h"
#include "CPState.h"
#include "TTSE.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
#include "time.h"
#include <sys/types.h>
#include <sys/stat.h>	
#include <errno.h>
#include <cerrno>
#include <float.h>

// An ugly hack to disable the timing function on PPC since the sysClkRate function not found
#if defined(__powerpc__)
#define CLOCKS_PER_SEC 1000
#endif

static bool transport_properties = true;

double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

TTSESinglePhaseTableClass::TTSESinglePhaseTableClass(){
	this->enable_writing_tables_to_files = true;
	SatL = NULL; 
	SatV = NULL;
}

TTSESinglePhaseTableClass::TTSESinglePhaseTableClass(Fluid *pFluid)
{
	this->pFluid = pFluid;
	// The default data location for the LUT
	this->root_path = std::string(format("TTSEData/%s",pFluid->get_name().c_str()));

	// Seed the generator for random number generation
	srand((unsigned int)time(NULL));
	this->enable_writing_tables_to_files = true;
	SatL = NULL;
	SatV = NULL;
}
void TTSESinglePhaseTableClass::set_size_ph(unsigned int Np, unsigned int Nh)
{
	this->Nh = Nh;
	this->Np = Np;

	h.resize(Nh);
	p.resize(Np);

	s.resize(Nh, std::vector<double>(Np, _HUGE));
	dsdh.resize(Nh, std::vector<double>(Np, _HUGE));
	dsdp.resize(Nh, std::vector<double>(Np, _HUGE));
	d2sdh2.resize(Nh, std::vector<double>(Np, _HUGE));	
	d2sdp2.resize(Nh, std::vector<double>(Np, _HUGE));
	d2sdhdp.resize(Nh, std::vector<double>(Np, _HUGE));

	T.resize(Nh, std::vector<double>(Np, _HUGE));
	dTdh.resize(Nh, std::vector<double>(Np, _HUGE));
	dTdp.resize(Nh, std::vector<double>(Np, _HUGE));
	d2Tdh2.resize(Nh, std::vector<double>(Np, _HUGE));	
	d2Tdp2.resize(Nh, std::vector<double>(Np, _HUGE));
	d2Tdhdp.resize(Nh, std::vector<double>(Np, _HUGE));

	rho.resize(Nh, std::vector<double>(Np, _HUGE));
	drhodh.resize(Nh, std::vector<double>(Np, _HUGE));
	drhodp.resize(Nh, std::vector<double>(Np, _HUGE));
	d2rhodh2.resize(Nh, std::vector<double>(Np, _HUGE));	
	d2rhodp2.resize(Nh, std::vector<double>(Np, _HUGE));
	d2rhodhdp.resize(Nh, std::vector<double>(Np, _HUGE));

	iL.resize(Np);
	iV.resize(Np);
	TL.resize(Np);
	TV.resize(Np);
	SL.resize(Np);
	SV.resize(Np);
	DL.resize(Np);
	DV.resize(Np);
}

void TTSESinglePhaseTableClass::set_size_Trho(unsigned int NT, unsigned int Nrho)
{
	this->NT = NT;
	this->Nrho = Nrho;

	T_Trho.resize(NT);
	rho_Trho.resize(Nrho);

	s_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dsdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dsdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2sdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2sdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2sdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));

	p_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dpdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dpdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2pdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2pdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2pdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));

	h_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dhdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dhdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2hdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2hdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2hdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
}

void TTSESinglePhaseTableClass::nearest_good_neighbor(int *i, int *j)
{
	// Left
	if (*i>0 && ValidNumber(rho[*i-1][*j]) && ValidNumber(T[*i-1][*j])){
		*i -= 1;
		return;
	}
	// Right
	else if (*i<(int)Nh-1 && ValidNumber(rho[*i+1][*j]) && ValidNumber(T[*i+1][*j])){
		*i += 1;
		return;
	}
	// Down
	else if (*j>0 && ValidNumber(rho[*i][*j-1]) && ValidNumber(T[*i][*j-1])){
		*j -= 1;
		return;
	}
	// Up
	else if (*j<(int)Np-1 && ValidNumber(rho[*i][*j+1]) && ValidNumber(T[*i][*j+1])){
		*j += 1;
		return;
	}
	else
	{
		throw ValueError(format("No neighbors found for %d,%d",i,j));
		return;
	}
}

void TTSESinglePhaseTableClass::nearest_good_neighbor_Trho(int *i, int *j)
{
	// Left
	if (*i>0 && ValidNumber(h_Trho[*i-1][*j]) && ValidNumber(p_Trho[*i-1][*j])){
		*i -= 1;
		return;
	}
	// Right
	else if (*i<(int)Nh-1 && ValidNumber(h_Trho[*i+1][*j]) && ValidNumber(p_Trho[*i+1][*j])){
		*i += 1;
		return;
	}
	// Down
	else if (*j>0 && ValidNumber(h_Trho[*i][*j-1]) && ValidNumber(p_Trho[*i][*j-1])){
		*j -= 1;
		return;
	}
	// Up
	else if (*j<(int)Np-1 && ValidNumber(h_Trho[*i][*j+1]) && ValidNumber(p_Trho[*i][*j+1])){
		*j += 1;
		return;
	}
	else
	{
		throw ValueError(format("No neighbors found for %d,%d",i,j));
		return;
	}
}

void TTSESinglePhaseTableClass::nearest_neighbor_ph(int i, int j, double *T0, double *rho0)
{
	// Left
	if (i>0 && ValidNumber(rho[i-1][j]) && ValidNumber(T[i-1][j])){
		*T0 = T[i-1][j];
		*rho0 = rho[i-1][j];
		return;
	}
	// Right
	else if (i<(int)Nh-1 && ValidNumber(rho[i+1][j]) && ValidNumber(T[i+1][j])){
		*T0 = T[i+1][j];
		*rho0 = rho[i+1][j];
		return;
	}
	// Down
	else if (j>0 && ValidNumber(rho[i][j-1]) && ValidNumber(T[i][j-1])){
		*T0 = T[i][j-1];
		*rho0 = rho[i][j-1];
		return;
	}
	// Up
	else if (j<(int)Np-1 && ValidNumber(rho[i][j+1]) && ValidNumber(T[i][j+1])){
		*T0 = T[i][j+1];
		*rho0 = rho[i][j+1];
		return;
	}
	else
	{
		*T0 = -1;
		*rho0 = -1;
		return;
	}
}
std::string join(std::vector<std::string> strings, char delim)
{
	std::string output = strings[0];
	for (unsigned int i = 1; i < strings.size(); i++)
	{
		output += format("%c%s",delim,strings[i].c_str());
	}
	return output;
}

void TTSESinglePhaseTableClass::matrix_to_file(std::string fName, std::vector< std::vector<double> > *A)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"wb");
	if (pFile != NULL)
	{
		for (unsigned int i = 0; i < Nh; i++)
		{
			fwrite(&(((*A)[i])[0]), sizeof(double),Np,pFile);
		}
		fclose(pFile);
	}
}

void TTSESinglePhaseTableClass::matrix_from_file(std::string fName, std::vector<std::vector<double> > *A)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"rb");
	if (pFile != NULL)
	{
		for (unsigned int i = 0; i < Nh; i++)
		{
			/*std::vector<double> row(Nh,0);
			fread(&(row[0]), sizeof(double), Np, pFile);
			double tgregt = 1;*/
			fread(&(((*A)[i])[0]), sizeof(double),Np,pFile);
		}
		fclose(pFile);
	}
}

void TTSESinglePhaseTableClass::vector_to_file(std::string fName, std::vector<double>*A)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"wb");
	if (pFile != NULL)
	{
		fwrite((const char *)&(*A).front(), sizeof(double),(*A).size(),pFile);
		fclose(pFile);
	}
}
void TTSESinglePhaseTableClass::vector_from_file(std::string fName, int N, std::vector<double> *vec)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"rb");
	if (pFile != NULL)
	{
		fread(&((*vec)[0]), sizeof(double), N, pFile);
		fclose(pFile);
	}
}
bool pathExists(std::string path)
{
	#if defined(_WIN32) // Defined for 32-bit and 64-bit windows
		struct _stat buf;
		// Get data associated with path using the windows libraries, 
		// and if you can (result == 0), the path exists
		if ( _stat( path.c_str(), &buf) == 0)
			return true;
		else
			return false;
	#endif
}
bool fileExists(const char *fileName)
{
	std::ifstream infile(fileName);
    return infile.good();
}
void make_dirs(std::string file_path)
{
	std::vector<std::string> pathsplit = strsplit(file_path,'/');
	std::string path = pathsplit[0];
	if (pathsplit.size()>0)
	{
		for (unsigned int i = 0; i < pathsplit.size(); i++)
		{
			if (!pathExists(path))
			{
			#ifdef _WIN32
				#if defined(_UNICODE)
					CreateDirectoryA((LPCSTR)path.c_str(),NULL);
				#else
					CreateDirectory((LPCSTR)path.c_str(),NULL);
				#endif
			#else 
				#if defined(__powerpc__)
				#else
					mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				#endif
			#endif
			}
			if (i < pathsplit.size()-1)
				path += format("/%s",pathsplit[i+1].c_str());
		}
	}
	else
	{
		throw ValueError(format("Could not make path [%s],file_path"));
	}
}

std::string get_file_contents(std::string filename)
{
	std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
	if (in)
	{
		std::string contents;
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	throw(errno);
}

bool TTSESinglePhaseTableClass::read_all_from_file(std::string root_path)
{
	std::string Fluid;
	double hmin,hmax,pmin,pmax;
	int Np, Nh;

	// Replace any '\' with '/' in the path
	for (unsigned int i = 0; i<root_path.length(); i++)
	{
		if (root_path[i] == '\\') root_path[i] = '/';
	}
	// If it ends with a '/', remove it for now
	if (root_path[root_path.size()-1] == '/')
		root_path = std::string(root_path,0,root_path.size()-1);

	if (!pathExists(root_path))
	{
		return false;
	}

	// Append a '/'
	root_path += format("/");

	if (!pathExists(root_path+format("Info_ph.txt")))
	{
		return false;
	}

	// Parse the info file to find the dimensions and boundaries and check that they are the same
	std::string info = get_file_contents((root_path+format("Info_ph.txt")).c_str());

	// Split into lines
	std::vector<std::string> lines = strsplit(info,'\r');

	for (unsigned int i = 0; i< lines.size(); i++)
	{
		// Split at ':'
		std::vector<std::string> line = strsplit(lines[i],':');
		if (line.size() == 1) // No : found
			continue;

		if (line[0].find("Fluid")!=  std::string::npos){
			Fluid = line[1];}
		else if (line[0].find("hmin")!=  std::string::npos) {	
			hmin = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("hmax")!=  std::string::npos) {	
			hmax = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("pmin")!=  std::string::npos) {	
			pmin = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("pmax")!=  std::string::npos) {	
			pmax = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("Np")!=  std::string::npos) {	
			Np = (int)strtol(line[1].c_str(),NULL,0);}
		else if (line[0].find("Nh")!=  std::string::npos) {	
			Nh = (int)strtol(line[1].c_str(),NULL,0);}
	}
	
	// Didn't work since at least one of the parameters was different
	// so we need to build the tables again
	if (!(!Fluid.compare(pFluid->get_name()) 
		  && Nh == this->Nh
		  && Np == this->Np
		  && fabs(pmin - this->pmin)<10*DBL_EPSILON
		  && fabs(pmax - this->pmax)<10*DBL_EPSILON 
		  && fabs(hmin - this->hmin)<10*DBL_EPSILON 
		  && fabs(hmax - this->hmax)<10*DBL_EPSILON
		)) return false;

	// Read all the data from the binary files
	vector_from_file(root_path + std::string("p_ph.ttse"),Np,&p);
	vector_from_file(root_path + std::string("h_ph.ttse"),Nh,&h);
	matrix_from_file(root_path + std::string("T_ph.ttse"),&T);
	matrix_from_file(root_path + std::string("dTdh_ph.ttse"),&dTdh);
	matrix_from_file(root_path + std::string("dTdp_ph.ttse"),&dTdp);
	matrix_from_file(root_path + std::string("d2Tdh2_ph.ttse"),&d2Tdh2);
	matrix_from_file(root_path + std::string("d2Tdp2_ph.ttse"),&d2Tdp2);
	matrix_from_file(root_path + std::string("d2Tdhdp_ph.ttse"),&d2Tdhdp);
	matrix_from_file(root_path + std::string("s_ph.ttse"),&s);
	matrix_from_file(root_path + std::string("dsdh_ph.ttse"),&dsdh);
	matrix_from_file(root_path + std::string("dsdp_ph.ttse"),&dsdp);
	matrix_from_file(root_path + std::string("d2sdh2_ph.ttse"),&d2sdh2);
	matrix_from_file(root_path + std::string("d2sdp2_ph.ttse"),&d2sdp2);
	matrix_from_file(root_path + std::string("d2sdhdp_ph.ttse"),&d2sdhdp);
	matrix_from_file(root_path + std::string("rho_ph.ttse"),&rho);
	matrix_from_file(root_path + std::string("drhodh_ph.ttse"),&drhodh);
	matrix_from_file(root_path + std::string("drhodp_ph.ttse"),&drhodp);
	matrix_from_file(root_path + std::string("d2rhodh2_ph.ttse"),&d2rhodh2);
	matrix_from_file(root_path + std::string("d2rhodp2_ph.ttse"),&d2rhodp2);
	matrix_from_file(root_path + std::string("d2rhoTdhdp_ph.ttse"),&d2rhodhdp);

	this->pratio = pow(pmax/pmin,1/((double)Np-1));
	this->logpratio = log(pratio); // For speed since log() is a slow function
	this->jpcrit_floor = (int)floor((log(pFluid->reduce.p)-logpmin)/logpratio);
	this->jpcrit_ceil = (int)ceil((log(pFluid->reduce.p)-logpmin)/logpratio);
	update_saturation_boundary_indices();

	return true;
}
void TTSESinglePhaseTableClass::write_all_to_file(std::string root_path)
{
	// Replace any '\' with '/' in the path
	for (unsigned int i = 0; i<root_path.length(); i++)
	{
		if (root_path[i] == '\\') root_path[i] = '/';
	}
	// If it ends with a '/', remove it for now
	if (root_path[root_path.size()-1] == '/')
		root_path = std::string(root_path,0,root_path.size()-1);

	if (!pathExists(root_path))
		make_dirs(root_path);

	// Append a '/'
	root_path += format("/");

	std::string header = std::string("Data for the TTSE method\nDO NOT CHANGE ANY OF THESE PARAMETERS FOR ANY REASON!\n\n");
		
	header += format("Fluid:%s\npmin:%23.19g\npmax:%23.19g\nNp:%25d\nhmin:%23.19g\nhmax:%23.19g\nNh:%25d\n",pFluid->get_name().c_str(),pmin,pmax,Np,hmin,hmax,Nh);
	
	clock_t t1,t2;
	t1 = clock();

	// Write the header information to a text file
	FILE *fp;
	fp = fopen((root_path+std::string("Info_ph.txt")).c_str(),"w");
	fprintf(fp,"%s",header.c_str());
	fclose(fp);

	// Write each of these files in binary mode
	vector_to_file(root_path + std::string("p_ph.ttse"),&p);
	vector_to_file(root_path + std::string("h_ph.ttse"),&h);
	matrix_to_file(root_path + std::string("T_ph.ttse"),&T);
	matrix_to_file(root_path + std::string("dTdh_ph.ttse"),&dTdh);
	matrix_to_file(root_path + std::string("dTdp_ph.ttse"),&dTdp);
	matrix_to_file(root_path + std::string("d2Tdh2_ph.ttse"),&d2Tdh2);
	matrix_to_file(root_path + std::string("d2Tdp2_ph.ttse"),&d2Tdp2);
	matrix_to_file(root_path + std::string("d2Tdhdp_ph.ttse"),&d2Tdhdp);
	matrix_to_file(root_path + std::string("s_ph.ttse"),&s);
	matrix_to_file(root_path + std::string("dsdh_ph.ttse"),&dsdh);
	matrix_to_file(root_path + std::string("dsdp_ph.ttse"),&dsdp);
	matrix_to_file(root_path + std::string("d2sdh2_ph.ttse"),&d2sdh2);
	matrix_to_file(root_path + std::string("d2sdp2_ph.ttse"),&d2sdp2);
	matrix_to_file(root_path + std::string("d2sdhdp_ph.ttse"),&d2sdhdp);
	matrix_to_file(root_path + std::string("rho_ph.ttse"),&rho);
	matrix_to_file(root_path + std::string("drhodh_ph.ttse"),&drhodh);
	matrix_to_file(root_path + std::string("drhodp_ph.ttse"),&drhodp);
	matrix_to_file(root_path + std::string("d2rhodh2_ph.ttse"),&d2rhodh2);
	matrix_to_file(root_path + std::string("d2rhodp2_ph.ttse"),&d2rhodp2);
	matrix_to_file(root_path + std::string("d2rhoTdhdp_ph.ttse"),&d2rhodhdp);
	t2 = clock();
	std::cout << "write time: " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;
}

double TTSESinglePhaseTableClass::build_ph(double hmin, double hmax, double pmin, double pmax, TTSETwoPhaseTableClass *SatL, TTSETwoPhaseTableClass *SatV)
{
	bool SinglePhase = false;

	this->hmin = hmin;
	this->hmax = hmax;
	this->pmin = pmin;
	this->pmax = pmax;
	this->logpmin = log(pmin);

	// If we can read them, we are done and don't need to rebuild
	if (read_all_from_file(root_path))
		return 0;

	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	double dh = (hmax - hmin)/(Nh - 1);
	pratio = pow(pmax/pmin,1/((double)Np-1));
	logpratio = log(pratio);

	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i<Nh; i++)
	{
		double hval = hmin + i*dh;
		h[i] = hval;
		for (unsigned int j = 0; j<Np; j++)
		{
			double pval = pmin*pow(pratio,(int)j);
			p[j] = pval;
			
			// Check whether the point is single phase
			// If pressure between ptriple point and pcrit, might be two-phase or single phase, otherwise definitely single phase
			if (pval <= pFluid->reduce.p && pval >= pFluid->params.ptriple)
			{
				if (SatL == NULL || SatV == NULL){
					// Not using TTSE method, use saturation (slow...)
					CPS.update(iP,pval,iQ,0.5);
					SinglePhase = (hval < CPS.hL() || hval > CPS.hV());
				}
				else{
					// Using the TTSE method, nice and fast
					SinglePhase = (hval < SatL->evaluate(iH,pval)  || hval > SatV->evaluate(iH,pval));
				}
			}
			else
			{
				SinglePhase = true;
			}
			
			// If enthalpy is outside the saturation region, continue and do the calculation as a function of p,h
			if (SinglePhase)
			{
				double T0=-1,rho0=-1,T,rho,rhoL,rhoV,TsatL,TsatV;

				// Find a good point around this point that is single-phase if any of its neighbors have been calculated
				nearest_neighbor_ph(i,j,&T0,&rho0);

				// If good T,rho was returned, use it as a guess value to calculate T,rho from p,h more quickly
				if (T0 > 0 && rho0 > 0)
				{
					// Get T,rho from p,h using our guess value
					CPS.pFluid->temperature_ph(pval,hval,&T,&rho,&rhoL,&rhoV,&TsatL,&TsatV,T0,rho0);
					CPS.flag_SinglePhase = true;
					CPS.update(iT, T, iD, rho);
				}
				else
				{
					// Probably here you are close to the saturation boundary since TTSE says it is single phase, but no neighbors to lean on
					CPS.flag_SinglePhase = true;

					double hsatLTTSE,hsatVTTSE;
					if (SatL == NULL || SatV == NULL){
						// Not using TTSE method, use saturation (slow...)
						CPS.update(iP,pval,iQ,0.5);
						hsatLTTSE = CPS.hL();
						hsatVTTSE = CPS.hV();
					}	
					else
					{
						hsatLTTSE = SatL->evaluate(iH,pval);
						hsatVTTSE = SatV->evaluate(iH,pval);
					}

					if (fabs(hval-hsatLTTSE)<10)
					{
						// Close to the saturated liquid boundary
						T0 = SatL->evaluate(iT,pval);
						rho0 = SatL->evaluate(iD,pval);
						CPS.pFluid->temperature_ph(pval,hval,&T,&rho,&rhoL,&rhoV,&TsatL,&TsatV,T0,rho0);
						CPS.flag_SinglePhase = true;
						CPS.update(iT, T, iD, rho);
					}
					else if (fabs(hval-hsatVTTSE)<10)
					{
						// Close to the saturated vapor boundary
						T0 = SatV->evaluate(iT,pval);
						rho0 = SatV->evaluate(iD,pval);
						CPS.pFluid->temperature_ph(pval,hval,&T,&rho,&rhoL,&rhoV,&TsatL,&TsatV,T0,rho0);
						CPS.flag_SinglePhase = true;
						CPS.update(iT, T, iD, rho);
					}
					else
					{
						CPS.update(iP, pval, iH, hval);
					}
					//std::cout << format("%d %d %g %g\n",i,j,pval,hval);
				}
				
				T = CPS.T();
				rho = CPS.rho();
				double cp = CPS.cp();

				double A = CPS.dpdT_constrho()*CPS.dhdrho_constT()-CPS.dpdrho_constT()*CPS.dhdT_constrho();

				this->T[i][j] = T;
				dTdh[i][j] = 1/cp;
				dTdp[i][j] = 1/A*CPS.dhdrho_constT();
				this->rho[i][j] = rho;
				drhodh[i][j] = 1/A*CPS.dpdT_constrho();
				drhodp[i][j] = -1/A*CPS.dhdT_constrho();
				s[i][j] = CPS.s();
				dsdh[i][j] = 1/T;
				dsdp[i][j] = -1/(T*rho);
				
				// Matrices for second derivatives of entropy as a function of pressure and enthalpy
				d2sdh2[i][j] = -1/(T*T)*dTdh[i][j];
				d2sdhdp[i][j] = -1/(T*T)*dTdp[i][j];
				d2sdp2[i][j] = 1/(T*T*rho)*dTdp[i][j]+1/(T*rho*rho)*drhodp[i][j];

				// These are common terms needed for a range of terms for T(h,p) as well as rho(h,p)
				double dAdT_constrho = CPS.d2pdT2_constrho()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrhodT()-CPS.d2pdrhodT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdT2_constrho();
				double dAdrho_constT = CPS.d2pdrhodT()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrho2_constT()-CPS.d2pdrho2_constT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdrhodT();

				//Matrices for temperature as a function of pressure and enthalpy
				double ddT_dTdp_h_constrho = 1/A*CPS.d2hdrhodT()-1/(A*A)*dAdT_constrho*CPS.dhdrho_constT(); //[check]
				double ddrho_dTdp_h_constT = 1/A*CPS.d2hdrho2_constT()-1/(A*A)*dAdrho_constT*CPS.dhdrho_constT(); //[check
				double ddT_dTdh_p_constrho = -1/(cp*cp)*(CPS.d2hdT2_constrho()-CPS.dhdp_constT()*CPS.d2pdT2_constrho()+CPS.d2hdrhodT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrhodT());
				double ddrho_dTdh_p_constT = -1/(cp*cp)*(CPS.d2hdrhodT()-CPS.dhdp_constT()*CPS.d2pdrhodT()+CPS.d2hdrho2_constT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrho2_constT());
				
				d2Tdh2[i][j]  = ddT_dTdh_p_constrho/CPS.dhdT_constp()+ddrho_dTdh_p_constT/CPS.dhdrho_constp();
				d2Tdhdp[i][j] = ddT_dTdp_h_constrho/CPS.dhdT_constp()+ddrho_dTdp_h_constT/CPS.dhdrho_constp();
				d2Tdp2[i][j]  = ddT_dTdp_h_constrho/CPS.dpdT_consth()+ddrho_dTdp_h_constT/CPS.dpdrho_consth();

				//// Matrices for density as a function of pressure and enthalpy
				double ddT_drhodp_h_constrho = -1/A*CPS.d2hdT2_constrho()+1/(A*A)*dAdT_constrho*CPS.dhdT_constrho();
				double ddrho_drhodp_h_constT = -1/A*CPS.d2hdrhodT()+1/(A*A)*dAdrho_constT*CPS.dhdT_constrho();
				double ddT_drhodh_p_constrho = 1/A*CPS.d2pdT2_constrho()-1/(A*A)*dAdT_constrho*CPS.dpdT_constrho();
				double ddrho_drhodh_p_constT = 1/A*CPS.d2pdrhodT()-1/(A*A)*dAdrho_constT*CPS.dpdT_constrho();
				
				d2rhodh2[i][j]  = ddT_drhodh_p_constrho/CPS.dhdT_constp()+ddrho_drhodh_p_constT/CPS.dhdrho_constp();
				d2rhodhdp[i][j] = ddT_drhodp_h_constrho/CPS.dhdT_constp()+ddrho_drhodp_h_constT/CPS.dhdrho_constp();
				d2rhodp2[i][j]  = ddT_drhodp_h_constrho/CPS.dpdT_consth()+ddrho_drhodp_h_constT/CPS.dpdrho_consth();

				/// Transport properties
				///
				if (transport_properties)
				{

				}
			}
			else
			{
				s[i][j] = _HUGE;
				dsdh[i][j] = _HUGE;
				dsdp[i][j] = _HUGE;
				d2sdh2[i][j] = _HUGE;
				d2sdhdp[i][j] = _HUGE;
				d2sdp2[i][j] = _HUGE;

				this->T[i][j] = _HUGE;
				dTdh[i][j] = _HUGE;
				dTdp[i][j] = _HUGE;
				d2Tdh2[i][j]  = _HUGE;
				d2Tdhdp[i][j] = _HUGE;
				d2Tdp2[i][j]  = _HUGE;

				this->rho[i][j] = _HUGE;
				drhodh[i][j] = _HUGE;
				drhodp[i][j] = _HUGE;
				d2rhodh2[i][j]  = _HUGE;
				d2rhodhdp[i][j] = _HUGE;
				d2rhodp2[i][j]  = _HUGE;
			}
		}
	}
	t2 = clock();
	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
	std::cout << elap << " to build single phase table" << std::endl;

	// Update the boundaries of the points within the single-phase regions
	update_saturation_boundary_indices();
	
	if (enable_writing_tables_to_files){
		write_all_to_file(root_path);
	}
	return elap;
}
double TTSESinglePhaseTableClass::build_Trho(double Tmin, double Tmax, double rhomin, double rhomax, TTSETwoPhaseTableClass *SatL, TTSETwoPhaseTableClass *SatV)
{
	bool SinglePhase = false;

	if (Tmin < 0 && Tmax < 0 && rhomin < 0 && rhomax < 0)
	{
		rhomin = 9e9;
		rhomax = 0;
		Tmin = 9e9;
		Tmax = 0;
		// Use single-phase table to figure out the range for T,rho
		for (unsigned int i = 0; i<Nh; i++)
		{
			for (unsigned int j = 0; j<Np; j++)
			{
				if (ValidNumber(rho[i][j]) && rho[i][j] > rhomax){
					rhomax = rho[i][j];
				}
				if (ValidNumber(rho[i][j]) && rho[i][j] < rhomin){
					rhomin = rho[i][j];
				}
				if (ValidNumber(T[i][j]) && T[i][j] > Tmax){
					Tmax = T[i][j];
				}
				if (ValidNumber(T[i][j]) && T[i][j] < Tmin){
					Tmin = T[i][j];
				}
			}
		}
	}
	this->Tmin = Tmin;
	this->Tmax = Tmax;
	this->rhomin = rhomin;
	this->rhomax = rhomax;

	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	double dT = (Tmax - Tmin)/((double)NT - 1);
	double drho = (rhomax - rhomin)/((double)Nrho - 1);

	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i<NT; i++)
	{
		double Tval = Tmin + i*dT;
		T_Trho[i] = Tval;
		for (unsigned int j = 0; j<Nrho; j++)
		{
			double rhoval = rhomin + j*drho;
			rho_Trho[j] = rhoval;
			
			// Check whether the point is single phase
			// If pressure between Ttriple point and Tcrit, might be two-phase or single phase, otherwise definitely single phase
			if (Tval <= pFluid->crit.T && Tval >= pFluid->params.Ttriple)
			{
				if (SatL == NULL || SatV == NULL){
					// Not using TTSE method, use saturation (slow...)
					CPS.update(iT,Tval,iQ,0.5);
					SinglePhase = (rhoval < CPS.rhoV() || rhoval > CPS.rhoL());
				}
				else{
					// Using the TTSE method, nice and fast
					double psatV = SatV->evaluate_T(Tval);
					double psatL = SatL->evaluate_T(Tval);
					SinglePhase = (rhoval < SatV->evaluate(iD,psatV)  || rhoval > SatL->evaluate(iD,psatL));
				}
			}
			else
			{
				SinglePhase = true;
			}
			
			// If enthalpy is outside the saturation region, continue and do the calculation as a function of T,rho
			if (SinglePhase)
			{
				CPS.update(iT,Tval,iD,rhoval);

				s_Trho[i][j] = CPS.s();
				dsdT_Trho[i][j] = CPS.dsdT_constrho();
				dsdrho_Trho[i][j] = CPS.dsdrho_constT();
				d2sdT2_Trho[i][j] = CPS.d2sdT2_constrho();
				d2sdTdrho_Trho[i][j] = CPS.d2sdrhodT();
				d2sdrho2_Trho[i][j] = CPS.d2sdrho2_constT();

				h_Trho[i][j] = CPS.h();
				dhdT_Trho[i][j] = CPS.dhdT_constrho();
				dhdrho_Trho[i][j] = CPS.dhdrho_constT();
				d2hdT2_Trho[i][j] = CPS.d2hdT2_constrho();
				d2hdTdrho_Trho[i][j] = CPS.d2hdrhodT();
				d2hdrho2_Trho[i][j] = CPS.d2hdrho2_constT();

				p_Trho[i][j] = CPS.p();
				dpdT_Trho[i][j] = CPS.dpdT_constrho();
				dpdrho_Trho[i][j] = CPS.dpdrho_constT();
				d2pdT2_Trho[i][j] = CPS.d2pdT2_constrho();
				d2pdTdrho_Trho[i][j] = CPS.d2pdrhodT();
				d2pdrho2_Trho[i][j] = CPS.d2pdrho2_constT();

				/// Transport properties
				///
				if (transport_properties)
				{

				}
			}
			else
			{
				s_Trho[i][j] = _HUGE;
				dsdT_Trho[i][j] = _HUGE;
				dsdrho_Trho[i][j] = _HUGE;
				d2sdT2_Trho[i][j] = _HUGE;
				d2sdTdrho_Trho[i][j] = _HUGE;
				d2sdrho2_Trho[i][j] = _HUGE;

				h_Trho[i][j] = _HUGE;
				dhdT_Trho[i][j] = _HUGE;
				dhdrho_Trho[i][j] = _HUGE;
				d2hdT2_Trho[i][j] = _HUGE;
				d2hdTdrho_Trho[i][j] = _HUGE;
				d2hdrho2_Trho[i][j] = _HUGE;

				p_Trho[i][j] = _HUGE;
				dpdT_Trho[i][j] = _HUGE;
				dpdrho_Trho[i][j] = _HUGE;
				d2pdT2_Trho[i][j] = _HUGE;
				d2pdTdrho_Trho[i][j] = _HUGE;
				d2pdrho2_Trho[i][j] = _HUGE;
			}
		}
	}
	t2 = clock();
	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
	std::cout << elap << " to build single phase table for T,rho" << std::endl;

	// Update the boundaries of the points within the single-phase regions
	update_saturation_boundary_indices();
	
	if (enable_writing_tables_to_files){
		write_all_to_file(root_path);
	}
	return elap;
}
//void TTSESinglePhaseTableClass::update_Trho_map()
//{
//	int ii,jj;
//	double Tmin,Tmax,rhomin,rhomax,rhoL,rhoV,TsatL,TsatV,dummy;
//	// Get the bounding values
//	Tmin = T[0][0];
//	rhomax = rho[0][0];
//	Tmax = T[Nh-1][Np-1];
//	rhomin = rho[Nh-1][0];
//	// Resize the arrays, using the same sizes as the base matrices
//	T_Trho.resize(Nh);
//	rho_Trho.resize(Np);
//	i_Trho.resize(Nh, std::vector<int>(Np, -1));
//	j_Trho.resize(Nh, std::vector<int>(Np, -1));
//	
//	for (unsigned int i = 0; i < Nh; i++)
//	{
//		double T = (Tmax-Tmin)/(Nh-1)*i+Tmin;
//		T_Trho[i] = T;
//
//		for (unsigned int j = 0; j < Np; j++)
//		{
//			double rho = (rhomax-rhomin)/(Np-1)*j+rhomin;
//			rho_Trho[j] = rho;
//
//			// T,rho --> p,h
//			double p = pFluid->pressure_Trho(T,rho);
//			double h = pFluid->enthalpy_Trho(T,rho);
//			double rhooV = pFluid->rhosatV(T);
//			double rhooL = pFluid->rhosatL(T);
//			double pV = pFluid->psatV_anc(T);
//			double pp = pFluid->pressure_Trho(T,rhooV);
//
//			// Find i,j from p,h
//			ii = (int)round(((h-hmin)/(hmax-hmin)*(Nh-1)));
//			jj = (int)round((log(p)-logpmin)/logpratio);
//
//			// Only keep values that are within the range for the table
//			if ( ii>=0 && ii < (int)Nh && jj>=0 && jj< (int)Np)
//			{
//				i_Trho[i][j] = ii;
//				j_Trho[i][j] = jj;
//			}
//			else
//			{
//				i_Trho[i][j] = -1;
//				j_Trho[i][j] = -1;
//			}
//		}
//	}
//}
void TTSESinglePhaseTableClass::update_saturation_boundary_indices()
{
	// Store some information about the phase boundaries so that we 
	// can use other inputs than p,h more easily

	for (unsigned int j = 0; j < Np; j++)
	{
		if (p[j] < pFluid->reduce.p)
		{
			iL[j] = -1;
			// Sweep left to right to find a phase boundary, use the first one that fails in the saturation region
			for (unsigned int i = 0; i < Nh; i++)
			{
				if (!ValidNumber(T[i][j]))
				{
					iL[j] = i;
					break;
				}
			}
			iV[j] = -1;
			// Sweep right to left to find a phase boundary, use the first one that fails in the saturation region
			for (int i = Nh-1; i > 0; i--)
			{
				if (!ValidNumber(T[i][j]))
				{
					iV[j] = i;
					break;	
				}
			}
			if (SatL != NULL && SatV != NULL)
			{
				TL[j] = SatL->evaluate(iT,p[j]);
				SL[j] = SatL->evaluate(iS,p[j]);
				DL[j] = SatL->evaluate(iD,p[j]);
				TV[j] = SatV->evaluate(iT,p[j]);
				SV[j] = SatV->evaluate(iS,p[j]);
				DV[j] = SatV->evaluate(iD,p[j]);
			}
			else
			{
				throw ValueError("SatL and SatV must be provided");
			}
		}
		else
		{
			iL[j] = -1;
			iV[j] = -1;
			TL[j] = _HUGE;
			SL[j] = _HUGE;
			DL[j] = _HUGE;
			TV[j] = _HUGE;
			SV[j] = _HUGE;
			DV[j] = _HUGE;
		}
	}
}
void TTSESinglePhaseTableClass::write_dotdrawing_tofile(char fName[])
{
	FILE *fp;
	fp = fopen(fName,"w");
	for (int j = Np-1; j>=0; j--)
	{
		for (unsigned int i = 0; i<Nh; i++)
		{
			if (ValidNumber(rho[i][j]))
			{
				fprintf(fp,".");
			}
			else
			{
				fprintf(fp,"X");
			}
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

double TTSESinglePhaseTableClass::check_randomly(long iParam, unsigned int N, std::vector<double> *h, std::vector<double> *p, std::vector<double> *EOSv, std::vector<double> *TTSE)
{	
	double val=0;
	h->resize(N);
	p->resize(N);
	EOSv->resize(N);
	TTSE->resize(N);
	
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;
		
		CPS.update(iH,h1,iP,p1);
		double sEOS = CPS.s();
		double cpEOS = CPS.cp();
		double TEOS = CPS.T();
		double rhoEOS = CPS.rho();

		// Store the inputs
		(*h)[i] = h1;
		(*p)[i] = p1;

		// Get the value from TTSE
		(*TTSE)[i] = evaluate(iParam,p1,h1);
		
		// Get the value from EOS
		switch (iParam)
		{
		case iS: 
			(*EOSv)[i] = sEOS; break;
		case iT:
			(*EOSv)[i] = TEOS; break;
		case iC:
			(*EOSv)[i] = cpEOS; break;
		case iD:
			(*EOSv)[i] = rhoEOS; break;
		default:
			throw ValueError();
		}
		
		std::cout << format("%g %g %g %g %g (h,p,EOS,TTSE,diff)\n",h1,p1,(*EOSv)[i],(*TTSE)[i],(*EOSv)[i]-(*TTSE)[i]);
	}
	return val;
}

double TTSESinglePhaseTableClass::evaluate_randomly(long iParam, unsigned int N)
{		
	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;

		if (p1 > pFluid->TTSESatL.pmax || h1 > pFluid->TTSESatV.evaluate(iH,p1) || h1 < pFluid->TTSESatL.evaluate(iH,p1))
		{
			// Get the value from TTSE
			evaluate(iParam,p1,h1);
		}
	}
	t2 = clock();
	return (double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6;
}

double TTSESinglePhaseTableClass::evaluate(long iParam, double p, double h)
{
	int i = (int)round(((h-hmin)/(hmax-hmin)*(Nh-1)));
	int j = (int)round((log(p)-logpmin)/logpratio);
	
	if (i<0 || i>(int)Nh-1 || j<0 || j>(int)Np-1)
	{
		throw ValueError(format("Input to TTSE [p = %0.16g, h = %0.16g] is out of range",p,h));
	}

	// If the value at i,j is too close to the saturation boundary, the nearest point i,j 
	// might be in the two-phase region which is not defined for single-phase table.  
	// Therefore, search around its neighbors for a better choice
	if (!ValidNumber(T[i][j])){
		nearest_good_neighbor(&i,&j);
	}

	// Distances from the node
	double deltap = p-this->p[j];
	double deltah = h-this->h[i];
	
	switch (iParam)
	{
	case iS:
		return s[i][j]+deltah*dsdh[i][j]+deltap*dsdp[i][j]+0.5*deltah*deltah*d2sdh2[i][j]+0.5*deltap*deltap*d2sdp2[i][j]+deltap*deltah*d2sdhdp[i][j]; break;
	case iT:
		return T[i][j]+deltah*dTdh[i][j]+deltap*dTdp[i][j]+0.5*deltah*deltah*d2Tdh2[i][j]+0.5*deltap*deltap*d2Tdp2[i][j]+deltap*deltah*d2Tdhdp[i][j]; break;
	case iD:
		return rho[i][j]+deltah*drhodh[i][j]+deltap*drhodp[i][j]+0.5*deltah*deltah*d2rhodh2[i][j]+0.5*deltap*deltap*d2rhodp2[i][j]+deltap*deltah*d2rhodhdp[i][j]; break;
	default:
		throw ValueError(format("Output key value [%d] to evaluate is invalid",iParam));
	}
}
bool isbetween(double x1, double x2, double x)
{
	return ((x>=x1 && x <= x2) || (x>=x2 && x <= x1));
}

double TTSESinglePhaseTableClass::evaluate_one_other_input(long iInput1, double Input1, long iOther, double Other)
{
	// My goal here is to find deltah
	int L,R,M,i;
	double p,dh1,dh2,a,b,c;
	std::vector<std::vector<double> > *mat;
	
	// Connect a pointer to the array of interest
	switch (iOther)
	{
	case iT:
		mat = &T; break;
	case iS:
		mat = &s; break;
	case iD:
		mat = &rho; break;
	}
	// One is pressure, we are getting enthalpy
	if (iInput1 == iP)
	{
		p = Input1;
		
		int j = (int)round((log(p)-logpmin)/logpratio);
		double deltap = p-this->p[j];

		if (j >= jpcrit_ceil) // Is is either supercritical pressure or just a little bit below critical pressure
		{
			// Very close to the boundary of the LUT, not 1-1 relationship between p-h and other
			// sets of inputs, need to allow for a bit of raggedness here
			if (   (iOther == iT && Other < 1.1*this->T[Np-1][j] && Other > this->T[Np-1][j])
			    || (iOther == iS && Other < this->s[Np-1][j]+0.1 && Other > this->s[Np-1][j])
			    || (iOther == iD && Other > 0.85*this->rho[Np-1][j] && Other < this->rho[Np-1][j])
			    )
			{
				i = Np-1;
			}
			else
			{
				// Do interval halving over the whole range since 
				// there can't be any saturation curve
				L = 0; R = Nh-1; M = (L+R)/2;
				while (R-L>1){
					if (isbetween((*mat)[M][j],(*mat)[R][j],Other)){ 
						L=M; M=(L+R)/2; continue;
					}
					else{ 
						R=M; M=(L+R)/2; continue;
					}
				}
				// Find which one of the bounds is closer
				if (fabs((*mat)[L][j]-Other)<fabs((*mat)[R][j]-Other)){
					i = L;
				}
				else{
					i = R;
				}
			}
		}
		else
		{
			if (   (iOther == iT && Other > SatV->evaluate(iT,p))
				|| (iOther == iS && Other > SatV->evaluate(iS,p))
				|| (iOther == iD && Other < SatV->evaluate(iD,p))
				)
			{
				// SUPERHEATED!!
				//
				// If it is within between the saturation curve and the first point in the SH region,
				// just use the first point in the superheated region
				if (   (iOther == iT && Other < this->T[iV[j]+1][j])
				    || (iOther == iS && Other < this->s[iV[j]+1][j])
					|| (iOther == iD && Other > this->rho[iV[j]+1][j])
					)
				{
					i = iV[j]+1;
				}
				// Very close to the boundary of the LUT, not 1-1 relationship between p-h and other
				// sets of inputs, need to allow for a bit of raggedness here
				else if (   (iOther == iT && Other < 1.1*this->T[Np-1][j] && Other > this->T[Np-1][j])
				         || (iOther == iS && Other < this->s[Np-1][j]+0.1 && Other > this->s[Np-1][j])
					     || (iOther == iD && Other > 0.85*this->rho[Np-1][j] && Other < this->rho[Np-1][j])
					     )
				{
					i = Np-1;
				}
				else
				{
					// SUPERHEATED!! (and away from saturation)
					// Make sure it is in the bounds
					switch (iOther)
					{
					case iT:
						if (Other > this->T[Nh-1][j]) 
							throw ValueError(format("Input T [%g] is greater than max T [%g] at LUT pressure [%g]",Other,this->T[Np-1][j],this->p[j]));
						break;
					case iD:
						if (Other < this->rho[Nh-1][j])
							throw ValueError(format("Input rho [%g] is less than minimum rho [%g] at LUT pressure [%g]",Other,this->rho[Np-1][j],this->p[j]));
						break;
					case iS:
						if (Other > this->s[Nh-1][j]) 
							throw ValueError(format("Input s [%g] is greater than max s [%g] at LUT pressure [%g]",Other,this->s[Np-1][j],this->p[j]));
						break;
					}
					
					L = iV[j]+1; R = Np-1; M = (L+R)/2;
					while (R-L>1)
					{
						if (isbetween((*mat)[M][j],(*mat)[R][j],Other))
						{ 
							L=M; M=(L+R)/2; continue;
						}
						else
						{ 
							R=M; M=(L+R)/2; continue;
						}
					}
					// Find which one of the bounds is closer
					if (fabs((*mat)[L][j]-Other)<fabs((*mat)[R][j]-Other))
						i = L;
					else
						i = R;
				}
			}
			else if (iL[j] < 2)
			{
				// We are at low pressure, so we don't know how to calculate, going to just use the i==1 element
				// if it is valid, or the i = 0 if not, otherwise, there are no values left and we have to fail
				if (ValidNumber(this->T[1][j])){
					i = 1;
				}
				else if (ValidNumber(this->T[0][j])){
					i = 0;
				}
				else{
					throw ValueError(format("Your inputs [%g,%g] do not yield a valid TTSE node",Input1,Other));
				}
			}
			
			else if (   (iOther == iT && Other < SatL->evaluate(iT,p))
				     || (iOther == iS && Other < SatL->evaluate(iS,p))
					 || (iOther == iD && Other > SatL->evaluate(iD,p))
					 )
			{
				// SUBCOOLED!!
				//
				// If it is within one spacing of the outlet variable of the saturation curve, 
				// just use the first point in the subcooled region
				if (   (iOther == iT && Other > this->T[iL[j]-1][j])
				    || (iOther == iS && Other > this->s[iL[j]-1][j])
					|| (iOther == iD && Other < this->rho[iL[j]-1][j])
					)
				{
					i = iL[j]-1;
				}
				else{
					// Make sure it is in the bounds of the LUT
					switch (iOther)
					{
					case iT:
						if (Other < this->T[0][j]) 
							throw ValueError(format("Input T [%g] is less than min T [%g] at LUT pressure [%g]",Other,this->T[0][j],this->p[j]));
						break;
					case iD:
						if (Other > this->rho[0][j]) 
							throw ValueError(format("Input rho [%g] is greater than max rho [%g] at LUT pressure [%g]",Other,this->rho[0][j],this->p[j]));
						break;
					case iS:
						if (Other < this->s[0][j])
							throw ValueError(format("Input s [%g] is less than min s [%g] at LUT pressure [%g]",Other,this->s[0][j],this->p[j]));
						break;
					}
					
					L = 0; R = iL[j]-1; M = (L+R)/2;
					// Its subcooled
					while (R-L>1)
					{
						if (isbetween((*mat)[M][j],(*mat)[R][j],Other))
						{ 
							L=M; M=(L+R)/2; continue;
						}
						else
						{ 
							R=M; M=(L+R)/2; continue;
						}
					}
					// Find which one of the bounds is closer
					if (fabs((*mat)[L][j]-Other)<fabs((*mat)[R][j]-Other))
						i = L;
					else
						i = R;
				}
			}
			else
			{
				// It's two-phase
				throw ValueError(format("It's two phase input"));
			}
		}

		// Now we calculate deltah
		switch (iOther)
		{
		// Quadratic in deltah
		// 0 = 0.5*deltah*deltah*d2Tdh2[i][j]+deltah*dTdh[i][j]+deltap*deltah*d2Tdhdp[i][j]+T[i][j]-T+deltap*dTdp[i][j]+0.5*deltap*deltap*d2Tdp2[i][j];
		// 0 = a*deltah^2+b*deltah+c
		case iT:
			a = 0.5*d2Tdh2[i][j];
			b = dTdh[i][j]+deltap*d2Tdhdp[i][j];
			c = T[i][j]-Other+deltap*dTdp[i][j]+0.5*deltap*deltap*d2Tdp2[i][j];
			break;
		case iS:
			a = 0.5*d2sdh2[i][j];
			b = dsdh[i][j]+deltap*d2sdhdp[i][j];
			c = s[i][j]-Other+deltap*dsdp[i][j]+0.5*deltap*deltap*d2sdp2[i][j];
			break;
		case iD:
			a = 0.5*d2rhodh2[i][j];
			b = drhodh[i][j]+deltap*d2rhodhdp[i][j];
			c = rho[i][j]-Other+deltap*drhodp[i][j]+0.5*deltap*deltap*d2rhodp2[i][j];
			break;
		}
		// Solutions from quadratic equation
		dh1 = (-b+sqrt(b*b-4*a*c))/(2*a);
		dh2 = (-b-sqrt(b*b-4*a*c))/(2*a);

		double hspacing = (hmax-hmin)/((double)Nh-1);
		// If only one is less than a multiple of enthalpy spacing, thats your solution
		if (fabs(dh1) < 10*hspacing && !(fabs(dh2) < 10*hspacing) )
			return dh1+this->h[i];
		else if (fabs(dh2) < 10*hspacing && !(fabs(dh1) < 10*hspacing) )
			return dh2+this->h[i];
		else{
			// Need to figure out which is the correct solution, try just the smaller one
			if (fabs(dh1)<fabs(dh2))
				return dh1+this->h[i];
			else
				return dh2+this->h[i];
		}
	}
	// One is enthalpy, we are getting pressure
	else if (iInput1 == iH)
	{
		throw ValueError("Sorry enthalpy and something else other than p is not valid for TTSE");
	}
	// Oops, neither enthalpy or pressure provided
	else
	{
		throw ValueError("Neither enthalpy nor pressure provided as the first parameter");
	}
}

double TTSESinglePhaseTableClass::evaluate_two_other_inputs(long iOutput, long iInput1, double Input1, long iInput2, double Input2)
{
	if (iInput1 == iT && iInput2 == iD)
	{
		int i = (int)round(((Input1-Tmin)/(Tmax-Tmin)*(NT-1)));
		int j = (int)round(((Input2-rhomin)/(rhomax-rhomin)*(Nrho-1)));
		
		if (i<0 || i>(int)NT-1 || j<0 || j>(int)Nrho-1)
		{
			throw ValueError(format("Input to TTSE [T = %0.16g, rho = %0.16g] is out of range",Input1,Input2));
		}

		// If the value at i,j is too close to the saturation boundary, the nearest point i,j 
		// might be in the two-phase region which is not defined for single-phase table.  
		// Therefore, search around its neighbors for a better choice
		if (!ValidNumber(s_Trho[i][j])){
			nearest_good_neighbor_Trho(&i,&j);
		}

		// Distances from the node
		double deltaT = Input1 - this->T_Trho[i];
		double deltarho = Input2 - this->rho_Trho[j];
		
		switch (iOutput)
		{
		case iS:
			return s_Trho[i][j]+deltaT*dsdT_Trho[i][j]+deltarho*dsdrho_Trho[i][j]+0.5*deltaT*deltaT*d2sdT2_Trho[i][j]+0.5*deltarho*deltarho*d2sdrho2_Trho[i][j]+deltaT*deltarho*d2sdTdrho_Trho[i][j]; break;
		case iP:
			return p_Trho[i][j]+deltaT*dpdT_Trho[i][j]+deltarho*dpdrho_Trho[i][j]+0.5*deltaT*deltaT*d2pdT2_Trho[i][j]+0.5*deltarho*deltarho*d2pdrho2_Trho[i][j]+deltaT*deltarho*d2pdTdrho_Trho[i][j]; break;
		case iH:
			return h_Trho[i][j]+deltaT*dhdT_Trho[i][j]+deltarho*dhdrho_Trho[i][j]+0.5*deltaT*deltaT*d2hdT2_Trho[i][j]+0.5*deltarho*deltarho*d2hdrho2_Trho[i][j]+deltaT*deltarho*d2hdTdrho_Trho[i][j]; break;
		default:
			throw ValueError(format("Output key value [%d] to evaluate is invalid",iOutput));
		}
	}
	else
	{

	}
}

double TTSESinglePhaseTableClass::evaluate_first_derivative(long iOF, long iWRT, long iCONSTANT, double p, double logp, double h)
{
	int i = (int)round(((h-hmin)/(hmax-hmin)*(Nh-1)));
	int j = (int)round((logp-logpmin)/logpratio);

	if (i<0 || i>(int)Nh-1 || j<0 || j>(int)Np-1)
	{
		throw ValueError(format("Input to TTSE deriv [p = %0.16g, h = %0.16g] is out of range",p,h));
	}
	
	// If the value at i,j is too close to the saturation boundary, the nearest point i,j 
	// might be in the two-phase region which is not defined for single-phase table.  
	// Therefore, search around its neighbors for a better choice
	if (!ValidNumber(T[i][j]))
	{
		nearest_good_neighbor(&i,&j);
	}

	// Distances from the node
	double deltah = h-this->h[i];
	double deltap = p-this->p[j];
	
	// This is a first-order expansion of the derivative around the node point.
	//
	// Derivatives for constant p
	if (iOF == iT && iWRT == iH && iCONSTANT == iP)
	{
		// Derivative of T w.r.t. h for p constant (for cp, the constant-pressure specific heat)
		return dTdh[i][j]+deltah*d2Tdh2[i][j]+deltap*d2Tdhdp[i][j];
	}
	else if (iOF == iS && iWRT == iH && iCONSTANT == iP)
	{
		// Derivative of s w.r.t. h for p constant
		return dsdh[i][j]+deltah*d2sdh2[i][j]+deltap*d2sdhdp[i][j];
	}
	else if (iOF == iD && iWRT == iH && iCONSTANT == iP)
	{
		// Derivative of density w.r.t. h for p constant
		return drhodh[i][j]+deltah*d2rhodh2[i][j]+deltap*d2rhodhdp[i][j];
	}

	// Derivatives for constant h
	else if (iOF == iT && iWRT == iP && iCONSTANT == iH)
	{
		// Derivative of T w.r.t. p for h constant
		return dTdp[i][j]+deltap*d2Tdp2[i][j]+deltah*d2Tdhdp[i][j];
	}
	else if (iOF == iS && iWRT == iP && iCONSTANT == iH)
	{
		// Derivative of s w.r.t. p for h constant
		return dsdp[i][j]+deltap*d2sdp2[i][j]+deltah*d2sdhdp[i][j];
	}
	else if (iOF == iD && iWRT == iP && iCONSTANT == iH)
	{
		// Derivative of density w.r.t. p for h constant
		return drhodp[i][j]+deltap*d2rhodp2[i][j]+deltah*d2rhodhdp[i][j];
	}
	else{
		throw ValueError(format("Your inputs [%d,%d,%d,%g,%g] are invalid to evaluate_first_derivative",iOF,iWRT,iCONSTANT,p,h));
	}
}

TTSETwoPhaseTableClass::TTSETwoPhaseTableClass(Fluid *pFluid, double Q)
{
	this->pFluid = pFluid;
	this->Q = Q;
}
void TTSETwoPhaseTableClass::set_size(unsigned int N)
{
	this->N = N;
	
	// Seed the generator
	srand((unsigned int)time(NULL));

	// Resize all the arrays
	h.resize(N);
	p.resize(N);
	logp.resize(N);
	T.resize(N);
	dTdp.resize(N);
	d2Tdp2.resize(N);
	rho.resize(N);
	logrho.resize(N);
	drhodp.resize(N);
	d2rhodp2.resize(N);
	s.resize(N);
	dsdp.resize(N);
	d2sdp2.resize(N);
	h.resize(N);
	dhdp.resize(N);
	d2hdp2.resize(N);
}

double TTSETwoPhaseTableClass::build(double pmin, double pmax, TTSETwoPhaseTableClass *other)
{
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	this->pmin = pmin;
	this->pmax = pmax;
	this->logpmin = log(pmin);
	this->logpmax = log(pmax);
	this->pratio = pow(pmax/pmin,1/((double)N-1));
	this->logpratio = log(pratio);

	double dlogp = (logpmax-logpmin)/(N-1);
	clock_t t1,t2;
	t1 = clock();
	// Logarithmic distribution of pressures
	for (unsigned int i = 0; i < N-1; i++)
	{
		// Calculate the pressure
		p[i] = exp(logpmin + i*dlogp);
		logp[i] = log(p[i]);
		// Update the class
		CPS.update(iP,p[i],iQ,Q);

		// Set the variables
		T[i] = CPS.T();
		dTdp[i] = CPS.dTdp_along_sat();
		d2Tdp2[i] = CPS.d2Tdp2_along_sat();
		h[i] = CPS.h();
		dhdp[i] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
		d2hdp2[i] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
		s[i] = CPS.s();
		dsdp[i] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
		d2sdp2[i] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
		rho[i] = CPS.rho();
		logrho[i] = log(CPS.rho());
		drhodp[i] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
		d2rhodp2[i] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();

		// If other is provided
		if (other != NULL)
		{
			other->pmin = pmin;
			other->pmax = pmax;
			other->logpmin = log(pmin);
			other->logpmax = log(pmax);

			other->p[i] = this->p[i];
			other->logp[i] = log(this->p[i]);

			// Set the variables
			other->T[i] = CPS.T();
			other->dTdp[i] = CPS.dTdp_along_sat();
			other->d2Tdp2[i] = CPS.d2Tdp2_along_sat();
			other->h[i] = (other->Q>0.5) ? CPS.hV() : CPS.hL();
			other->dhdp[i] = (other->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
			other->d2hdp2[i] = (other->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
			other->s[i] = (other->Q>0.5) ? CPS.sV() : CPS.sL();
			other->dsdp[i] = (other->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
			other->d2sdp2[i] = (other->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
			other->rho[i] = (other->Q>0.5) ? CPS.rhoV() : CPS.rhoL();
			other->logrho[i] = log(other->rho[i]);
			other->drhodp[i] = (other->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
			other->d2rhodp2[i] = (other->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();
		}
	}
	//// At the last point
	//CPS.flag_SinglePhase = true; // Don't have it check the state or do a saturation call
	//p[N-1] = CPS.pFluid->reduce.p;
	//logp[N-1] = log(p[N-1]);
	//CPS.update(iT,CPS.pFluid->reduce.T,iD,CPS.pFluid->reduce.rho);
	//T[N-1] = CPS.T();
	//dTdp[N-1] = CPS.dTdp_along_sat();
	//d2Tdp2[N-1] = CPS.d2Tdp2_along_sat();
	//h[N-1] = CPS.h();
	//dhdp[N-1] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
	//d2hdp2[N-1] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
	//s[N-1] = CPS.s();
	//dsdp[N-1] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
	//d2sdp2[N-1] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
	//rho[N-1] = CPS.rho();
	//logrho[N-1] = log(CPS.rho());
	//drhodp[N-1] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
	//d2rhodp2[N-1] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();

	t2 = clock();
	std::cout << double(t2-t1)/CLOCKS_PER_SEC << " to build both two phase tables" << std::endl;
	return double(t2-t1)/CLOCKS_PER_SEC;
}

double TTSETwoPhaseTableClass::evaluate(long iParam, double p)
{
	double logp = log(p);
	int i = (int)round(((logp-logpmin)/(logpmax-logpmin)*(N-1)));
	// If the value is just a little bit below the range, clip 
	// it back to the range of the LUT
	if (i == -1) i = 0;
	// If the pressure is just barely above the critical pressure
	// or between the critical pressure and the next lowest point,
	// just just the next lowest point to avoid some of the 
	// derivatives that are infinite at the critical point
	if (i == N || i == N-1 ) i = N-2;
	// If it is really out of the range, throw an error
	if (i<0 || i>(int)N-1)
	{
		throw ValueError(format("p [%g] is out of range[%g,%g], yielded index of: %d",p,pmin,pmax,i));
	}
		
	double log_PI_PIi = logp-this->logp[i];
	double pi = this->p[i];
	
	switch (iParam)
	{
	case iS:
		return s[i]+log_PI_PIi*pi*dsdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2sdp2[i]*pi*pi;
	case iT:
		return T[i]+log_PI_PIi*pi*dTdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2Tdp2[i]*pi*pi;
	case iH:
		//return h[i]+(pi-this->p[i])*dhdp[i]+0.5*(pi-this->p[i])*(pi-this->p[i])*d2hdp2[i];
		return h[i]+log_PI_PIi*pi*dhdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2hdp2[i]*pi*pi;
	case iD:
		{
		// log(p) v. log(rho) gives close to a line for most of the curve
		double dRdPI = pi/rho[i]*drhodp[i];
		return exp(logrho[i]+log_PI_PIi*(1.0+0.5*log_PI_PIi*(1-dRdPI))*dRdPI+0.5*log_PI_PIi*log_PI_PIi*d2rhodp2[i]*pi*pi/rho[i]);
		}
	default:
		throw ValueError();
	}
}
double TTSETwoPhaseTableClass::evaluate_T(double T)
{
	int L,R,M;
	double a,b,c,pi,log_PI_PIi1,log_PI_PIi2,logp_spacing;

	logp_spacing = this->logp[2]-this->logp[1];

	// Do interval halving over the whole range to find the nearest temperature
	L = 0; R = N - 2; M = (L+R)/2;
	if (isbetween(this->T[N-2],pFluid->reduce.T,T)){
		L = N-2;
	}
	else
	{
		while (R - L>1)
		{
			if (isbetween(this->T[M],this->T[R],T)){ 
				L=M; M=(L+R)/2; continue;
			}
			else{ 
				R=M; M=(L+R)/2; continue;
			}
		}
	}
	// T = T[i]+log_PI_PIi*pi*dTdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2Tdp2[i]*pi*pi;
	// T = T[i]+log_PI_PIi*pi*dTdp[i]+ 0.5*log_PI_PIi^2*pi*dTdp[i]+0.5*log_PI_PIi*log_PI_PIi*d2Tdp2[i]*pi*pi;
	// 0 = 0.5*log_PI_PIi^2*pi*dTdp[i]+0.5*log_PI_PIi^2*d2Tdp2[i]*pi*pi+log_PI_PIi*pi*dTdp[i]+T[i]-T;
	pi = this->p[L];
	a = 0.5*(pi*dTdp[L]+d2Tdp2[L]*pi*pi);
	b = pi*dTdp[L];
	c = this->T[L]-T;

	// Solutions from quadratic equation
	log_PI_PIi1 = (-b+sqrt(b*b-4*a*c))/(2*a);
	log_PI_PIi2 = (-b-sqrt(b*b-4*a*c))/(2*a);

	// Get the pressures
	double p1 = exp(log_PI_PIi1+this->logp[L]);
	double p2 = exp(log_PI_PIi2+this->logp[L]);

	// If only one is less than spacing of enthalpy, thats your solution
	if (fabs(log_PI_PIi1)<2*logp_spacing && !(fabs(log_PI_PIi2)<2*logp_spacing))
		return p1;
	else if (fabs(log_PI_PIi2)<2*logp_spacing && !(fabs(log_PI_PIi1)<2*logp_spacing))
		return p2;
	else
		throw ValueError(format("More than one solution found[%g,%g] in evaluate_T for TTSE for input %g",p1,p2,T));
}
double TTSETwoPhaseTableClass::evaluate_sat_derivative(long iParam, double p)
{
	double logp = log(p);
	int i = (int)round(((logp-logpmin)/(logpmax-logpmin)*(N-1)));
	// If the value is just a little bit out of the range, clip 
	// it back to the range of the LUT
	if (i == -1) i = 0;
	if (i == N) i = N-1;
	// If it is really out of the range, throw an error
	if (i<0 || i>(int)N-1)
	{
		throw ValueError(format("p [%g] is out of range",p));
	}
		
	double pi = this->p[i];
	
	switch (iParam)
	{
	case iT:
		{
			/// First order expansion of dTdp around point of interest
			return dTdp[i]+(p-pi)*d2Tdp2[i];
		}
	case iH:
		{	
			/// First order expansion of dhdp around point of interest
			return dhdp[i]+(p-pi)*d2hdp2[i];
		}
	case iS:
		{	
			/// First order expansion of dsdp around point of interest
			return dsdp[i]+(p-pi)*d2sdp2[i];
		}
	case iD:
		{
			/// First order expansion of drhodp around point of interest
			return drhodp[i]+(p-pi)*d2rhodp2[i];
		}
	default:
		throw ValueError(format("Cannot use the key [%d] provided in evaluate_sat_derivative",iParam));
	}
}

double TTSETwoPhaseTableClass::evaluate_randomly(long iParam, unsigned int N)
{		
	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;

		// Get the value from TTSE
		evaluate(iParam,p1);
	}
	t2 = clock();
	return (double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6;
}


double TTSETwoPhaseTableClass::check_randomly(long iParam, unsigned int N, std::vector<double> *p, std::vector<double> *EOSv, std::vector<double> *TTSE)
{	
	double val=0;
	p->resize(N);
	EOSv->resize(N);
	TTSE->resize(N);
	
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		
		CPS.update(iP,p1,iQ,this->Q);
		double hEOS = CPS.h();
		double sEOS = CPS.s();
		double TEOS = CPS.T();
		double rhoEOS = CPS.rho();

		// Store the inputs
		(*p)[i] = p1;

		// Get the value from TTSE
		(*TTSE)[i] = evaluate(iParam,p1);
		
		// Get the value from EOS
		switch (iParam)
		{
		case iS: 
			(*EOSv)[i] = sEOS; break;
		case iT:
			(*EOSv)[i] = TEOS; break;
		case iH:
			(*EOSv)[i] = hEOS; break;
		case iD:
			(*EOSv)[i] = rhoEOS; break;
		default:
			throw ValueError();
		}
		
		std::cout << format("%g %g %g %g TTSE (p,EOS,TTSE, diff)\n",p1,(*EOSv)[i],(*TTSE)[i],((*EOSv)[i]-(*TTSE)[i]));
	}
	return val;
}

