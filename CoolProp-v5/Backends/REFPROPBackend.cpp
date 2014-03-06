/*
 * AbstractBackend.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

#include <string>
#include "CoolProp.h"

#include "REFPROPBackend.h"

namespace CoolProp::Backends {

REFPROPBackend::REFPROPBackend(std::string fluid_name) {
	// Do the REFPROP instantiation for this fluid

	// Try to add this fluid to REFPROP - might want to think about making array of 
	// components and setting mole fractions if they change a lot.
	this->set_REFPROP_fluid(fluid_name);

	// Set all constants that can be accessed from REFPROP
	// Tcrit, pcrit, accentric...


}

REFPROPBackend::~REFPROPBackend() {
	// TODO Auto-generated destructor stub
}

void REFPROPBackend::set_REFPROP_fluid(std::string fluid_name)
{
	long ierr=0;
	char hf[refpropcharlength*ncmax], herr[errormessagelength+1];
	std::string sRef, components_joined;
	std::string RefString;
	std::string fdPath = get_REFPROP_fluid_path();

	// Check platform support
	if(!REFPROPFluidClass::refpropSupported()){
		throw NotImplementedError("You cannot use the REFPROPFluidClass.");
	}

	// Load REFPROP if it isn't loaded yet
	load_REFPROP();
	
	// If the name of the refrigerant doesn't match 
	// that of the currently loaded refrigerant
	if (LoadedREFPROPRef.compare(Ref))
	{
		// If the fluid name starts with the string "REFPROP-MIX:"
		if (Ref.find("REFPROP-MIX:") == 0)
		{
			// Keep everything after the "REFPROP-MIX:"
			components_joined = Ref.substr(12,Ref.size()-12);
		
			// Sample sRef is "R32[0.697615]&R125[0.302385]" -  this is R410A
			// Or you could do "R410A.mix" to use the full mixture model for this predefined mixture
			
			// Try to process predefined mixtures with .mix or .MIX in the file name
			if (components_joined.find(".mix") != std::string::npos || components_joined.find(".MIX") != std::string::npos)
			{
				char hf[255];
				char hfiles[10000];
				char herr[255];
				double xx[ncmax];
				strcpy(hf,components_joined.c_str());

				SETMIXdll(hf, hfmix, hrf, 
						  &i, hfiles, xx,
						  &ierr, herr,
						  255,
						  255,
						  3, // lengthofreference
						  10000,
						  255);
				// c-string needs to be 0-terminated
				for (unsigned int j = 0; j < 255*ncmax; j++)
				{
					if (hfiles[j] == 32) // empty char
					{
						hfiles[j] = 0;
						break;
					}
				}
                // Resize the vector of mole fractions
				x.resize(i);

				RefString = std::string(hfiles,strlen(hfiles)+1);
				for (int j = 0; j < i; j++)
				{
					x[j] = xx[j];
				}
			}
			else
			{
				// Split the components_joined into the components
				std::vector<std::string> components_split = strsplit(components_joined,'&');

				if (components_split.size() == 1)
				{
					throw ValueError(format("REFPROP mixture specified composition desired [%s], but only one component found",components_joined.c_str()).c_str());
				}

				// Flush out the refrigerant string for REFPROP
				RefString.clear();

				// Resize the vector of mole fractions
				x.resize(components_split.size());

				for (unsigned int j=0;j<components_split.size();j++)
				{	
					// Get component name and mole fraction (as strings)
					std::vector<std::string> comp_fraction = strsplit(components_split[j],'[');

					if (comp_fraction.size() != 2)
					{
						throw ValueError(format("Could not parse name[molefraction] [%s]",components_split[j].c_str()).c_str());
					}
					
					// Build the refrigerant string
					if (j == 0){
						RefString = fdPath + comp_fraction[0]+".fld";
					}
					else{
						RefString += "|" + fdPath + comp_fraction[0]+".fld";
					}
					// Convert the mole fraction (as string) to a number
					x[j] = strtod(comp_fraction[1].c_str(),NULL);

					// Update the number of components
					i = j+1;
				}
			}
		}
		// Name starts with REFPROP-
		else if (Ref.find("REFPROP-") == 0)
		{
			// Keep everything after the "REFPROP-"
			sRef = Ref.substr(8,Ref.size()-8);

			if (!sRef.compare("Air") || !sRef.compare("R507A") || !sRef.compare("R404A") || !sRef.compare("R410A") || !sRef.compare("R407C") || !sRef.compare("SES36"))
			{
				i=1;
				RefString = fdPath + std::string(sRef)+std::string(".ppf");
				x[0]=1.0;     //Pseudo-Pure fluid
			}
			else
			{
				i=1;
				RefString = fdPath + std::string(sRef)+std::string(".fld");
				x[0]=1.0;     //Pure fluid
			}
		}
		else
		{
			throw ValueError(format("REFPROP fluid string [%s] is invalid", Ref.c_str()));
		}

		ierr=999;
		// Set path to fluid files
//		// std::string rpPath (refpropfluidpath);
//		if (rpPath.length()>0)
//		{
//			printf("Setting REFPROP path to: %s\n",rpPath.c_str());
//			char refproppath[refpropcharlength+1];
//			strcpy(refproppath,rpPath.c_str());
//			SETPATHdll(refproppath);
//			free(refproppath);
//		}

		char* hfm = (char*) calloc(refpropcharlength+8, sizeof(char));
		strcpy(hfm,fdPath.c_str());
		strcat(hfm,hfmix);
		strcpy(hf,RefString.c_str());

		//...Call SETUP to initialize the program
		SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
		  refpropcharlength*ncmax,refpropcharlength,
		  lengthofreference,errormessagelength);

		if (ierr > 0){
			//...Call SETUP with capital letters
			for(int i = 0; i < strlen(hrf); i++)
			{
				hrf[i] = toupper(hrf[i]);
			}
			for(int i = 0; i < strlen(hfm); i++)
			{
				hfm[i] = toupper(hfm[i]);
			}
			for(int i = 0; i < strlen(hf); i++)
			{
				hf[i] = toupper(hf[i]);
			}
			SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
			  refpropcharlength*ncmax,refpropcharlength,
			  lengthofreference,errormessagelength);
		}

		if (ierr > 0){
			//...Call SETUP with lower case letters
			for(int i = 0; i < strlen(hrf); i++)
			{
				hrf[i] = tolower(hrf[i]);
			}
			for(int i = 0; i < strlen(hfm); i++)
			{
				hfm[i] = tolower(hfm[i]);
			}
			for(int i = 0; i < strlen(hf); i++)
			{
				hf[i] = tolower(hf[i]);
			}
			SETUPdll(&i, hf, hfm, hrf, &ierr, herr,
			  refpropcharlength*ncmax,refpropcharlength,
			  lengthofreference,errormessagelength);
		}

		free (hfm);

		if (ierr > 0){
			throw ValueError(format("REFPROP: %s",herr).c_str());
			return false;
		}
		else if (ierr < 0)
		{
			set_warning(herr);
		}
		//Copy the name of the loaded refrigerant back into the temporary holder
		LoadedREFPROPRef = std::string(Ref);
		
		unsigned int jmax;
		for (jmax = 0; jmax < ncmax; jmax++)
		{
			if (jmax == x.size())
			{
				break;
			}
			if (x[jmax] < 1e-13)
			{
				break;
			}
		}
		x.resize(jmax);
		LoadedREFPROPx = x;
		return true;
	}
	else
	{
		x = LoadedREFPROPx;
	}
	return true;
}
	
}
void REFPROPBackend::update(int input_pair, double value1, double value2)
{
	switch( input_pair)
	{
		case PT_INPUTS:
		{
			T = value1; p = value2/1000.0; // Want p in [kPa]

			// Use flash routine to find properties
			TPFLSHdll(&_T,&_p,&(x[0]),&d,&dl,&dv,xliq,xvap,&q,&e,&h,&s,&cv,&cp,&w,&ierr,herr,errormessagelength); 
			if (ierr > 0) { throw ValueError(format("%s",herr).c_str()); } else if (ierr < 0) {set_warning(format("%s",herr).c_str());}

			// Set all cache values that can be set
			break;
		}
	default:
	{
		throw ValueError(format("This set of inputs [%d,%d] is not yet supported",iName1,iName2));
	}
	throw NotImplementedError("REFPROPBackend is not yet implemented");
	}
}

} /* namespace CoolProp::Backends */
