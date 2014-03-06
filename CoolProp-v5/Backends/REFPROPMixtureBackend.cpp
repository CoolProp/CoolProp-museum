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
//#include "CoolProp.h"

#if defined(__ISWINDOWS__)
#include <windows.h>
#elif defined(__ISLINUX__)
#include <dlfcn.h>
#elif defined(__ISAPPLE__)
#include <dlfcn.h>
#endif

#include "REFPROP_lib.h"
#include "REFPROP.h"
#include "CoolPropTools.h"

#include <stdlib.h>
#include "string.h"
#include <stdio.h>
#include <iostream>

// Some constants for REFPROP... defined by macros for ease of use 
#define refpropcharlength 255
#define filepathlength 255
#define lengthofreference 3
#define errormessagelength 255
#define ncmax 20		// Note: ncmax is the max number of components
#define numparams 72 
#define maxcoefs 50

// Check windows
#if _WIN32 || _WIN64
   #if _WIN64
     #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

// Check GCC
#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

#include "REFPROPBackend.h"

#if defined(__ISWINDOWS__)
HINSTANCE RefpropdllInstance=NULL;
#elif defined(__ISLINUX__)
void *RefpropdllInstance=NULL;
#elif defined(__ISAPPLE__)
void *RefpropdllInstance=NULL;
#else
void *RefpropdllInstance=NULL;
#endif

bool load_REFPROP()
{
	// If REFPROP is not loaded
	if (RefpropdllInstance==NULL)
	{
		// Load it
		#if defined(__ISWINDOWS__)
			#if defined(ENV64BIT)
				// 64-bit code here.
				TCHAR refpropdllstring[100] = TEXT("refprp64.dll");
				RefpropdllInstance = LoadLibrary(refpropdllstring);
			#elif defined (ENV32BIT)
				// 32-bit code here.
				TCHAR refpropdllstring[100] = TEXT("refprop.dll");
				RefpropdllInstance = LoadLibrary(refpropdllstring);
			#else
				// INCREASE ROBUSTNESS. ALWAYS THROW AN ERROR ON THE ELSE.
				#error "Must define either ENV32BIT or ENV64BIT"
			#endif

			
		#elif defined(__ISLINUX__)
			RefpropdllInstance = dlopen ("librefprop.so", RTLD_LAZY);
		#elif defined(__ISAPPLE__)
			RefpropdllInstance = dlopen ("librefprop.dylib", RTLD_LAZY);
		#else
			throw NotImplementedError("We should not reach this point.");
			RefpropdllInstance = NULL;
		#endif

		if (RefpropdllInstance==NULL)
		{
			#if defined(__ISWINDOWS__)
//				int  dw            = ::GetLastError();
//				char lpBuffer[256] = _T("?");
//				if(dwLastError != 0) {
//				    ::FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM,       // Has to be a system error
//				                     NULL,                            // No formatter
//				                     dw,                              // Get error message for this int
//				                     MAKELANGID(LANG_NEUTRAL,SUBLANG_DEFAULT),  // Use system language
//				                     lpBuffer,                        // Write output
//				                     STR_ELEMS(lpBuffer)-1,           // Length of output
//				                     NULL);
//				}
//				printf(lpBuffer);
//				printf("\n");
				              printf("Could not load refprop.dll \n\n");
				throw AttributeError("Could not load refprop.dll, make sure it is in your system search path. In case you run 64bit and you have a REFPROP license, try installing the 64bit DLL from NIST.");
			#elif defined(__ISLINUX__)
				fputs (dlerror(), stderr);
				              printf("Could not load librefprop.so \n\n");
				throw AttributeError("Could not load librefprop.so, make sure it is in your system search path.");
			#elif defined(__ISAPPLE__)
				fputs (dlerror(), stderr);
				              printf("Could not load librefprop.dylib \n\n");
				throw AttributeError("Could not load librefprop.dylib, make sure it is in your system search path.");
			#else
				throw NotImplementedError("Something is wrong with the platform definition, you should not end up here.");
			#endif
			return false;
		}

		#if defined(__ISWINDOWS__)
		
		// Get data associated with path using the windows libraries, 
		// and if you can (result == 0), the path exists
		#ifdef __MINGW32__
			struct stat buf;
			if ( stat( "c:\\Program Files\\REFPROP\\fluids", &buf) != 0){
				throw ValueError("REFPROP fluid files must be copied to c:\\Program Files\\REFPROP\\fluids");
			}
		#else
			struct _stat buf;
			if ( _stat( "c:\\Program Files\\REFPROP\\fluids", &buf) != 0){
				throw ValueError("REFPROP fluid files must be copied to c:\\Program Files\\REFPROP\\fluids");
			}
		#endif
		#endif

		if (setFunctionPointers()!=COOLPROP_OK)
		{
			              printf("There was an error setting the REFPROP function pointers, check types and names in header file.\n");
			throw AttributeError("There was an error setting the REFPROP function pointers, check types and names in header file.");
			return false;
		}
		return true;
	}
	return true;
}

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
