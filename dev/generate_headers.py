"""
In this module, we do some of the preparatory work that is needed to get
CoolProp ready to build.  This includes setting the correct versions in the 
headers, generating the fluid files, etc.
"""
from __future__ import division, print_function
from datetime import datetime
import subprocess
import os
import sys

def version_to_file(root_dir):
    print('*** Generating version.h ***')
    # Get the version from the version.txt file
    version = open(os.path.join(root_dir,'version.txt'),'r').read().strip()
    
    # Format the string to be written
    string_for_file = '//Generated by the preprocess.py script on {t:s}\n\nstatic char version [] ="{v:s}";'.format(t = str(datetime.now()),v = version)
    
    # Include path relative to the root
    include_dir = os.path.join(root_dir,'include')
    
    # The name of the file to be written into
    file_name = os.path.join(include_dir,'version.h')
    
    # Write to file
    f = open(file_name,'w')
    f.write(string_for_file)
    f.close()
    
def gitrev_to_file(root_dir):
    """
    If a git repo, use git to update the gitrevision.  If not a git repo, read 
    the gitrevision from the gitrevision.txt file.  Otherwise, fail.
    """
    print('*** Generating gitrevision.h ***')
        
    try:
        try:
            subprocess.check_call('git --version', shell=True)
            print('git is accessible at the command line')
        except subprocess.CalledProcessError:
            print('git was not found')
            return
        p = subprocess.Popen('git rev-parse HEAD', 
                             stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE,
                             shell = True)
        stdout, stderr = p.communicate()

        # Include path relative to the root
        include_dir = os.path.join(root_dir,'include')
        
        if p.returncode != 0:
            print('tried to update git revision, but it failed for some reason (building from zip file?)')
            gitstring = '//Generated by the preprocess.py script on {t:s}\n\n std::string gitrevision = \"???????\";'
            f = open(os.path.join(include_dir,'gitrevision.h'),'w')
            f.write(gitstring)
            f.close()
            return
        else:
            rev = stdout.strip()
            print('git revision is', rev)
            
            gitstring = '//Generated by the preprocess.py script on {t:s}\n\nstd::string gitrevision = \"{rev:s}\";'.format(t = str(datetime.now()), rev = rev)
            
            try:
                is_hash = rev.find(' ') == -1 # python 2.x
            except TypeError:
                is_hash = ' ' not in str(rev) # python 3.x
                                
            if not is_hash:
                raise ValueError('No hash returned from call to git, got '+rev+' instead')
            
            f = open(os.path.join(include_dir,'gitrevision.h'),'w')
            f.write(gitstring)
            f.close()
                
    except (subprocess.CalledProcessError,OSError):
        pass
    
if __name__=='__main__':
    version_to_file(root_dir = '..')
    gitrev_to_file(root_dir = '..')
    import JSON_to_CPP
    JSON_to_CPP.TO_CPP(root_dir = '..')