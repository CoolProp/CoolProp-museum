
import subprocess, os
import glob2 as glob

exports = ['-s','EXPORTED_FUNCTIONS=\"[\'_main\',\'_F2K\',\'_HAProps\',\'_Props1\',\'_PropsS\']\"']
optimization = '-O2'

def compile_sources():
    for f in glob.glob(os.path.join('..','..','CoolProp-v5','**','*.cpp')):
        
        # Don't compile CoolPropDLL in the compile pass to avoid duplicate symbols
        if f.find('..\..\CoolProp-v5\main.cpp') > -1: 
            continue 
        
        call = ['em++',optimization,f,'-I../../CoolProp-v5','-c','-DEXTERNC']+ exports
        print 'Calling:',' '.join(call)
        subprocess.check_output(' '.join(call), shell = True)

def link():
    call = ['em++',optimization,'-o','coolprop.js']+glob.glob('*.o')+['-DEXTERNC']  +  exports
    print 'Calling:',' '.join(call)
    subprocess.check_output(' '.join(call), shell = True)

def closure_compiler():
    call = ['java','-Xmx1024m','-jar','compiler.jar','--js','coolprop.js','--js_output_file','coolprop2.js','--compilation_level','ADVANCED_OPTIMIZATIONS','--language_in','ECMASCRIPT5']
    print 'Using the closure compiler, this will take a while...   (from https://developers.google.com/closure/compiler/)'
    print 'Calling:',' '.join(call)
    subprocess.check_output(' '.join(call), shell = True)

def cleanup():
    for file in glob.glob('*.o'):
        print 'removing',file
        os.remove(file)

def run():
    os.startfile('index.html')

if __name__=='__main__':
    compile_sources()
    link()
    #closure_compiler()
    cleanup()
    #run()