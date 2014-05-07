import json as pyjson
from datetime import datetime
import struct

values = [
    ('all_fluids.json','../src/Fluids/all_fluids_JSON.h','all_fluids_JSON'),
    ('mixtures/mixture_excess_term.json', '../src/Backends/mixture_excess_term_JSON.h', 'mixture_excess_term_JSON'),
    ('mixtures/mixture_reducing_parameters.json', '../src/Backends/mixture_reducing_parameters_JSON.h', 'mixture_reducing_parameters_JSON')
]

def to_chunks(l, n):
    if n<1:
        n=1
    return [l[i:i+n] for i in range(0, len(l), n)]
        
for infile,outfile,variable in values:
    # Check you haven't messed up the JSON file and it will still load
    pyjson.loads(open(infile,'r').read())
    
    # see another idea:
    # http://stackoverflow.com/questions/7366391/embedding-a-text-file-in-an-exe-which-can-be-accessed-using-fopen
    
    json = open(infile,'r').read()
    
    # convert each character to hex and add a terminating NULL character to end the 
    # string, join into a comma separated string
    h = [hex(struct.unpack("b",b)[0]) for b in json] + ['0x00']
    
    # Break up the file into lines of 16 hex characters
    chunks = to_chunks(h, 16)
    
    # Put the lines back together again
    # The chunks are joined together with commas, and then EOL are used to join the rest
    hex_string = ',\n'.join([', '.join(chunk) for chunk in chunks])
    
    # Generate the output string
    output  = '// File generated by the script dev/JSON_to_C++.py on '+ str(datetime.now()) + '\n\n'
    output += '// JSON file encoded in binary form\n'
    output += 'const unsigned char '+variable+'_binary[] = {\n' + hex_string + '\n};'+'\n\n'
    output += '// Combined into a single std::string \n'
    output += 'std::string {v:s}({v:s}_binary, {v:s}_binary + sizeof({v:s}_binary)/sizeof({v:s}_binary[0]));'.format(v = variable)
    
    f = open(outfile,'w')
    f.write(output)
    f.close()
