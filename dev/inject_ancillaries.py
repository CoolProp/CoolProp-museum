import os, json, glob

def inject_ancillaries():
    master = []

    for file in glob.glob(os.path.join('fluids','*.json')):
        path, file_name = os.path.split(file)
        fluid_name = file_name.split('.')[0]
        # Load the fluid file
        fluid = json.load(open(os.path.join('fluids', fluid_name+'.json'), 'r'))
        
        # Load the ancillary
        anc = json.load(open(os.path.join('ancillaries',fluid_name+'_anc.json'),'r'))
        # Apply the ancillary by merging dictionaries
        fluid.update(anc)
        # Write fluid back to file
        fp = open(os.path.join('fluids', fluid_name+'.json'),'w')
        fp.write(json.dumps(fluid, indent = 2))
    
if __name__=='__main__':
    inject_ancillaries()