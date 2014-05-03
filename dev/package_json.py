import os, json, glob

def package_json():

    master = []
    
    for file in glob.glob(os.path.join('fluids','*.json')):
        
        path, file_name = os.path.split(file)
        fluid_name = file_name.split('.')[0]
        
        # Load the fluid file
        fluid = json.load(open(os.path.join('fluids', fluid_name+'.json'), 'r'))
        
        master += [fluid]

    fp = open('all_fluids_verbose.json','w')
    fp.write(json.dumps(master, indent = 2))
    fp.close()
    
    fp = open('all_fluids.json','w')
    fp.write(json.dumps(master))
    fp.close()
    
if __name__=='__main__':
    package_json()