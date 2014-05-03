
# Fluids from Reeves, 1964, "Melting curves of Pressure-Transmitting Fluids"
Simon_curves = """
{
    "n-Propane" : {
        "T_0" : 85.3, "a" : 7.180e8, "c" : 1.283, "p_0" : 0.0, "T_max" : 168.63, "BibTeX" : "Reeves-JCP-1964"
    },
    "n-Pentane" : {
        "T_0" : 143.5, "a" : 6.600e8, "c" : 1.649, "p_0" : 0.0, "T_max" : 156.2, "BibTeX" : "Reeves-JCP-1964"
    },
    "Isopentane" : {
        "T_0" : 112.5, "a" : 5.916e8, "c" : 1.563, "p_0" : 0, "T_max" : 212.16, "BibTeX" : "Reeves-JCP-1964"
    },
    "Propylene" : [
        {
            "T_0" : 86.0, "a" : 3.196e8, "c" : 2.821, "p_0" : 0, "T_max" : 109.6, "BibTeX" : "Reeves-JCP-1964"
        },
        {
            "T_0" : 109.6, "a" : 3.064e8, "c" : 3.871, "p_0" : 4.450e8, "T_max" : 145.3, "BibTeX" : "Reeves-JCP-1964"
        }
    ]
}
"""

polynomial_in_Tr = """
{
    "Argon" : {
        "T_t" : 83.8058, "a" : [-7476.2665, 9959.0613], "t" : [1.05,1.275], "p_t" : 68891, "T_max" : 254.006, "BibTeX" : "Tegeler-JPCRD-1999"
    },
    "Fluorine" : {
        "T_t" : 53.4811, "a" : [988043.478261], "t" : [2.1845], "p_t" : 252, "T_max" : 55.3811
    },
    "Nitrogen" : {
        "T_t" : 63.151, "a" : [12798.61], "t" : [1.78963], "p_t" : 12523, "T_max" : 283.751
    },
    "Ethane" : {
        "T_t" : 90.368, "a" : [2.23626315e8, 1.05262374e8], "t" : [1.0, 2.55], "p_t" : 1.14, "T_max" : 110.17, "BibTeX" : "Buecker-JCRD-2006"
    },
    "Isobutane" : {
        "T_t" : 113.73, "a" : [1.9536371309e9], "t" : [6.12], "p_t" : 0.0219, "T_max" : 124.93, "BibTeX" : "Buecker-JPCRD-2006B"
    },
    "Ethylene" : [
    {
        "T_t" : 90.368, "a" : [2947001.84], "t" : [2.045], "p_t" : 122.65, "T_max" : 103.989
    },
    {
        "T_t" : 103.989, "a" : [6.82693421], "t" : [1.089], "p_t" : 46.8e6, "T_max" : 110.369
    }
    ],
    "n-Butane" : {
        "T_t" : 134.895, "a" : [5.585582364e8], "t" : [2.206], "p_t" : 0.653, "T_max" : 163.9, "BibTeX" : "Buecker-JPCRD-2006B"
    }
}
"""

polynomial_in_theta = """
{
    "Methanol" : {
        "T_t" : 175.61, "a" : [5.330770e9, 4.524780e9, 3.888861e10], "t" : [1, 1.5, 4], "p_t" : 0.187, "T_max" : 245.91
    },
    "CarbonDioxide" : {
        "T_t" : 216.592, "a" : [1955.5390, 2055.4593], "t" : [1, 2], "p_t" : 51795, "T_max" : 327.6
    }
}
"""

import json, numpy as np, matplotlib.pyplot as plt, pandas
    
# for fluid,value in (json.loads(Simon_curves)).iteritems():
#     if isinstance(value, list):
#         print fluid
#         continue
#     Tmin = value['T_0']
#     Tmax = value['T_max']
#     T = np.linspace(Tmin, Tmax, 200)
#     T_0 = value['T_0']
#     p_0 = value['p_0']
#     a = value['a']
#     c = value['c']
#     
#     p = p_0 + a*(T-T_0)**c - 1
#     
#     plt.plot(T, p)
#     
#     df = pandas.read_csv('melting_curves/'+fluid+'.mlt',names=['T','p','rho'])
#     
#     plt.plot(df['T'], df['p'], 'o')
#     
#     plt.show()
    

for fluid, value in (json.loads(polynomial_in_Tr)).iteritems():
    if isinstance(value, list):
        print fluid
        continue
    print fluid
    a = value['a']
    t = value['t']
    T_t = value['T_t']
    p_t = value['p_t']
    
    Tmin = T_t
    Tmax = value['T_max']
    T = np.linspace(Tmin, Tmax, 200)
    
    RHS = 0
    for i in range(len(a)):
        RHS += a[i]*((T/T_t)**t[i] - 1)
    
    p = p_t*(RHS + 1)
    
    plt.plot(T, p)
    
    df = pandas.read_csv('melting_curves/' + fluid + '.mlt', names=['T','p','rho'])
    
    plt.plot(df['T'], df['p'], 'o')
    
    plt.show()
    
# print json.loads(polynomial_in_theta)

# plt.show()
    
#CycloHexane.mlt 401.7