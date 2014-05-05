
# Data from Mulero, JPCRD, 2012
# CAS codes added from CoolProp 4.2, and wikipedia where necessary

# raw_data = """Acetone 0.0633 1.160
# Ammonia 0.1028 1.211 -0.09453 5.585
# Argon 0.037 1.25
# Benzene 0.07298 1.232 -0.0007802 0.8635 -0.0001756 0.3065
# n-Butane 0.05138 1.209
# Butene 0.05644 1.248
# Dodecane 0.0154 4.180 0.0480 1.170
# Decafluorobutane 0.04429 1.242
# CarbonDioxide 0.07863 1.254
# CarbonMonoxide 0.02843 1.148
# CarbonylSulfide 0.07246 1.407
# Cyclohexane 0.06485 1.263
# Decane 0.05473 1.290
# Deuterium 0.009376 1.258
# D2O (Heavy water) -0.1423 2.645 0.2094 1.214
# DimethylEther 0.063157 1.2595
# Ethane 0.07602 1.320 -0.02912 1.676
# Ethanol 0.05 0.952
# Ethylene 0.0477 1.170
# Fluorine 0.03978 1.218
# Helium 0.0004656 1.040 0.001889 2.468 -0.002006 2.661
# Heptane 0.07765 1.319 -0.02599 1.600
# Hexane 0.210952 1.0962 -0.158485 1.05893
# Hydrogen -1.4165 0.63882 0.746383 0.659804 0.675625 0.619149
# HydrogenSulfide 0.078557 1.2074
# Isobutane -0.01639 2.102 0.06121 1.304
# Isobutene 0.0545 1.230
# Isohexane 0.05024 1.194
# Isopentane 0.051 1.209
# Krypton 0.0447 1.245
# Methane 0.03825 1.191 -0.006024 5.422 -0.0007065 0.6161
# Methanol 0.22421 1.3355 -0.21408 1.677 0.083233 4.4402
# Neon 0.012254 1.4136 0.02728 1.4517 -0.025715 1.6567
# Nitrogen 0.02898 1.246
# NitrousOxide 0.07087 1.204
# n-Nonane 0.05388 1.262
# Octane 0.34338 1.6607 -0.50634 1.9632 0.2238 2.3547
# Oxygen 0.03843 1.225
# Parahydrogen 0.005314 1.060
# Pentane 0.08015 1.408 0.004384 1.031 -0.03437 1.818
# Perfluoropentane 0.04394 1.254
# Propane 0.05334 1.235 -0.01748 4.404
# Propylene 0.05268 1.186
# Propyne 0.05801 1.205
# R11 0.06212 1.247
# R12 -0.000124 0.4318 0.05662 1.263
# R13 0.05045 1.269
# R14 0.0423 1.24
# R21 0.06924 1.259
# R22 3.0587 1.41809 -2.99856 1.42291
# R23 -0.32359 1.6055 0.37702 1.5232
# R32 0.07147 1.246
# R41 0.05049 1.242
# R113 0.0556 1.24
# R114 0.05239 1.258
# R115 0.04771 1.246
# R116 0.047593 1.2666 -0.0073402 1.9892
# R123 0.056151 1.2367
# R124 0.05175 1.197
# R125 0.05252 1.237
# R134a 0.05801 1.241
# R141b 7.3958e-5 0.066331 0.059941 1.2214
# R142b 0.05685 1.237
# R143a 0.05416 1.255
# R152a 0.05808 1.2115
# R161 0.05385 1.111
# R218 0.04322 1.224
# R227ea 0.06127 1.192 -0.009516 0.9795 -0.00192 1.421
# R236ea 0.306974 1.12614 -0.247277 1.09899
# R236fa 0.05389 1.249
# R245ca 0.069297 1.2795 -0.022419 3.1368
# R245fa 0.073586 1.0983 0.0103 0.60033 -0.02663 0.72765
# R365mfc 0.0534 1.210
# R1234yf 0.06274 1.394
# RC318 0.0507 1.250
# SulfurDioxide 0.0803 0.928 0.0139 1.570 -0.0114 0.364
# SulfurHexafluoride 0.0538 1.271 -4.064e-5 0.2116
# Toluene 0.06897 1.291
# Trifluoroiodomethane 0.05767 1.298
# Water -0.1306 2.471 0.2151 1.233
# Xenon -0.11538 1.0512 0.16598 1.098"""
# CAS_not_in_CP42 = {
#  'R245ca' : '679-86-7',
#  'Trifluoroiodomethane' : '2314-97-8',
#  'R115' :  '76-15-3',
#  'Perfluoropentane' : '678-26-2',
#  'Decafluorobutane' : '355-25-9',
#  'D2O' : '7789-20-0'
# }
# import CoolProp.CoolProp as CP    
# for line in raw_data.split('\n'):
#     name = line.split(' ')[0]
#     CAS = CP.get_fluid_param_string(name, 'CAS')
#     if CAS.find('-') == -1:
#         CAS = CAS_not_in_CP42[name]
#     print CAS + ' ' + line
    
    
Tc_dict = {'Argon':150.687,
'Benzene':562.02,
'1-Butene':419.29, #Butene from Mulero
'CarbonMonoxide':132.86,
'Cyclohexane':553.64,
'D2O':643.847,
'n-Decane':617.7,
'Ethane':305.322,
'Fluorine':144.414,
'Helium':5.1953,
'Isobutane':407.81,
'Isobutene':418.09,
'Isopentane':460.35,
'Methanol':513.38,
'Nitrogen':126.192,
'n-Nonane':594.55,
'Oxygen':154.581,
'Parahydrogen':32.938,
'Perfluorobutane':386.326,
'n-Propane':369.89,
'Propyne':402.38,
'R11':471.11,
'R113':487.21,
'R114':418.83,
'R115':353.1,
'R123':456.831,
'R1234yf':367.85,
'R124':395.425,
'R125':339.173,
'R134a':374.21,
'R14':227.51,
'R142b':410.26,
'R143a':345.857,
'R152a':386.411,
'R21':451.48,
'R218':345.02,
'R227ea':374.9,
'R32':351.255,
'R365mfc':460,
'RC318':388.38,
'SulfurDioxide':430.64,
'Toluene':591.75,
'Trifluoroiodomethane':396.44,
'Water':647.096,
'Acetone':508.1,
'n-Butane':425.125,
'CarbonDioxide':304.128,
'DimethylEther':400.378,
'n-Dodecane':658.1, #Decane in Mulero
'Ethylene':282.35,
'n-Heptane':540.13,
'n-Hexane':507.82,
'HydrogenSulfide':373.1,
'Isohexane':497.7,
'Krypton':209.48,
'NitrousOxide':309.52,
'n-Pentane':469.7,
'R116':293.03,
'R13':302,
'R23':299.293,
'R245ca':447.57,
'Xenon':289.733,
'CarbonylSulfide':378.77,
'Hydrogen':33.145,
'Methane':190.564,
'n-Octane':569.32,
'Perfluoropentane':420.555,
'Propylene':364.211,
'R12':385.12,
'R141b':477.5,
'R161':375.3,
'R22':369.295,
'R236ea':412.44,
'R236fa':398.07,
'R245fa':427.16,
'R41':317.28,
'SulfurHexafluoride':318.723,
'Ammonia':405.4,
'Deuterium':38.34,
'Ethanol':513.9,
'Neon':44.4918,
'Decafluorobutane':113.3+273.15 # According to http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=19#GeneralData, not in Mulero
}

Mulero_data = """67-64-1 Acetone 0.0633 1.160
7664-41-7 Ammonia 0.1028 1.211 -0.09453 5.585
7440-37-1 Argon 0.037 1.25
71-43-2 Benzene 0.07298 1.232 -0.0007802 0.8635 -0.0001756 0.3065
106-97-8 n-Butane 0.05138 1.209
106-98-9 1-Butene 0.05644 1.248
112-40-3 n-Dodecane 0.0154 4.180 0.0480 1.170
355-25-9 Decafluorobutane 0.04429 1.242
124-38-9 CarbonDioxide 0.07863 1.254
630-08-0 CarbonMonoxide 0.02843 1.148
463-58-1 CarbonylSulfide 0.07246 1.407
110-82-7 Cyclohexane 0.06485 1.263
124-18-5 n-Decane 0.05473 1.290
7782-39-0 Deuterium 0.009376 1.258
7789-20-0 D2O -0.1423 2.645 0.2094 1.214
115-10-6 DimethylEther 0.063157 1.2595
74-84-0 Ethane 0.07602 1.320 -0.02912 1.676
64-17-5 Ethanol 0.05 0.952
74-85-1 Ethylene 0.0477 1.170
7782-41-4 Fluorine 0.03978 1.218
7440-59-7 Helium 0.0004656 1.040 0.001889 2.468 -0.002006 2.661
142-82-5 n-Heptane 0.07765 1.319 -0.02599 1.600
110-54-3 n-Hexane 0.210952 1.0962 -0.158485 1.05893
1333-74-0 Hydrogen -1.4165 0.63882 0.746383 0.659804 0.675625 0.619149
7783-06-4 HydrogenSulfide 0.078557 1.2074
75-28-5 Isobutane -0.01639 2.102 0.06121 1.304
115-11-7 Isobutene 0.0545 1.230
107-83-5 Isohexane 0.05024 1.194
78-78-4 Isopentane 0.051 1.209
7439-90-9 Krypton 0.0447 1.245
74-82-8 Methane 0.03825 1.191 -0.006024 5.422 -0.0007065 0.6161
67-56-1 Methanol 0.22421 1.3355 -0.21408 1.677 0.083233 4.4402
7440-01-9 Neon 0.012254 1.4136 0.02728 1.4517 -0.025715 1.6567
7727-37-9 Nitrogen 0.02898 1.246
10024-97-2 NitrousOxide 0.07087 1.204
111-84-2 n-Nonane 0.05388 1.262
111-65-9 n-Octane 0.34338 1.6607 -0.50634 1.9632 0.2238 2.3547
7782-44-7 Oxygen 0.03843 1.225
1333-74-0p Parahydrogen 0.005314 1.060
109-66-0 n-Pentane 0.08015 1.408 0.004384 1.031 -0.03437 1.818
678-26-2 Perfluoropentane 0.04394 1.254
74-98-6 n-Propane 0.05334 1.235 -0.01748 4.404
115-07-1 Propylene 0.05268 1.186
74-99-7 Propyne 0.05801 1.205
75-69-4 R11 0.06212 1.247
75-71-8 R12 -0.000124 0.4318 0.05662 1.263
75-72-9 R13 0.05045 1.269
75-73-0 R14 0.0423 1.24
75-43-4 R21 0.06924 1.259
75-45-6 R22 3.0587 1.41809 -2.99856 1.42291
75-46-7 R23 -0.32359 1.6055 0.37702 1.5232
75-10-5 R32 0.07147 1.246
593-53-3 R41 0.05049 1.242
76-13-1 R113 0.0556 1.24
76-14-2 R114 0.05239 1.258
76-15-3 R115 0.04771 1.246
76-16-4 R116 0.047593 1.2666 -0.0073402 1.9892
306-83-2 R123 0.056151 1.2367
2837-89-0 R124 0.05175 1.197
354-33-6 R125 0.05252 1.237
811-97-2 R134a 0.05801 1.241
1717-00-6 R141b 7.3958e-5 0.066331 0.059941 1.2214
75-68-3 R142b 0.05685 1.237
420-46-2 R143a 0.05416 1.255
75-37-6 R152a 0.05808 1.2115
353-36-6 R161 0.05385 1.111
76-19-7 R218 0.04322 1.224
431-89-0 R227ea 0.06127 1.192 -0.009516 0.9795 -0.00192 1.421
431-63-0 R236ea 0.306974 1.12614 -0.247277 1.09899
690-39-1 R236fa 0.05389 1.249
679-86-7 R245ca 0.069297 1.2795 -0.022419 3.1368
460-73-1 R245fa 0.073586 1.0983 0.0103 0.60033 -0.02663 0.72765
406-58-6 R365mfc 0.0534 1.210
754-12-1 R1234yf 0.06274 1.394
115-25-3 RC318 0.0507 1.250
7446-09-5 SulfurDioxide 0.0803 0.928 0.0139 1.570 -0.0114 0.364
2551-62-4 SulfurHexafluoride 0.0538 1.271 -4.064e-5 0.2116
108-88-3 Toluene 0.06897 1.291
2314-97-8 Trifluoroiodomethane 0.05767 1.298
7732-18-5 Water -0.1306 2.471 0.2151 1.233
7440-63-3 Xenon -0.11538 1.0512 0.16598 1.098"""

import glob, numpy as np, json, os
if __name__=='__main__':
    for row in Mulero_data.split('\n'):
        row = row.split(' ')
        cas = row.pop(0)
        name = row.pop(0)
        a = row[0:len(row):2]
        a = [float(_) for _ in a]
        n = row[1:len(row):2]
        n = [float(_) for _ in n]
        if name not in Tc_dict:
            raise ValueError( 'could not find Tc for ' + name)
            continue
        Tc = Tc_dict[name]
        
        # The dictionary of values for the surface tension
        j_st = dict(Tc = Tc, 
                    a = np.array(a).tolist(), 
                    n = n,
                    BibTeX = 'Mulero-JPCRD-2012',
                    description = 'sigma = sum(a_i*(1-T/Tc)^n_i)'
                    )
        
        fname = 'fluids/'+name+'.json'
        if not os.path.exists(fname):
            print fname+' does not exist'
            continue
            
        j = json.load(open(fname,'r'))
        
        j['ANCILLARIES']['surface_tension'] = j_st
        
        fp = open(fname, 'w')
        fp.write(json.dumps(j, indent = 2, sort_keys = True))
        fp.close()