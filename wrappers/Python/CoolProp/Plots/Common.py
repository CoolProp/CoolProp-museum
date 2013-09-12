# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:13:20 2013

@author: logan
"""

import matplotlib


class BasePlot(object):
    AXIS_LABLES = {'T': ["Temperature", r"[$K$]"],
                   'P': ["Pressure", r"[$kPa$]"],
                   'S': ["Entropy", r"[$kJ/kg K$]"],
                   'H': ["Enthalpy", r"[$kJ/kg$]"],
                   'V': [],
                   'RHO': ["Density", r"[$kg/m^3$]"]}

    COLOR_MAP = {'T': 'Darkred',
                 'P': 'DarkCyan',
                 'H': 'DarkGreen',
                 'D': 'DarkBlue',
                 'S': 'DarkOrange',
                 'Q': 'black',}

    LINE_IDS = {'TS': ['P', 'D'], #'H'],
                'PH': ['S', 'T', 'D'],
                'HS': ['P'], #'T', 'D'],
                'PS': ['H', 'T', 'D'],
                'PD': ['T', 'S', 'H'],
                'TD': ['P'], #'S', 'H'],
                'PT': ['D', 'P', 'S'],}

    def __init__(self, fluid_ref, graph_type, **kwargs):
        if not isinstance(graph_type, str):
            raise TypeError("Invalid graph_type input, expeceted a string")

        graph_type = graph_type.upper()
        if len(graph_type) >= 2 and graph_type[1:len(graph_type)] == 'RHO':
            graph_type = graph_type[0] + graph_type[1:len(graph_type)]

        if graph_type.upper() not in self.LINE_IDS.keys():
            raise ValueError("You have to specify the kind of plot, use " \
                              + str(self.LINE_IDS.keys()))

        self.fluid_ref = fluid_ref
        self.graph_type = graph_type.upper()

        self.figure = kwargs.get('fig', matplotlib.pyplot.figure())
        self.axis = kwargs.get('axis', matplotlib.pyplot.gca())

    def __get_sat_lines(kind=None, kmin=None, kmax=None, x=[0.,1.]):
        """
        Calculates bubble and dew line in the quantities for your plot.

        You can specify if you need evenly spaced entries in either
        pressure or temperature by supplying kind='p' and kind='T'
        (default), respectively.

        Limits can be set with kmin (default: minimum from EOS) and
        kmax (default: critical value).

        Returns lines[] - a 2D array of dicts containing 'x' and 'y' coordinates
        for bubble and dew line. Additionally, the dict holds the keys
        'kmax', 'label' and 'opts', those can be used for plotting as well.
        """

        xName,yName,plot = _plotXY(plot)

        if (str(kind).lower()=='p'):
            kind = 'P'
        else:
            kind = 'T'

        (kmin,kmax) = _satBounds(Ref, kind, xmin=kmin, xmax=kmax)
        k0          = numpy.linspace(kmin,kmax,1000)

        iName       = 'Q'
        iVal        = x
        kVal        = [k0 for i in iVal]

        if xName!=kind:
            (Xx,Yx) = _getI_YX(Ref,iName,kind,xName,iVal,kVal)
        else:
            (Xx,Yx) = (kVal,kVal)

        if yName!=kind:
            (Xy,Yy) = _getI_YX(Ref,iName,kind,yName,iVal,kVal)
        else:
            (Xy,Yy) = (kVal,kVal)

        # Merge the two lines, capital Y holds important information. We merge on X values
        # Every entry, eg. Xy, contains two arrays of values.
        lines = []
        for j in range(len(Yx)): # two dimensions: i = {0,1}
            line = {
              'x' : Yx[j],
              'y' : Yy[j],
              'kmax' : kmax
              }
            if iVal[j]==0.:
                line['label'] = 'bubble line'
                line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1.00 }
            elif iVal[j]==1.:
                line['label'] = 'dew line'
                line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1.00 }
            else:
                line['label'] = _getIsoLineLabel(iName,iVal[j]),
                line['opts'] = { 'color':getIsoLineColour(iName), 'lw':0.75, 'alpha':0.5}

            lines.append(line)

        return lines