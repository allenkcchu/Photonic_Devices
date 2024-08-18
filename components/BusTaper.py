# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 17:01:14 2023

@author: achu
"""

import sys
sys.path.append("C:\\Users\\achu\\Documents\\PIC\\7268 - NGC AOX P2\\aim_photonics_pdk_v6p5_act-master\\EPDA_Assisted_Tapeout\\GDSfactory")
sys.path.append("C:\\Users\\achu\\Documents\\PIC\\7268 - NGC AOX P2\\Layout")
import numpy as np
import matplotlib.pyplot as plt

import gdsfactory as gf
from gdsfactory.typings import LayerSpec, CrossSectionSpec
from gdsfactory.cross_section import cross_section, Transition, CrossSection, Section
from gdsfactory.components.wire import wire_corner
from AIMPhotonics_ACT1.tech import LAYER, SEstrip, M1AMXS, M2AMXS
import AIMPhotonics_ACT1 as AIM

@gf.cell
def BusTaper(length = 7, Lc = 0, ribIn = 0.48, ribMid = 0.4, slab = 0.7, ribRef = False):
    c = gf.Component()
    
    rib_in = Section(width=ribIn, offset=0, layer=(0,1), port_names=('o1','o2'))
    etch_in = Section(width=ribIn+2*(0+0.5), offset=0, layer=(0,2))
    slab_in = Section(width=ribIn+2*0, offset=0, layer=(0,3))
    rib_mid = Section(width=ribMid, offset=0, layer=(0,1), port_names=('o1','o2'))
    etch_mid = Section(width=ribMid+2*(slab+0.5), offset=0, layer=(0,2))
    slab_mid = Section(width=ribMid+2*slab, offset=0, layer=(0,3))
    
    x_in = CrossSection(sections=[rib_in, etch_in, slab_in])
    x_mid = CrossSection(sections=[rib_mid, etch_mid, slab_mid])
    
    xs_in = Transition(cross_section1 = x_in, 
                       cross_section2 = x_mid, 
                       width_type = 'linear')
    xs_out = Transition(cross_section1 = x_mid, 
                        cross_section2 = x_in, 
                        width_type = 'linear')
    
    ctmp = gf.Component('tmp')
    buffer_in = c << gf.path.straight(1).extrude(SEstrip)
    taper_in = ctmp << gf.path.straight(length).extrude(xs_in)
    CouplingL = ctmp << gf.path.straight(Lc).extrude(x_mid)
    taper_out = ctmp << gf.path.straight(length).extrude(xs_out)
    buffer_out = c << gf.path.straight(1).extrude(SEstrip)
    
    taper_in.connect('o1', buffer_in.ports['o2'])
    CouplingL.connect('o1', taper_in.ports['o2'])
    taper_out.connect('o1', CouplingL.ports['o2'])
    buffer_out.connect('o1', taper_out.ports['o2'])
    c.add_port(name='BusLoc',port=taper_out.ports['o1'])
    
    c << gf.geometry.boolean([ctmp.extract(layers=[(0,3)])], [], 'A+B', layer=LAYER.SEAM)
    c << gf.geometry.boolean([ctmp.extract(layers=[(0,2)])], [ctmp.extract(layers=[(0,1)])], 'A-B', layer=LAYER.REAM)
    if ribRef:
        c << gf.geometry.boolean([ctmp.extract(layers=[(0,1)])], [], 'A+B', layer=LAYER.Waveguide)
    
    c.add_port(name='o1',port=buffer_in.ports['o1'])
    c.add_port(name='o2',port=buffer_out.ports['o2'])
    
    return c

if __name__ == "__main__":
    BusTaper(Lc = 0, ribRef=True).show(show_ports=True)