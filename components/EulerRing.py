# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:51:35 2023

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

from M321_Orion.components.VarEuler180 import VarEuler180
from M321_Orion.components.CmptMetalTaper import CmptMetalTaper

Vcc = CmptMetalTaper(length=4.5,width1=1,width2=2,layer=LAYER.M2AM)
GND = CmptMetalTaper(length=4.5,width1=1,width2=2,layer=LAYER.M1AM)

@gf.cell
def EulerRing(radius: float = 4, 
              Lc: float = 0, 
              width1: float = 0.48, 
              width2: float = 1, 
              slab: float = 1.5, 
              etch: float = 0.5, 
              p: float = 0.5,
              doping: bool = True, 
              ribRef: bool = False) -> gf.Component:
    c = gf.Component()
    
    s0 = Section(width=width1, offset=0, layer=(0,1), port_names=('o1','o2'))
    s1 = Section(width=width1+2*slab, offset=0, layer=(0,2))
    s2 = Section(width=width1+2*etch, offset=0, layer=(0,3))
    print(radius)
    if Lc == 0:
        bend1 = c << VarEuler180(radius = radius, 
                                 width1 = width1, 
                                 width2 = width2, slab = slab, etch = etch, p=p, doping=doping, metal='M1', ribRef=ribRef)
        
        bend2 = c << VarEuler180(radius = radius, 
                                 width1 = width1, 
                                 width2 = width2, slab = slab, etch = etch, p=p, doping=doping, metal='M2', ribRef=ribRef)
        bend2.mirror_x()
        bend1.connect('o1',bend2.ports['o1'])
        
    else:
        bend1 = c << VarEuler180(radius = radius, 
                                 width1 = width1, 
                                 width2 = width2, slab = slab, etch = etch, p=p, doping=doping, ribRef=ribRef)
        
        bend2 = c << VarEuler180(radius = radius, 
                                 width1 = width1, 
                                 width2 = width2, slab = slab, etch = etch, p=p, doping=doping, ribRef=ribRef)
        bend2.mirror_x()
        ctmp = gf.Component('tmp')
        xs = CrossSection(sections=[s0,s1,s2])
        straight1 = ctmp << gf.path.straight(Lc).extrude(xs)
        straight2 = ctmp << gf.path.straight(Lc).extrude(xs)
        straight1.connect('o2',bend1.ports['o1'])
        straight2.connect('o2',bend1.ports['o2'])
        bend2.connect('o1',straight1.ports['o1'])
        c << gf.geometry.boolean([ctmp.extract(layers=[(0,2),])], [], 'A+B', layer=LAYER.SEAM)
        c << gf.geometry.boolean([ctmp.extract(layers=[(0,3),])], [ctmp.extract(layers=[(0,1),])], 'A-B', layer=LAYER.REAM)
        if ribRef:
            c << gf.geometry.boolean([ctmp.extract(layers=[(0,1)])], [], 'A+B', layer=LAYER.Waveguide)
    
    c.add_ports(bend1.get_ports_list())
    if doping:
        
        Vcc = c << CmptMetalTaper(length=4.5,width1=1,width2=2,layer=LAYER.M2AM)
        Vcc.move([bend1.xmax+5-Vcc.xmin, bend1.ports['gnd'].y+2-Vcc.ports['e1'].y])
        Vgnd = c << CmptMetalTaper(length=4.5,width1=1,width2=2,layer=LAYER.M1AM)
        Vgnd.mirror_y()
        Vgnd.move([bend1.xmax+5-Vgnd.xmin, bend1.ports['gnd'].y-Vgnd.ports['e1'].y])
        
        routeList = list()
        connections = (
            [[bend1.ports['gnd']],Vgnd.ports['e1']],
            [[bend2.ports['vcc']],Vcc.ports['e1']],
            )
        stepList = [
            None,
            [{"dy":int(radius)},{"dx":3*int(radius)}],
            ]
        xsList = [M1AMXS, M2AMXS]
        layerList = [LAYER.M1AM,LAYER.M2AM]
        for N, connection in enumerate(connections):
            routeList.append(
                gf.routing.get_bundle(
                    connection[0], connection[1],
                    cross_section=xsList[N],
                    layer=layerList[N],
                    bend=wire_corner,
                    width=1,
                    steps=stepList[N]
                )
            )
    
        for routes in routeList:
            for route in routes:
                c.add(route.references)
    
        c.add_port(name='00_heater+',port=Vgnd.ports['e2'],port_type='electrical')
        c.add_port(name='00_heater-',port=Vcc.ports['e2'],port_type='electrical')
    
    return c

if __name__ == "__main__":
    EulerRing(radius = 3.5, width1=0.48, width2=0.6 ,etch = 2, slab=3, p=0.75, doping=False, ribRef=False).show(show_ports=True)