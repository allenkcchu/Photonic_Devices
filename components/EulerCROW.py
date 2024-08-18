# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:17:13 2023

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
from gdsfactory.pdk import get_active_pdk
from AIMPhotonics_ACT1.tech import LAYER, SEstrip, M1AMXS, M2AMXS
import AIMPhotonics_ACT1 as AIM

from Photonic_Devices.components.EulerRing import EulerRing
from Photonic_Devices.components.BusTaper import BusTaper

@gf.cell
def EulerCROW(order = 2, radius = 4, Lc = 0, gap1 = 0.15, gap2 = 0.15, bus_width = 0.4, ring_width = 0.5, 
              slab = 1.5, etch = 0.7, p=0.5, length = 15, doping=True):
    c = gf.Component()
    WGwidth = 0.48
    
    bus = BusTaper(length = length,
                   Lc = Lc,ribIn=WGwidth,ribMid=bus_width,slab=etch, ribRef=True)
    ring = EulerRing(radius = radius, Lc=Lc,
                     width1 = WGwidth, 
                     width2 = ring_width, slab = slab, etch = etch, p=p, doping=doping, ribRef=True)
    
    ringtmp = gf.Component('ringtmp')
    bustmp = gf.Component('bustmp')
    
    if order == 1:
        bus1 = bustmp << bus
        bus2 = bustmp << bus
        ring1 = ringtmp << ring
    
        bus1.move([ring1.ports['o1'].x-bus1.ports['BusLoc'].x,
                   ring1.ports['o1'].y-WGwidth/2-gap1-bus_width/2-bus1.ports['BusLoc'].y])
        bus2.move([ring1.ports['o2'].x-bus2.ports['BusLoc'].x,
                   ring1.ports['o2'].y+WGwidth/2+gap1+bus_width/2-bus2.ports['BusLoc'].y])
        
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.REAM]),
                                  bustmp.extract(layers=[LAYER.REAM])], 
                                 [ringtmp.extract(layers=[LAYER.Waveguide]),
                                  bustmp.extract(layers=[LAYER.Waveguide])], 'A-B', layer=LAYER.REAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.SEAM]),
                                  bustmp.extract(layers=[LAYER.SEAM])],
                                 [], 'A+B', layer=LAYER.SEAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.NNNAM]),],
                                 [], 'A+B', layer=LAYER.NNNAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.CBAM]),],
                                 [], 'A+B', layer=LAYER.CBAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.M1AM]),],
                                 [], 'A+B', layer=LAYER.M1AM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.M2AM]),],
                                 [], 'A+B', layer=LAYER.M2AM)
        
        c.add_port('00_input',port = bus1.ports['o1'])
        c.add_port('01_through',port = bus1.ports['o2'])
        c.add_port('03_drop',port = bus2.ports['o1'])
        c.add_port('02_add',port = bus2.ports['o2'])
        c.add_ports(ring1.get_ports_list(port_type='electrical',width=2),suffix='_ring1')
    elif order == 2:
        bus1 = bustmp << bus
        bus2 = bustmp << bus
        ring1 = ringtmp << ring
        ring2 = ringtmp << ring
        
        ring2.move([ring1.ports['o2'].x-ring2.ports['o1'].x,
                    ring1.ports['o2'].y+WGwidth+gap2-ring2.ports['o1'].y])
    
        bus1.move([ring1.ports['o1'].x-bus1.ports['BusLoc'].x,
                   ring1.ports['o1'].y-WGwidth/2-gap1-bus_width/2-bus1.ports['BusLoc'].y])
        bus2.move([ring2.ports['o2'].x-bus2.ports['BusLoc'].x,
                   ring2.ports['o2'].y+WGwidth/2+gap1+bus_width/2-bus2.ports['BusLoc'].y])
        
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.REAM]),
                                  bustmp.extract(layers=[LAYER.REAM])], 
                                 [ringtmp.extract(layers=[LAYER.Waveguide]),
                                  bustmp.extract(layers=[LAYER.Waveguide])], 'A-B', layer=LAYER.REAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.SEAM]),
                                  bustmp.extract(layers=[LAYER.SEAM])],
                                 [], 'A+B', layer=LAYER.SEAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.NNNAM]),],
                                 [], 'A+B', layer=LAYER.NNNAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.CBAM]),],
                                 [], 'A+B', layer=LAYER.CBAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.M1AM]),],
                                 [], 'A+B', layer=LAYER.M1AM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.M2AM]),],
                                 [], 'A+B', layer=LAYER.M2AM)
        
        c.add_port('00_input',port = bus1.ports['o1'])
        c.add_port('01_through',port = bus1.ports['o2'])
        c.add_port('03_drop',port = bus2.ports['o2'])
        c.add_port('02_add',port = bus2.ports['o1'])
        c.add_ports(ring1.get_ports_list(port_type='electrical',width=2),suffix='_ring1')
        c.add_ports(ring2.get_ports_list(port_type='electrical',width=2),suffix='_ring2')
    elif order == 3:
        bus1 = bustmp << bus
        bus2 = bustmp << bus
        ring1 = ringtmp << ring
        ring2 = ringtmp << ring
        ring3 = ringtmp << ring
        
        ring2.move([ring1.ports['o2'].x-ring2.ports['o1'].x,
                    ring1.ports['o2'].y+WGwidth+gap2-ring2.ports['o1'].y])
        ring3.move([ring2.ports['o2'].x-ring3.ports['o1'].x,
                    ring2.ports['o2'].y+WGwidth+gap2-ring3.ports['o1'].y])
    
        bus1.move([ring1.ports['o1'].x-bus1.ports['BusLoc'].x,
                   ring1.ports['o1'].y-WGwidth/2-gap1-bus_width/2-bus1.ports['BusLoc'].y])
        bus2.move([ring3.ports['o2'].x-bus2.ports['BusLoc'].x,
                   ring3.ports['o2'].y+WGwidth/2+gap1+bus_width/2-bus2.ports['BusLoc'].y])
        
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.REAM]),
                                  bustmp.extract(layers=[LAYER.REAM])], 
                                 [ringtmp.extract(layers=[LAYER.Waveguide]),
                                  bustmp.extract(layers=[LAYER.Waveguide])], 'A-B', layer=LAYER.REAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.SEAM]),
                                  bustmp.extract(layers=[LAYER.SEAM])],
                                 [], 'A+B', layer=LAYER.SEAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.NNNAM]),],
                                 [], 'A+B', layer=LAYER.NNNAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.CBAM]),],
                                 [], 'A+B', layer=LAYER.CBAM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.M1AM]),],
                                 [], 'A+B', layer=LAYER.M1AM)
        c << gf.geometry.boolean([ringtmp.extract(layers=[LAYER.M2AM]),],
                                 [], 'A+B', layer=LAYER.M2AM)
        
        c.add_port('00_input',port = bus1.ports['o1'])
        c.add_port('01_through',port = bus1.ports['o2'])
        c.add_port('03_drop',port = bus2.ports['o1'])
        c.add_port('02_add',port = bus2.ports['o2'])
        c.add_ports(ring1.get_ports_list(port_type='electrical',width=2),suffix='_ring1')
        c.add_ports(ring2.get_ports_list(port_type='electrical',width=2),suffix='_ring2')
        c.add_ports(ring3.get_ports_list(port_type='electrical',width=2),suffix='_ring3')
    return c

if __name__ == "__main__":
    EulerCROW(order = 3, radius = 4.5, slab = 1.5, etch = 0.7, Lc = 0, bus_width = 0.48, ring_width=0.6, p=0.75,
              gap1 = 0.25, gap2 = 0.4, doping=False).show(show_ports=True)
    
    
    # Generate bend for simulation
    
    # orderList = [3]
    # pList = [0.75] 
    # # rList = np.arange(2,6.1,0.1)
    # rList = [3.5]
    # # LcList = np.arange(0,2.1,0.2)
    # LcList = [0]
    # gap1List = [0.25]
    # gap2List = [0.4]
    # # gap2List = [0.4]
    # busWidth = 0.48
    # ringWidthList = [0.6]
    # length = 30
    
    # for order in orderList:
    #     for p in pList:
    #         for radius in rList:
    #             for ringWidth in ringWidthList:
    #                 for Lc in LcList:
    #                     for gap1 in gap1List:
    #                         for gap2 in gap2List:
    #                             gf.clear_cache()
    #                             c = gf.Component(f"crow{order}_p{int(p*100)}_r{int(radius*1000)}_Lc{int(Lc*1000)}_busW{int(busWidth*1000)}_gapbr{int(gap1*1e3)}_ringW{int(ringWidth*1000)}_gaprr{int(gap2*1e3)}")
    #                             # crow = c << EulerCROW(order = order,
    #                             #                       radius = radius, 
    #                             #                       Lc = Lc, 
    #                             #                       gap1 = gap1, 
    #                             #                       gap2 = gap2, 
    #                             #                       bus_width = busWidth, 
    #                             #                       ring_width = ringWidth, 
    #                             #                       slab = 1.5, 
    #                             #                       etch = 0.7, 
    #                             #                       length = length,
    #                             #                       p=p, doping=False)
    #                             crow = c << EulerCROW(order = order,
    #                                                   radius = radius, 
    #                                                   Lc = Lc, 
    #                                                   gap1 = gap1, 
    #                                                   gap2 = gap2, 
    #                                                   bus_width = busWidth, 
    #                                                   ring_width = ringWidth, 
    #                                                   slab = 6, 
    #                                                   etch = 5, 
    #                                                   length = length,
    #                                                   p=p, doping=False)
                                
    #                             PortRef = c << gf.components.bbox(bbox=([c.xmin,crow.ports['00_input'].y],
    #                                                                     [c.xmax,crow.ports['03_drop'].y]),layer=(0,0))
                                
    #                             c.show()
    #                             import os
    #                             # SimPath = 'Z:\\USERS\\Allen\\Design\\CROW\\Simulation_Structures\\'
    #                             # SimPath = 'C:\\Users\\achu\\Documents\\PIC\\7268 - NGC AOX P2\Layout\\Simulation\\'
    #                             c.write_gds(os.path.join(SimPath,f"crow{order}_XLSlab_p{int(p*100)}_r{int(radius*1000)}_Lc{int(Lc*1000)}_busW{int(busWidth*1000)}_gapbr{int(gap1*1e3)}_ringW{int(ringWidth*1000)}_gaprr{int(gap2*1e3)}.gds"))