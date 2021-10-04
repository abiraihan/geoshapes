# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 18:20:37 2021

@author: ABIR RAIHAN
"""

class ShapeError(Exception):
    
    def __init__(self, geoms:str):
        
        self.geoms = geoms
        self.message = f"\n{self.geoms} is not acceptable.\nOnly 'Square', 'Rectangle', 'Square-Rectangle' are acceptable geometry shape to process"
        super().__init__(self.message)

class RangeError(Exception):
    
    def __init__(self, distance, reason):
        
        self.distance = distance
        self.reason = reason
        self.message = f"\nReason -- {self.reason} : Out of range. Assigned distance should be at least less than {self.distance} to identify a zones"
        super().__init__(self.message)

class InputError(Exception):
    
    def __init__(self, conditions, reasons):
        
        self.conditions = conditions
        self.reasons = reasons
        self.message = f"\nReason -- {self.reasons} : |'{self.conditions}'| isn't acceptable. Check Input parameters"
        super().__init__(self.message)