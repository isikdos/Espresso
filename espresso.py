#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 22:08:31 2019

@author: iaricanli
"""

import copy

T = True
F = False
D = "_"

"""
Generate the inputs for the algorithm. A list of dictionarys.
Each element of the list represents a different boolean expression input --
say if we are trying to reduce something for a 7 segment display over a period
of multiple timesteps.

Each dictionary represents a boolean function. The key is the boolean input
represented in integer form (aka A^!B^!C^D -> 1001 -> 9) and the value is 
either True, False, or Don't Care.

Reading from a file not yet supported.

Do whatever you want in here.

Arguments:
    LEN : integer
        - The dimensionality of the desired output truth table. AKA, 
          the number of boolean variables.
          
Return:
    tt (truth table): list (dict (int, Booleanish))
        - The only question here really is what is a booleanish?
          in addition to True and False there is a third concept of Don'tCare
          which is being represented here as "_". It fails "is True" but 
          passes "== True". this is abused heavily.
"""
def Get_Truth_Table(LEN):
    tt = list()
    e = dict()

    e = {
            0: T,
            1: T,
            2: D,   
            3: T,
            4: F,
            5: F,
            6: T,
            7: T,
            }

    tt.append(e)
    return tt

"""
This is a helper function, existant in case I want to
expand the code and allow different iteration over the inputs.
If this returned randomly, the code would no longer necessarily
output predictable results. 
"""
def _Provide_Index(dim):
    return range(0, dim)

"""
Performs the expansion of a cube through the Nspace
Attempting to expand through ever dimension, one at a time.

While it does this, it maps the boolean expressions to minterms
and it maps the minterms to the boolean expressions. Thus providing
the program a quick view of the rows and columns regarding the results found.

Arguments:
    boolean array: dict (int, truth)
        - The key maps to the integer representation of the inputs
        - The value points to whether that is mapped to True, False, or DC
    
    idx: int
        - The space in the boolean array we are beginning at, where the
            expansion begins from.
            
    dim: int
        - The total number of dimensions we are operating in.
       
    # REFERENCED BY VALUE
    minterms2bln: dict (int, boolean_expression)
        - Maps the minterms (tracked by integer -- same as in the idx we see
            above), and keeps track of which minterms are related to what
            boolean expressions.
            
Return:
    covered_minterms: set
        - The defined set of minterms covered by the boolean expression
    bln_expr: boolean_expression
        - The boolean expression (defined by a point and then a mask)
            that covered the aforementioned minterms.

"""
def Expand_Cube(boolean_array, idx, dim, minterms2bln):
    bln_expr = boolean_expression(idx, 0)
    
    # Define the space of the cube
    space = [idx]
    covered_minterms = {idx}
    if idx in minterms2bln:
        minterms2bln[idx].add(bln_expr)
    else:
        minterms2bln[idx] = {bln_expr}

    # Iterate over the indices however we decide
    for i in _Provide_Index(dim): 
        
        # Control variable to exit a loop
        _continue = False
        
        # Convert the index into the representitive integer
        dim2int = 2**i
        
        # The space being explored
        new_space = list()
        for index in space:
            # MAGIC LINE
            # We need to turn 1s into 0s and 0s into 1s, depending on the index
            new_index = index ^ dim2int

            # We're expanding the cube,  verify that we're expanding it into
            # valid space. If it is valid, add the expanding indices into list
            if new_index in boolean_array and boolean_array[new_index]:               
                new_space.append(new_index)
            else:
                # If the new space doesn't pan out _perfectly_, keep going to
                # the the next index
                _continue = True
                break

        # We don't want to extend into the space of the selected index
        # if it didn't pan out. So skip this one and move on to the next 
        # dimension.
        if not _continue:
            # We like the new dimension, and are going to cover all the new 
            # elements into it.

            space.extend(new_space)
            
            for ns in new_space:
                # If the value at the boolean array is specifically
                # True and not just Don't Care, add it to the Covered Minterms
                if boolean_array[ns] is T:
                    covered_minterms.add(ns)
                    if ns in minterms2bln:
                        minterms2bln[ns].add(bln_expr)
                    else:
                        minterms2bln[ns] = {bln_expr}
            
            # Allow the Mask to contain the information regarding the dimension
            # that was just covered.
            bln_expr.mask += dim2int
     
    return covered_minterms, bln_expr

class boolean_expression(object):
    def __init__(self, idx, mask):
        self.idx = idx
        self.mask = mask

    def __eq__(self, b):
        return self.idx == b.idx and self.mask == b.mask
    
    def __hash__(self):
        return hash((self.idx, self.mask))

    def __str__(self):
        return "boolean_expression({0}, {1})".format(self.idx, self.mask)
    
    def __repr__(self):
        return self.__str__()

def Expand(truth_table, dim):
    #
    # Iterate over every boolean output
    # 
    
    expr_per_output = list()
    
    for boolean_array in truth_table:
        bln2minterms= dict()
        minterms2bln = dict()        

        for idx, bln in boolean_array.items():
            if bln is T:
                covered_minterms, bln_expr = Expand_Cube(boolean_array, 
                                                         idx, 
                                                         dim,
                                                         minterms2bln)   
                bln2minterms[bln_expr] = covered_minterms
        
        # bln2minterms and minterms2bln
        # are two dictionaries that are dually referent
        # in order to keep computations fast.
        
        expr_per_output.append((bln2minterms, minterms2bln))
        
    return expr_per_output
    
def Intersect( list_of_maps ):
    #
    # Finds intersections between boolean statements and
    # the minterms they cover
    #
    lom = list()

    # Iterate over every solution-set per output
    itr_list_of_maps = copy.deepcopy(list_of_maps)
    for bln2minterms, minterms2bln in itr_list_of_maps:
        # First we're going to look for any case where a minterm
        # maps to only one boolean expression.
        required_blns = set()
        todelete = list()
        itr_minterms2bln = copy.deepcopy(minterms2bln)
        for minterm, set_of_blns in itr_minterms2bln.items():

            if len(set_of_blns) == 1:
                # WE found one!
                # Take it
                required_bln = set_of_blns.pop()

                # Now find all the minterms related to the boolean
                minterms_correlated_to_bln = bln2minterms[required_bln]
                # Iterate over them
                for correlated_minterm in minterms_correlated_to_bln:
                    # and remove the boolean from their knowledge
                    minterms2bln[correlated_minterm].remove(required_bln)

                    
                # Then delete the entire boolean from the booly-books
                del bln2minterms[required_bln]
                todelete.append(minterm)
                # And remember what we've done on this day, this evil day.
                required_blns.add( required_bln )
     
        for i in todelete:
            del minterms2bln[i]
        
        # Now we get rid of booleans as we determine that they are "the best candidate
        while len(minterms2bln):
            # We are looking at only a SINGLE minterm.
            # Scanning a subspace to decrease overall computation time
            # and keep everything in linear time.
            minterm = Select_Minterm(minterms2bln)
            most = 0
            best_candidate = None
            # We determine the "Best candidate" as the boolean expression
            # with the greatest number of related minterms
            for bln in minterms2bln[minterm]:

                if len(bln2minterms[bln]) > most:
                    best_candidate = bln
                    most = len(bln2minterms[bln])
            
            required_blns.add( best_candidate )
            # Now find all the minterms related to the boolean
            minterms_correlated_to_bln = bln2minterms[best_candidate]
            # Iterate over them
            todelete = list()
            for correlated_minterm in minterms_correlated_to_bln:
                # Delete all minterms correlated to the highest-scoring boolean
                for related_bln in minterms2bln[correlated_minterm]:
                    todelete.append((related_bln, correlated_minterm))
                # Forreal, delete them
                del minterms2bln[correlated_minterm]
            
            for related_bln, correlated_minterm in todelete:
                bln2minterms[related_bln].remove(correlated_minterm)

            # The ndelete the aforementioned best candidate
            del bln2minterms[best_candidate]
    
                    
        lom.append(required_blns)
        
    return lom
            
"""
This is a helper function, existant in case I want to
expand the code and allow different iteration over the inputs.
If this returned randomly, the code would no longer necessarily
output predictable results. 
"""
def Select_Minterm(minterms2bln):
    return list(minterms2bln.keys())[0]
            

def main(dim):

    #
    # Define truth table
    #
    
    truth_table = Get_Truth_Table(dim)

    # Perform the Expand operation on every output set
    list_of_maps = Expand(truth_table, dim) 
    list_of_covering_blns = Intersect(list_of_maps)

    return list_of_covering_blns

