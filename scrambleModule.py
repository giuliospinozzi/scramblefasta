#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import random
from numpy import ceil

#########################################################################################
####### GLOBAL VARS
#########################################################################################

nucleotides = list("ACGTN")

#########################################################################################
### FUNCTIONS
#########################################################################################

def doDeletion(input_string, del_type, del_parameters):
    """
    deletions: <del_type> | <del_parameters>:
        from_last | len >0
        from_first | len >0
        Nbp | [N size {1..seqlen}, starting position 0-based >0 && <len(input_string)]
    """
    instring_list = list(input_string)
    output_string = ""
    if del_type == "from_last":
        output_string = ''.join( instring_list[:len(instring_list)-del_parameters] )
    elif del_type == "from_first":
        output_string = ''.join( instring_list[del_parameters:] )
    elif del_type == "Nbp":
        output_string = ''.join( instring_list[:del_parameters[1]] + instring_list[del_parameters[1]+del_parameters[0]:] )
    else:
        print "[AP]\tError, deletion type not defined."
        sys.exit()
    if ''.join(instring_list) == output_string:
        print "it's a bad story, man!"
    return output_string

def doInsertion(input_string, ins_type, ins_parameters, random_seq = True):
    """
    if random(ACGTN) do:
        insertions: <ins_type> | <ins_parameters>:
            from_last | len >0 
            from_first | len >0
            Nbp | [N size {1..seqlen}, starting position 0-based >0 && <len(input_string)]
    if defined dstring do:
        insertions: <ins_type> | <ins_parameters>:
            from_Nbp (not used because so far it is the only one) | [user defined string, starting position 0-based >0 && <len(input_string)]
    """
    instring_list = list(input_string)
    output_string = ""
    if random_seq:
        if ins_type == "from_last":
            output_string = ''.join(instring_list) + ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters))
        elif ins_type == "from_first":
            output_string = ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters)) + ''.join(instring_list)
        elif ins_type == "Nbp":
            output_string = ''.join( instring_list[:ins_parameters[1]] ) + ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters[0])) + ''.join(instring_list[ins_parameters[1]:] )
        else:
            print "[AP]\tError, insertion type not defined."
            sys.exit()
    else: # if not random
        output_string = ''.join( instring_list[:ins_parameters[1]] ) + ins_parameters[0] + ''.join(instring_list[ins_parameters[1]:])
    if ''.join(instring_list) == output_string:
        print "it's a bad story, man!"
    return output_string

def doMutation(input_string, mut_type, mut_parameters):
    """
    (random) mutation excluding same string to replace: <mut_type> | <mut_parameters>:
        Nbp | [start bp, end bp 0-based] where this is a closed interval thus you are mutating from the staring the the ending bp included.

    NB: each base in the interval MUST be different, else you could obtain an unexpected mutation in 2 positions in the same span: from AAAA, span 2 from start, you may obtain CCAA but also ACAA -> the first one is ok but the second one corresponds to the case of SINGLE mutation at base 2.
    """
    instring_list = list(input_string)
    output_string = ''.join( instring_list[:mut_parameters[0]] )
    if mut_type == "Nbp":
        if mut_parameters[1]-mut_parameters[0]<=0:
            print "[AP]\tWarning: your mutated string has an input interval NON positive! mut_parameters[1] - mut_parameters[0] <= 0::", mut_parameters[1], mut_parameters[0]
        else:
            for bp_index in range(mut_parameters[0],mut_parameters[1]):
                char_tomutate = str(instring_list[bp_index]).capitalize()
                other_nucleotides = [x for x in nucleotides if x != char_tomutate]
                random_string = str(random.choice(other_nucleotides))
                output_string += random_string
            output_string += ''.join(instring_list[mut_parameters[1]:])
    else:
        print "[AP]\tError, mutation type not defined."
        sys.exit()
    if ''.join(instring_list) == output_string:
        print "it's a bad story, man!"
    return output_string



def getArgs(func_name, input_string, task_type):
    """
    This functions aims to provide the input of the functions above,
    defining some rules to assign xxx_parameters
    """
    
    string_len = len(input_string)
    
    ### Settings
    
    # from_first (len) limit
    from_first_limit = min(4, string_len-1)
    
    # from_last (len) limit
    from_last_limit = min(4, string_len-1)
    
    # Nbp
    N_bp_lim_inf = 0
    N_bp_lim_sup = string_len -1
    N_bp_span_min = 1
    N_bp_span_max = 1 # this way each function acts on one nucleotide only!
    #                 # Here an example, if you want a N_bp_span_max related to the length of the sequence: N_bp_span_max = int(ceil(string_len/100.0))
    
    
    ### Assignment
    
    if task_type == 'from_first':
        return input_string, task_type, random.randint(1, from_first_limit)
        
    elif task_type == 'from_last':
        return input_string, task_type, random.randint(1, from_last_limit)
        
    elif task_type == 'Nbp':
        
        span = random.randint(N_bp_span_min, N_bp_span_max)
        start = None
        
        if (func_name == 'doDeletion') or (func_name == 'doMutation'):
            start = random.randint(N_bp_lim_inf, N_bp_lim_sup - span)
            if (func_name == 'doDeletion'):
                return input_string, task_type, [span, start]
            elif (func_name == 'doMutation'):
                return input_string, task_type, [start, start+span]
        
        elif (func_name == 'doInsertion'):    
            start = random.randint(N_bp_lim_inf, N_bp_lim_sup)
            return input_string, task_type, [span, start]
