#! /usr/bin/env python

"""db.py

Usage:
    db.py split FILE PATH
    db.py concat PATH FILE...
    db.py display FILE...
    db.py explicit FILE...
    db.py implicit FILE...
    db.py dirk FILE...

Options:
    -h, --help    show help
    FILE          json file to split or join
    PATH          path to the output directory (for split) of output file (for concat)
"""

import json
import os, sys

def hash_label(label):
    id = label.lower()
    replacements = [
                (" ","_"),
                ("(",""),
                (")",""),
                (",",""),
                ("-",""),
                ("/",""),
            ]
    for old,new in replacements :
        id = id.replace(old,new)
    return id

def split(filename,outputdir="."):
    """ split json array in filename into multiple files

    name of new files build from the label of each element
    if no label, use a hash of each element as a name
    """
    with open( filename ,'r') as f:
        data = json.load(f)
    for elm in data:
        basename = hash_label(elm['label']) if 'label' in elm else hex(abs(hash(str(elm))))
        json.dump(
            elm,
            open(os.path.join(outputdir, basename+".json") ,'w'),
            indent=4
        )

def concat(filenames,outputfile="data.json"):
    """ concatenate all filesnames into one json file
    """
    data = []
    for filename in filenames:
        data.append(json.load(open(filename,'r')))
    json.dump(data,open(outputfile,'w'),indent=4)

def display(filename):
    """ display Butcher tableau
    """
    import sympy as sp
    with open( filename ,'r') as f:
        data = json.load(f)
    A = sp.Matrix(data['A'])
    b = sp.Matrix(data['b'])
    c = sp.Matrix(data['c'])

    lenCol = [0]*(1+len(A.tolist()))
    lenCol[0] = max([ len(str(ci)) for ci in c ])
    for i in range(A.shape[1]):
        lenCol[i+1] = max([ len(str(aij)) for aij in A.col(i) ])
    for i,bi in enumerate(b):
        lenCol[i+1] = max(lenCol[i+1],len(str(bi)))

    print(data['label'])

    for i in range(A.shape[1]):
        print("{{:^{}}}".format(lenCol[0]).format(str(c[i])), end=" │ ")
        print(" ".join([ "{{:^{}}}".format(lenCol[j+1]).format(str(aij)) for j,aij in enumerate(A.row(i)) ]))
    print("─"*lenCol[0],end="─┼─")
    print("─"*(sum(lenCol[1:])+len(lenCol[1:])))
    print(" "*lenCol[0],end=" │ ")
    print(" ".join([ "{{:^{}}}".format(lenCol[j+1]).format(str(bj)) for j,bj in enumerate(b) ]))

def is_explicit(A):
    """ test if a matrix A of a Butcher tableau represent an explicit RK method
    """
    import sympy as sp
    return ( sum([
        sum([ sp.Abs(A[i,j]) for j in range(i,A.shape[0]) ])
        for i in range(0,A.shape[1])
        ]) == 0 )

def is_implicit(A):
    """ test if a matrix A of a Butcher tableau represent an implicit RK method
    """
    return not is_explicit(A)

def is_dirk(A):
    """ test if a matrix A of a Butcher tableau represent a DIRK method
    """
    import sympy as sp
    return ( sum([
            sum([ sp.Abs(A[i,j]) for j in range(i+1,A.shape[0]) ])
            for i in range(0,A.shape[1])
            ]) == 0 ) and ( sum([ sp.Abs(A[i,i]) for i in range(0,A.shape[0]) ]) != 0 )

def extract_A(filename):
    """ extract matrix A in Butcher tableau in filename
    """
    import sympy as sp
    with open(filename,'r') as f:
        data = json.load(f)
    return sp.Matrix(data['A'])

def filter_rk(fn, filenames):
    """ filter in filename with fn function
    """
    r = filter( lambda f:fn(extract_A(f)), filenames )
    print(*list(r),sep=" ")

def explicit(filenames):
    filter_rk(is_explicit,filenames)

def implicit(filenames):
    filter_rk(is_implicit,filenames)

def dirk(filenames):
    filter_rk(is_dirk,filenames)

if __name__ == '__main__':
    from docopt import docopt

    args = docopt( __doc__ , sys.argv[1:] )

    if args['concat'] :
        concat(args['FILE'],args['PATH'])
    elif args['split'] :
        split(args['FILE'][0],args['PATH'])
    elif args['display'] :
        for file in args['FILE']:
            display(file)
            print()
    elif args['explicit'] :
        explicit(args['FILE'])
    elif args['implicit'] :
        implicit(args['FILE'])
    elif args['dirk'] :
        dirk(args['FILE'])
    else:
        print("error...")
