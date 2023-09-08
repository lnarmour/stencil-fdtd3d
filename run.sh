#!/bin/bash

function check_dep {
    dep=$1
    hard_req=$2
    if [[ -n $hard_req ]]; then lvl='error'; else lvl='warning'; fi
    which $dep > /dev/null 2>&1
    if [[ $? != 0 ]]; then
        echo "$lvl: '$dep' does not appear to be installed"
        if [[ -n $hard_req ]]; then 
            exit 1
        fi
    fi
}

function check_py_dep {
    dep=$1
    hard_req=$2
    if [[ -n $hard_req ]]; then lvl='error'; else lvl='warning'; fi
    python3 -c "import $dep" > /dev/null 2>&1
    if [[ $? != 0 ]]; then
        echo "$lvl: python package '$dep' does not appear to be installed"
        if [[ -n $hard_req ]]; then 
            exit 1
        fi
    fi
}

# check args
if [[ $# != 7 ]]; then
    echo "usage: $0 T Nx Ny Nz Sx Sy Sz"
    echo ""
    echo "Parameters:"
    echo "    T        : number of time steps"
    echo "    Nx,Ny,Nz : data grid size"
    echo "    Sx,Sy,Sz : coordinates of energy source (0<=Si<Ni)"
    echo ""
    exit 1
fi

# check dependencies
check_dep gcc 1
check_dep make 1
check_dep python3 1
check_py_dep numpy 1
check_py_dep matplotlib 1
check_dep ffmpeg

# compile C code
make -C codegen

# run stencil
./codegen/fdtd3d $@

# generate images (and mp4 if ffmpeg installed)
python3 codegen/plot.py $@

if [[ -f sim.mp4 ]]; then
    echo "generated sim.mp4"
fi
