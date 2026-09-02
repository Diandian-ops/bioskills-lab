#!/usr/bin/env pymol
cmd.load("1UBQ_1CRN.sup.pdb", "structure1")
cmd.load("pdbs/1CRN.pdb", "structure2")
hide all
set all_states, off
show cartoon, structure1 and ( i.    5 or i.    6 or i.    7 or i.   11 or i.   12 or i.   13 or i.   14 or i.   15 or i.   16 or i.   17 or i.   18 or i.   19 or i.   20 or i.   21 or i.   22 or i.   23 or i.   24 or i.   25 or i.   26 or i.   27 or i.   28 or i.   29 or i.   30 or i.   31 or i.   36 or i.   41 or i.   42 or i.   43 or i.   48 or i.   49 or i.   70 or i.   71)
show cartoon, structure2 and ( i.    3 or i.    4 or i.    5 or i.    6 or i.    7 or i.    8 or i.   11 or i.   12 or i.   15 or i.   16 or i.   18 or i.   19 or i.   20 or i.   21 or i.   22 or i.   23 or i.   24 or i.   25 or i.   26 or i.   27 or i.   28 or i.   29 or i.   30 or i.   31 or i.   32 or i.   33 or i.   34 or i.   35 or i.   37 or i.   38 or i.   40 or i.   41)
color blue, structure1
color red, structure2
set ribbon_width, 6
set stick_radius, 0.3
set sphere_scale, 0.25
set ray_shadow, 0
bg_color white
set transparency=0.2
zoom polymer and ((structure1) or (structure2))

