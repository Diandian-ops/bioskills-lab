#!/usr/bin/env pymol
cmd.load("1UBQ_1CRN.sup.pdb", "structure1")
cmd.load("pdbs/1CRN.pdb", "structure2")
hide all
set all_states, off
remove not n. CA and not n. C3'
bond structure1 and i.    5, structure1 and i.    6
bond structure1 and i.    6, structure1 and i.    7
bond structure1 and i.    7, structure1 and i.   11
bond structure1 and i.   11, structure1 and i.   12
bond structure1 and i.   12, structure1 and i.   13
bond structure1 and i.   13, structure1 and i.   14
bond structure1 and i.   14, structure1 and i.   15
bond structure1 and i.   15, structure1 and i.   16
bond structure1 and i.   16, structure1 and i.   17
bond structure1 and i.   17, structure1 and i.   18
bond structure1 and i.   18, structure1 and i.   19
bond structure1 and i.   19, structure1 and i.   20
bond structure1 and i.   20, structure1 and i.   21
bond structure1 and i.   21, structure1 and i.   22
bond structure1 and i.   22, structure1 and i.   23
bond structure1 and i.   23, structure1 and i.   24
bond structure1 and i.   24, structure1 and i.   25
bond structure1 and i.   25, structure1 and i.   26
bond structure1 and i.   26, structure1 and i.   27
bond structure1 and i.   27, structure1 and i.   28
bond structure1 and i.   28, structure1 and i.   29
bond structure1 and i.   29, structure1 and i.   30
bond structure1 and i.   30, structure1 and i.   31
bond structure1 and i.   31, structure1 and i.   36
bond structure1 and i.   36, structure1 and i.   41
bond structure1 and i.   41, structure1 and i.   42
bond structure1 and i.   42, structure1 and i.   43
bond structure1 and i.   43, structure1 and i.   48
bond structure1 and i.   48, structure1 and i.   49
bond structure1 and i.   49, structure1 and i.   70
bond structure1 and i.   70, structure1 and i.   71
bond structure2 and i.    3, structure2 and i.    4
bond structure2 and i.    4, structure2 and i.    5
bond structure2 and i.    5, structure2 and i.    6
bond structure2 and i.    6, structure2 and i.    7
bond structure2 and i.    7, structure2 and i.    8
bond structure2 and i.    8, structure2 and i.   11
bond structure2 and i.   11, structure2 and i.   12
bond structure2 and i.   12, structure2 and i.   15
bond structure2 and i.   15, structure2 and i.   16
bond structure2 and i.   16, structure2 and i.   18
bond structure2 and i.   18, structure2 and i.   19
bond structure2 and i.   19, structure2 and i.   20
bond structure2 and i.   20, structure2 and i.   21
bond structure2 and i.   21, structure2 and i.   22
bond structure2 and i.   22, structure2 and i.   23
bond structure2 and i.   23, structure2 and i.   24
bond structure2 and i.   24, structure2 and i.   25
bond structure2 and i.   25, structure2 and i.   26
bond structure2 and i.   26, structure2 and i.   27
bond structure2 and i.   27, structure2 and i.   28
bond structure2 and i.   28, structure2 and i.   29
bond structure2 and i.   29, structure2 and i.   30
bond structure2 and i.   30, structure2 and i.   31
bond structure2 and i.   31, structure2 and i.   32
bond structure2 and i.   32, structure2 and i.   33
bond structure2 and i.   33, structure2 and i.   34
bond structure2 and i.   34, structure2 and i.   35
bond structure2 and i.   35, structure2 and i.   37
bond structure2 and i.   37, structure2 and i.   38
bond structure2 and i.   38, structure2 and i.   40
bond structure2 and i.   40, structure2 and i.   41
show stick, structure1 and ( i.    5 or i.    6 or i.    7 or i.   11 or i.   12 or i.   13 or i.   14 or i.   15 or i.   16 or i.   17 or i.   18 or i.   19 or i.   20 or i.   21 or i.   22 or i.   23 or i.   24 or i.   25 or i.   26 or i.   27 or i.   28 or i.   29 or i.   30 or i.   31 or i.   36 or i.   41 or i.   42 or i.   43 or i.   48 or i.   49 or i.   70 or i.   71)
show stick, structure2 and ( i.    3 or i.    4 or i.    5 or i.    6 or i.    7 or i.    8 or i.   11 or i.   12 or i.   15 or i.   16 or i.   18 or i.   19 or i.   20 or i.   21 or i.   22 or i.   23 or i.   24 or i.   25 or i.   26 or i.   27 or i.   28 or i.   29 or i.   30 or i.   31 or i.   32 or i.   33 or i.   34 or i.   35 or i.   37 or i.   38 or i.   40 or i.   41)
color blue, structure1
color red, structure2
set ribbon_width, 6
set stick_radius, 0.3
set sphere_scale, 0.25
set ray_shadow, 0
bg_color white
set transparency=0.2
zoom polymer and ((structure1) or (structure2))

