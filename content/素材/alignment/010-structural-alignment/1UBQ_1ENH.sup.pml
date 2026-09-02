#!/usr/bin/env pymol
cmd.load("1UBQ_1ENH.sup.pdb", "structure1")
cmd.load("pdbs/1ENH.pdb", "structure2")
hide all
set all_states, off
remove not n. CA and not n. C3'
bond structure1 and i.    7, structure1 and i.    8
bond structure1 and i.    8, structure1 and i.    9
bond structure1 and i.    9, structure1 and i.   10
bond structure1 and i.   10, structure1 and i.   11
bond structure1 and i.   11, structure1 and i.   12
bond structure1 and i.   12, structure1 and i.   13
bond structure1 and i.   13, structure1 and i.   14
bond structure1 and i.   14, structure1 and i.   15
bond structure1 and i.   15, structure1 and i.   19
bond structure1 and i.   19, structure1 and i.   21
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
bond structure1 and i.   31, structure1 and i.   32
bond structure1 and i.   32, structure1 and i.   33
bond structure1 and i.   33, structure1 and i.   34
bond structure1 and i.   34, structure1 and i.   35
bond structure1 and i.   35, structure1 and i.   36
bond structure1 and i.   36, structure1 and i.   37
bond structure1 and i.   37, structure1 and i.   38
bond structure1 and i.   38, structure1 and i.   39
bond structure1 and i.   39, structure1 and i.   40
bond structure1 and i.   40, structure1 and i.   41
bond structure1 and i.   41, structure1 and i.   42
bond structure1 and i.   42, structure1 and i.   43
bond structure1 and i.   43, structure1 and i.   44
bond structure1 and i.   44, structure1 and i.   45
bond structure1 and i.   45, structure1 and i.   46
bond structure1 and i.   46, structure1 and i.   47
bond structure1 and i.   47, structure1 and i.   48
bond structure2 and i.    7, structure2 and i.    8
bond structure2 and i.    8, structure2 and i.    9
bond structure2 and i.    9, structure2 and i.   10
bond structure2 and i.   10, structure2 and i.   11
bond structure2 and i.   11, structure2 and i.   14
bond structure2 and i.   14, structure2 and i.   15
bond structure2 and i.   15, structure2 and i.   18
bond structure2 and i.   18, structure2 and i.   19
bond structure2 and i.   19, structure2 and i.   23
bond structure2 and i.   23, structure2 and i.   25
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
bond structure2 and i.   35, structure2 and i.   36
bond structure2 and i.   36, structure2 and i.   37
bond structure2 and i.   37, structure2 and i.   38
bond structure2 and i.   38, structure2 and i.   39
bond structure2 and i.   39, structure2 and i.   40
bond structure2 and i.   40, structure2 and i.   41
bond structure2 and i.   41, structure2 and i.   42
bond structure2 and i.   42, structure2 and i.   43
bond structure2 and i.   43, structure2 and i.   44
bond structure2 and i.   44, structure2 and i.   45
bond structure2 and i.   45, structure2 and i.   47
bond structure2 and i.   47, structure2 and i.   49
bond structure2 and i.   49, structure2 and i.   50
bond structure2 and i.   50, structure2 and i.   52
bond structure2 and i.   52, structure2 and i.   53
bond structure2 and i.   53, structure2 and i.   54
bond structure2 and i.   54, structure2 and i.   55
show stick, structure1 and ( i.    7 or i.    8 or i.    9 or i.   10 or i.   11 or i.   12 or i.   13 or i.   14 or i.   15 or i.   19 or i.   21 or i.   22 or i.   23 or i.   24 or i.   25 or i.   26 or i.   27 or i.   28 or i.   29 or i.   30 or i.   31 or i.   32 or i.   33 or i.   34 or i.   35 or i.   36 or i.   37 or i.   38 or i.   39 or i.   40 or i.   41 or i.   42 or i.   43 or i.   44 or i.   45 or i.   46 or i.   47 or i.   48)
show stick, structure2 and ( i.    7 or i.    8 or i.    9 or i.   10 or i.   11 or i.   14 or i.   15 or i.   18 or i.   19 or i.   23 or i.   25 or i.   26 or i.   27 or i.   28 or i.   29 or i.   30 or i.   31 or i.   32 or i.   33 or i.   34 or i.   35 or i.   36 or i.   37 or i.   38 or i.   39 or i.   40 or i.   41 or i.   42 or i.   43 or i.   44 or i.   45 or i.   47 or i.   49 or i.   50 or i.   52 or i.   53 or i.   54 or i.   55)
color blue, structure1
color red, structure2
set ribbon_width, 6
set stick_radius, 0.3
set sphere_scale, 0.25
set ray_shadow, 0
bg_color white
set transparency=0.2
zoom polymer and ((structure1) or (structure2))

