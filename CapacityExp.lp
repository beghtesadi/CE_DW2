\ SCIP STATISTICS
\   Problem name     : t_vaccine
\   Variables        : 31 (0 binary, 0 integer, 0 implicit integer, 31 continuous)
\   Constraints      : 23
\   Obj. scale       : 1
\   Obj. offset      : 0
Minimize
 Obj: +1.76197314026085 t_wt2_0 +2.23802685973915 t_wt2_1 +1.76197314026085 t_wt3_0 +2.23802685973915 t_wt3_1
      +1.76197314026085 t_wt3_2 +2.23802685973915 t_wt3_3 +55.0616606331516 Alpha2_0 +69.9383393668485 Alpha2_1
      +55.0616606331516 Alpha3_0 +19.5827350227176 Alpha3_1 +55.0616606331516 Alpha3_2 +19.5827350227176 Alpha3_3
      +33.0369963798909 Alpha2_0 +41.9630036201091 Alpha2_1 +19.5827350227176 Alpha3_3 +33.0369963798909 Alpha2_0
      +19.5827350227176 Alpha2_1 +15.4172649772824 Alpha2_0 +19.5827350227176 Alpha2_1 +15.4172649772824 Alpha2_0
      +15.4172649772824 Alpha2_0
Subject to
 sh2_0_0: +4 t_Z0 +4 t_wt2_0 +100 Alpha2_0 +50 Alpha2_0 >= +3
 sh2_0_1: +4 t_Z1 +4 t_wt2_0 +100 Alpha2_0 +50 Alpha2_0 >= +3
 sh2_1_0: +5 t_Z2 +5 t_wt2_1 +50 Alpha2_1 >= +4
 sh2_1_1: +5 t_Z3 +5 t_wt2_1 +50 Alpha2_1 >= +4
 sh3_0_0: +4 t_Z0 +150 Alpha2_0 +150 Alpha3_0 +99 Alpha2_0 +50 Alpha2_0 >= +3
 sh3_1_0: +5 t_Z1 +5 t_wt3_1 +150 Alpha2_0 +50 Alpha3_1 +100 Alpha2_0 +50 Alpha2_0 +50 Alpha2_0 >= +4
 sh3_2_0: +4 t_Z2 +4 t_wt3_2 +150 Alpha2_1 +150 Alpha3_2 +100 Alpha2_1 +50 Alpha2_1 >= +3
 sh3_3_0: +5 t_Z3 +5 t_wt3_3 +150 Alpha2_1 +50 Alpha3_3 +100 Alpha2_1 +50 Alpha3_3 +50 Alpha2_1 >= +4
 z0_2: +1 t_Z0 -1 Alpha2_0 -1 Alpha2_0 -1 Alpha2_0 >= +0
 z0_3: +1 t_Z0 -1 t_wt3_0 -1 Alpha3_0 >= +0
 z1_2: +1 t_Z1 -1 Alpha2_0 -1 Alpha2_0 -1 Alpha2_0 >= +0
 z1_3: +1 t_Z1 -1 Alpha3_1 >= +0
 z2_2: +1 t_Z2 -1 Alpha2_1 -1 Alpha2_1 -1 Alpha2_1 -1 Alpha2_1 >= +0
 z2_3: +1 t_Z2 -1 Alpha3_2 >= +0
 z3_2: +1 t_Z3 -1 Alpha2_1 -1 Alpha2_1 -1 Alpha2_1 -1 Alpha2_1 >= +0
 z3_3: +1 t_Z3 -1 Alpha3_3 >= +0
 zn: +0.194034334187542 t_Z0 +0.24645895087767 t_Z1 +0.24645895087767 t_Z2 +0.313047764057118 t_Z3 <= +0.25
 w2_0: +1 t_wt2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0 +1 Alpha2_0
  = +1
 w2_1: +1 t_wt2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 +1 Alpha2_1 = +1
 w3_0: +1 t_wt3_0 +1 Alpha3_0 +1 Alpha3_0 = +1
 w3_1: +1 t_wt3_1 +1 Alpha3_1 +1 Alpha3_1 = +1
 w3_2: +1 t_wt3_2 +1 Alpha3_2 +1 Alpha3_2 = +1
 w3_3: +1 t_wt3_3 +1 Alpha3_3 +1 Alpha3_3 +1 Alpha3_3 = +1
Bounds
 0 <= t_Z0 <= 1
 0 <= t_Z1 <= 1
 0 <= t_Z2 <= 1
 0 <= t_Z3 <= 1
End
