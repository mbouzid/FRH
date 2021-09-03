# FRH
Fix-and-relax heuristic applied to an order acceptance scheduling problem under TOU tarrifs and taxed carbon emissions periods. 
The heuristic strategy involves the sequence-dependent setup times.

# Usage

```./FRH [instance.dat] [observation_window_size] [overlapping_steps] [MAX|MIN] ```

e.g.

```./FRH Dataslack_15orders_Tao1R1_1.dat 20 3 MAX```

# Prerequisites

-- IBM Cplex v.12.9

## Input

The input file format is .dat (from Cplex Studio).   

## Output

Console output. Returns the completion times and the sequence in the following format.

```
C[i]==x;
C[j]==y;
u[i][j]==1;
```

The objective value and the solving time is also given.



