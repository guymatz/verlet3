# verlet3

## Assignment

Write a Fortran program that implements the Velocity Verlet algorithm
for calculating classical trajectories for the reaction A + BC
$\rightarrow$ AB + C, with A, B, and C being hydrogen atoms.
Use the subroutine `jpca15` in file `jpca15.f` (Rampino S, The
Journal of Physical Chemistry A 120,
[4683-4692](http://dx.doi.org/10.1021/acs.jpca.5b10018), 2016) to
calculate the interaction potential between the three atoms. Try
different initial conditions and collect at least one non-reactive
and one reactive trajectory.

### Guidelines and tips

#### The `jpca15` subroutine

The forces can easily be calculated as the derivative of the
potential (see `verlet-1.md`). The subroutine `jpca15.f` takes as
input argument a vector with the three interatomic distances (AB, AC,
and BC) and returns as output the potential energy and the vector of
the derivatives of the potential with respect to the interatomic
distances (AB, AC, and BC). Distances are in bohr and energies are in
eV.

#### Visualizing the trajectories.

Save the trajectories in XYZ format (see `verlet-2.md`) and
visualize them through a molecular visualize (e.g.,
[VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)).

## Notes

### Audio recordings
1. https://recorder.google.com/980c897d-2a57-4002-a5ad-3eff0bfe72cd
2. https://recorder.google.com/234570c9-b3eb-49eb-87f5-a2d1ac9dcce6

### Initial Conditions
I am starting with the following initial conditions (in data/reactive.dat):

```
6000 1.0                                     ! number of iterations, tau
3                                            ! number of atoms
                                             ! mass, location xyz, velocity xyz
1.0080 -10.0  0.0  0.0  1.0  0.0  0.0        ! m,    x, y, z,      vx, vy, vz
1.0080  10.0  0.0  0.0  0.0  0.0  0.0        ! m,    x, y, z,      vx, vy, vz
1.0080 11.40065  0.0  0.0  0.0  0.0  0.0        ! m,    x, y, z,      vx, vy, vz
```

### RUN
For command-line usage:
fpm run -- -h

fpm run -- -v -f app/atoms.dat -x

### TEST
fpm test
