## TODO

# By priority order in short term
- Correctly take into account boundary conditions (ghost cells added before interfacing).
- Implement Constrained Transport.
- Implement time integrator.
- Use real physical initial conditions aswell as realistic values for P and T.
- Implement diagnostics and unit testing.
- Post-processing
- Computation of new dt on each loop for CFL condition ?

# And then
- Extend the model (3D, Energy Equation instead of hard coded P, Non-Ideal MHD, more sofisticated Riemann solvers).
- Clean up the code.
- Integrate in PHARE.