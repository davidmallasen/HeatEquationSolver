# HeatEquationSolver
Heat Equation Solver using MPI with Idle Period Propagation monitoring

## Building and testing
To compile and execute `heat_equation_solver.c` on Beskow run:
~~~
cc heat_equation_solver.c -o heat_equation_solver.out
srun -n 256 heat_equation_solver.out
~~~

There are some documented parameters at the beginning of `heat_equation_solver.c` to configure the run.

To check a simple 4x4 test compare `heat_equation_result_4.txt` with the output of:
~~~
python3 heat_equation_tester.py
~~~
This serial python script is a modification of the one found in [this url](https://scipython.com/book/chapter-7-matplotlib/examples/the-two-dimensional-diffusion-equation/) and can be parametrized similarly to `heat_equation_solver.c`.

To check a larger execution visually, save the debug output:
~~~
srun -n 256 heat_equation_solver.out > heat_equation_result.txt
~~~
And visualize using the python script:
~~~
python3 visualize_heat_equation.py
~~~
