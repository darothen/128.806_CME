## Computational Methods Exercise 

**Author**: Daniel Rothenberg <darothen_at_mit_dot_edu>
**Latest Version**: Spring, 2015


This is a problem set for MIT course 12.806/12.306/10.571, "[Atmospheric Physics and Chemistry](http://ocw.mit.edu/courses/chemical-engineering/10-571j-atmospheric-physics-and-chemistry-spring-2006/)", designed to give students a chance to practice solving the types of problems they'll encounter later on in the course. Enclosed is a solution implemented in an IPython notebook.

### To build:

The writeup can be compiled to LaTeX by running `make`. The python code is self-contained; for instance, the following code block instantiates and runs the model using a pre-packaged, time-varying emissions function:

```python
    from integrator import RK4
    from carbon_model import CarbonModel, emissions_22p6

    y0 = np.array([612, 730, 140, 37000, 580, 1500, 5300, 1.0])
    emis_func = lambda t: emissions_22p6(t, 100.)
    model = CarbonModel(*y0, emis_func=emis_func)

    o = model.integrate(RK4, 150.)
    o.atm_ppm.plot()
```

### To build solution:

1. Run the IPython notebook `solution_with_notes.ipynb` to generate static content

2. Generate the basic markdown file from the ipython notebook by running 

```shell
ipython nbconvert --to markdown solution_with_notes.ipynb
```

3. Edit the markdown:

    - Replace all the indented code blocks with backtick style Python/shell blocks
    - Remove superfluous configuration code
    - Tweak environments so that all the LaTeX is nice

4. Convert the markdown to tex using pandoc via the provided Makefile

```shell
make pdf
```

5. Compile using xelatex

```shell
xelatex -shell-escape solution_with_notes.tex
```
