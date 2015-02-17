## Computational Methods Exercise 

**Author**: Daniel Rothenberg <darothen_at_mit_dot_edu>
**Latest Version**: Spring, 2015


This is a problem set for MIT course 12.806/12.306/10.571, "[Atmospheric Physics and Chemistry](http://ocw.mit.edu/courses/chemical-engineering/10-571j-atmospheric-physics-and-chemistry-spring-2006/)", designed to give students a chance to practice solving the types of problems they'll encounter later on in the course. Enclosed is a solution implemented in an IPython notebook.

### To build:

The writeup can be compiled to LaTeX by running `make`. The python code is self-contained; for instance, the following code block instantiates and runs the model using a pre-packaged, time-varying emissions function:

```python
    from carbon_model import CarbonModel, emissions_22p6

    y0 = np.array([612, 730, 140, 37000, 580, 1500, 5300, 1.0])
    emis_func = lambda t: emissions_22p6(t, 100.)
    model = CarbonModel(*y0, emis_func=emis_func)

    o = model.integrate(odeint, 140.)
    o.atm_ppm.plot()
```