
Compartmental Model of the Global Carbon Cycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model was implemented separately (see the appendices). Your mileage
may vary in the following portion of the exercise, but in general I'd
recommend wrapping your model so that it only takes 1-2 lines of code to
initialize the model and then run it, returning all the output you need
to analyze the results.

Problem 1
^^^^^^^^^

In the framework of a simple box-model like the one we've written here,
one can analyze the lifetime of species in each reservoir by collecting
the all the sources and sinks and re-writing the equation in the form

.. math:: \frac{dC}{dt} = \text{Sources} - \text{Sinks}

where "Sources" and "Sinks" are all positive and negative terms,
respectively, in the original system of ODEs. Neglecting the source term
and considering the fact that most of the sinks are going to be
represented in the form :math:`k_i X`, we can represent the lifetime as

.. math:: \tau = \frac{X}{\sum\limits_{i=1}^n k_i}

:math:`k_i` is, in our case, each of the arrows leading *away* from a
given reservoir in Figure 22.6 of Seinfeld and Pandis. This formula
yields the following lifetimes for each reservoir:

.. code-block:: python

    from IPython.display import Latex, display

    reservoirs = {
        # name -> [amount (PgC), removal rate (PgC/yr)
        'M1': [612., 100. + 57. + 19],
        'M2': [730., 12. + 57. + 58.],
        'M3': [140., 100. + 18.],
        'M4': [37000., 40. + 70.],
        'M5': [580., 50. + 50.],
        'M6': [1500., 50.],
    }


    format_str = r"$\tau$(M$_{n:1d}$) = {tau:>5.1f} years"
    for n in xrange(1, 7):
        [X, ks] = reservoirs['M%1d' % n]
        k = np.sum(ks)
        tau = X/ks
        display(Latex(format_str.format(n=int(n), tau=tau)))

.. math::
    \tau(M_1) &=   3.5 \mathrm{years} \\
    \tau(M_2) &=   5.7 \mathrm{years} \\
    \tau(M_3) &=   1.2 \mathrm{years} \\
    \tau(M_4) &= 336.4 \mathrm{years} \\
    \tau(M_5) &=   5.8 \mathrm{years} \\
    \tau(M_6) &=  30.0 \mathrm{years} \\


Problem 2
^^^^^^^^^

There are two different surface ocean reservoirs because the equilibrium
dissolution of |CO2| in seawater is dependent on the water's
temperature and pH.

|CO2| hydrolizes when it disolves in seawater:

.. raw:: latex

   \begin{align*}
    \cee{CO2(g) + H2O &<=> CO2.H2O} \\
    \cee{CO2.H2O & <=> H+ + HCO3-} \\
    \cee{HCO3- &<=> H+ + CO3^{2-}}
   \end{align*}

This produces both carbonate and bicarbonate ions. The abundance of
dissolved salts in the ocean affects pH locally between a range of about
7.5 and 8.4, which is further impacted by temperature. Also, this
reaction system is buffered; as :math:`\ce{CO2}` dissolves in seawater,
that reservoir can take up less and less :math:`\ce{CO2}` from the
atmosphere. The complexity of this reaction sequence is simplified in
the compartmental model using a simple parameterization developed by
`Ver et al
(1999) <http://earth.geology.yale.edu/~ajs/1999/07-09.1999.11Ver.pdf>`__,
which modifies the loss mechanism of atmospheric :math:`\ce{CO2}` via
dissolution in seawater as :math:`F = kM^\beta`, where :math:`\beta` is
a positive constant which accounts for this complex chemistry. Using the
two reservoirs and suitable associated :math:`\beta`'s, the
compartmental model can attempt to simulate oceanic carbon sink.

Problem 3
^^^^^^^^^

Assume that all the carbon in the atmosphere is present in the form of
|CO2|. The (volume) mixing ratio of a gas in the air is given
as the ratio of a trace constituents molar concentration to that of the
full gas. Note that this is a useful quantity to know because as the
density of the air changes, the mixing ratio will remain the same, so
one can always back out the mass of a given trace constituent with
minimal extra information.

To a good approximation, the mean molecular weight of air can be
computed just from the relative abundance and molecular weights of
N:sub:`2`, |CO2|, and Ar (which total 99.3%
of the atmosphere by mass or volume):

.. math::

    M_a = (0.78 \cdot 28 \mathrm{g/mol}) + (0.21 \cdot 32 \mathrm{g/mol}) + (0.01 \cdot 40 \mathrm{g/mol}) \approx 28.97 \mathrm{g/mol}

(note that we rounded up the abundance of argon).

The molar mass of carbon is 12.01 g/mol; it's just about 16 g/mol for
elemental oxygen, so carbon dioxide has a molar mass of 44.01 g/mol. If
the pre-industrial atmospheric carbon burden was 612 Pg, then we can
multiply by the molar ratio of carbon to carbon dioxide to yield the
mass of |CO2| in the atmosphere,

.. math::

    \mathrm{Mass(\ce{CO2})} &= 612 \mathrm{ Pg(C)} \times
        \frac{44.01 \mathrm{g(\ce{CO2})/mol}}{12.01 \mathrm{ g(C)/mol}} \\
                            &= 2244 \mathrm{ Pg(\ce{CO2})}

Since we know the mass of the atmosphere is roughly
5.15 x 10^18 kg, the mass mixing ratio of |CO2| is just

.. math::

    r(\ce{CO2}) = \frac{224.3\mathrm{ Pg}}{5.15\times10^{18}\mathrm{ kg}} = 435.73 \mathrm{ ppm(m)}

To convert equivalently between mass mixing ratio and molar/volume
mixing ratio, we have to account for difference in number concentrations
of the molecules of each constituent, which are going to differ because
they have molar mass. However, this computation can be radically
shortened by noting that the mass mixing ratio is approximately equal to
the molar mixing ratio times the ratio of the molar weight of the trace
gas to the average molar weight of the entire gas,

.. math::

    \chi(\ce{\mathrm{X}})\frac{M(\mathrm{X})}{M_a} = \mathrm{r(X)}

Substituting in our intermediate values, we get

.. math::

    \chi(\ce{CO2}) = 435.73 \mathrm{ ppm(m)} \times \frac{28.97\mathrm{ g/mol}}{44.01\mathrm{ g/mol}} \approx  286.8 \mathrm{ ppm(v)}

which is a very reasonable estimate for pre-industrial |CO2|
based on numerous proxy records.
