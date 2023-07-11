=======
Classes and functions
=======

.. autoclass:: beansp.Beans

  .. autofunction:: beansp.Beans.__init__

  .. autofunction:: beansp.Beans.do_run

  .. autofunction:: beansp.Beans.lnprior

  .. autofunction:: beansp.Beans.lnprob

  .. autofunction:: beansp.Beans.plot_model

    Produces a plot like so:

    .. image:: ./plot_model_example.png
      :width: 600

    The persistent flux measurements (*red dots*, left-hand *y*-axis) are
    shown, joined by lines implying the use of linear interpolation for
    flux inbetween.
    Fluence of the observed bursts are indicated (*gray circles*,
    right-hand *y*-axis) along with the predicted bursts (*blue stars*).
    The time of the referenc eburst is indicated (*black vertical line*).
    For the purposes of simulation the code assumes the accretion rate is
    constant between the predicted bursts, which is indicated by the
    stepped line. 

  .. autofunction:: beansp.Beans.do_analysis
