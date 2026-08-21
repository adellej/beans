.. highlight:: shell

===============================
Latest results and bibliography
===============================

Since the beansp code has gone through several generations, we collect here
the latest results for each application that has been published. We
consider these as the best estimates of the system parameters according to
the data comparison with the chosen model (meaning, if the model is
incorrect, then the system parameters will similarly be untrustworthy).

1. SAX J1808.4-3658
-------------------

The 2002 outburst of this source was the original target for the analysis
which preceded development of the beansp code. Four bursts were detected
with *RXTE*/PCA, and the corresponding persistent flux history (along with
a bolometric correction to the 3-25 keV flux of 2.21) is used as
the input to the code.

Shown below are the parameter limits for the latest run, ``base_newalpha``
from 2026 August. This analysis is based on the new alpha values stored in
``data/1808_bursts_newalpha.txt``; ``base_newalpha.ini`` and posterior
samples can be found in the data accompanying :ref:`Galloway et al.  (2026) <gal26>`.

The prior on :math:`Z_{\rm CNO}` label “beta“
indicates the same Beta-function as adopted by earlier analyses with beansp.
The :math:`M_{\rm NS}`, :math:`R_{\rm NS}` prior label “steiner” indicates the same
2-dimensional prior derived from the analysis of 
`Steiner et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.476..421S>`_
for neutron stars in globular clusters, as also used by the same analyses.

..
    Table assembled manually to explore format

.. table:: System parameters derived from the 2002 outburst of SAX J1808.4–3658
   :widths: auto

   =================== ================================= ======================= ==================================
   Parameter           Units                             Prior                   base_newalpha
   =================== ================================= ======================= ==================================
   :math:`X`           \                                 :math:`U[10^{-5},0.76]` :math:`0.41_{-0.10}^{+0.13}` 
   :math:`Z_{\rm CNO}` \                                 beta                    :math:`0.0151_{-0.0067}^{+0.0119}` 
   :math:`Q_{\rm b}`   MeV/nucleon                       :math:`U[10^{-6},5.0]`  :math:`0.6_{-0.2}^{+0.3}` 
   :math:`d`           kpc                               :math:`U[1,20]`         :math:`3.2\pm0.3` 
   :math:`\xi_b`       \                                 :math:`U[0.01,2]`       :math:`1.03_{-0.17}^{+0.20}` 
   :math:`\xi_p`       \                                 :math:`U[0.01,10]`      :math:`1.6\pm0.4` 
   :math:`M_{\rm NS}`  :math:`M_\odot`                   steiner                 :math:`1.8_{-0.4}^{+0.5}` 
   :math:`R_{\rm NS}`  km                                steiner                 :math:`12.2\pm1.2` 
   :math:`g`           :math:`10^{14}\ {\rm cm\,s^{-2}}` \                       :math:`2.2_{-0.4}^{+0.7}` 
   :math:`1+z`         \                                 \                       :math:`1.338_{-0.082}^{+0.111}` 
   =================== ================================= ======================= ==================================


2. IGR J17498-2921
------------------

Eight bursts were observed during the 2011 outburst of this, the second
application for beansp.

Shown below are the parameter limits for the latest run, ``base_alpha``
from 2026 August. This analysis is based on the new alpha values stored in
``data/17498_bursts_newalpha.txt``; the config file ``base_alpha.ini`` and posterior
samples can be found in the data accompanying :ref:`Galloway et al.  (2026) <gal26>`.

.. table:: System parameters derived from the 2011 outburst of IGR J17498-2921
    :widths: auto

    =================== ================================= ======================= ===============================
    Parameter           Units                             Prior                   base_alpha                     
    =================== ================================= ======================= ===============================
    :math:`X`           \                                 :math:`U[10^{-5},0.76]` :math:`0.17_{-0.06}^{+0.16}`   
    :math:`Z_{\rm CNO}` \                                 beta                    :math:`0.004_{-0.002}^{+0.019}`
    :math:`Q_b`         MeV/nucleon                       :math:`U[10^{-6},5]`    :math:`2.6_{-0.5}^{+0.9}`      
    :math:`d`           kpc                               :math:`U[1,20]`         :math:`4.9\pm0.4`              
    :math:`\xi_b`       \                                 :math:`U[0.01,2]`       :math:`1.83_{-0.19}^{+0.12}`   
    :math:`\xi_p`       \                                 :math:`U[0.01,10]`      :math:`1.6_{-0.3}^{+0.2}`      
    :math:`M_{\rm NS}`  :math:`M_\odot`                   steiner                 :math:`2.1_{-0.5}^{+0.3}`      
    :math:`R_{\rm NS}`  km                                steiner                 :math:`12.1_{-0.8}^{+1.1}`     
    :math:`\cos i`      \                                 \                       :math:`0.61_{-0.16}^{+0.42}`   
    :math:`g`           :math:`10^{14}\ {\rm cm\,s^{-2}}` \                       :math:`2.5_{-0.6}^{+0.8}`      
    :math:`1+z`         \                                 \                       :math:`1.42\pm0.11`            
    =================== ================================= ======================= ===============================


3. SRGA J144459.2-604207
------------------------

We assembled a set of 14 daily burst epochs covering the 2024 outburst of 
SRGA J144459.2−604207, and attempted to constrain the system properties of
this object by analysing a subset of the epochs with beansp, in "ensemble"
mode. 

Shown below are the parameter limits for the run ``base_sel_mr``
based on the epoch averages stored in
``data/144459_ensbl_selected.txt``; the config file ``base_sel_mr.ini``
and posterior samples can be found in the data accompanying :ref:`Galloway
et al.  (2026) <gal26>`.

.. table:: System parameters derived from the 2024 outburst of SRGA J144459.2−604207
    :widths: auto

    =================== ================================= ======================== ============================
    Parameter           Units                             Prior                    base_sel_mr_prod            
    =================== ================================= ======================== ============================
    :math:`X`           \                                 :math:`U[10^{-5},0.76]`  :math:`0.54_{-0.06}^{+0.08}`
    :math:`Z_{\rm CNO}` \                                 :math:`U[10^{-5},0.056]` :math:`0.0109\pm0.0019`     
    :math:`Q_b`         MeV/nucleon                       :math:`U[10^{-6},5]`     :math:`0.6_{-0.2}^{+0.4}`   
    :math:`d`           kpc                               :math:`U[1,20]`          :math:`10.9_{-1.7}^{+2.1}`  
    :math:`\xi_b`       \                                 :math:`U[0.01,2]`        :math:`1.2\pm0.4`           
    :math:`\xi_p`       \                                 :math:`U[0.01,10]`       :math:`1.2\pm0.4`           
    :math:`M_{\rm NS}`  :math:`M_\odot`                   :math:`U[1.15,2.5]`      :math:`1.9\pm0.4`           
    :math:`R_{\rm NS}`  km                                :math:`U[9,17]`          :math:`11.8_{-1.9}^{+2.2}`  
    :math:`\cos i`      \                                 \                        :math:`0.4_{-0.2}^{+0.4}`   
    :math:`g`           :math:`10^{14}\ {\rm cm\,s^{-2}}` \                        :math:`2.3_{-0.6}^{+1.4}`   
    :math:`1+z`         \                                 \                        :math:`1.36_{-0.08}^{+0.18}`
    =================== ================================= ======================== ============================

4. Simulated data
-----------------

A key obstacle leaving questions
over the reliability of the beansp burst matching algorithm is the lack of well
characterised calibration data to test on. We seek to establish whether, in the
case where we know the system parameters giving rise to a set of burst data,
that we can recover them accurately via beansp. A secondary objective is
to understand the systematic uncertainties that may arise between different
numerical models.

In the most recent work we have approached these questions with two
datasets; first, KEPLER simulations designed to replicate the bursts
during the 2002 outburst of SAX J1808.4-3658, by Johnston et al.
(2018).
Parameter constraints for the beansp run designed to replicate those
analysis are shown below.

.. table:: System parameters for KEPLER simulations of SAX J1808.4-3658
    :widths: auto

    =================== ================================= ===== ======================= ===============================
    Parameter           Units                             Input Prior                   johnston18                     
    =================== ================================= ===== ======================= ===============================
    :math:`X`           \                                 0.44  :math:`U[10^{-5},0.76]` :math:`0.56_{-0.10}^{+0.11}`   
    :math:`Z_{\rm CNO}` \                                 0.02  beta                    :math:`0.015_{-0.005}^{+0.008}`
    :math:`Q_b`         MeV/nucleon                       0.3   :math:`U[10^{-6},5]`    :math:`0.5\pm0.3`              
    :math:`d`           kpc                               3.5   :math:`U[1,20]`         :math:`3.3_{-0.3}^{+0.4}`      
    :math:`\xi_b`       \                                 (1.0) :math:`U[0.01,2]`       :math:`0.82_{-0.19}^{+0.34}`   
    :math:`\xi_p`       \                                 1.1   :math:`U[0.01,10]`      :math:`1.4\pm0.3`              
    :math:`M_{\rm NS}`  :math:`M_\odot`                   1.4   steiner                 :math:`2.2_{-0.4}^{+0.3}`      
    :math:`R_{\rm NS}`  km                                11.2  steiner                 :math:`11.3_{-1.0}^{+1.1}`     
    :math:`\cos i`      \                                 \     \                       :math:`0.21_{-0.06}^{+0.14}`   
    :math:`g`           :math:`10^{14}\ {\rm cm\,s^{-2}}` \     \                       :math:`3.2_{-0.8}^{+0.9}`      
    :math:`1+z`         \                                 \     \                       :math:`1.49_{-0.13}^{+0.14}`   
    =================== ================================= ===== ======================= ===============================

We also performed "ensemble"-mode runs based on SETTLE-simulated data intended to
mimic the properties of bursts observed from SRGA J144459.2-604207 during
it's 2024 outburst. Below are the results from one of the simulation runs.

.. table:: System parameters for SETTLE simulations of SRGA J144459.2-604207
    :widths: auto

    =================== ================================= ===== ======================== ===============================
    Parameter           Units                             Input Prior                    sim1_alpha_mr             
    =================== ================================= ===== ======================== ===============================
    :math:`X`           \                                 0.20  :math:`U[10^{-5},0.76]`  :math:`0.22_{-0.05}^{+0.04}`   
    :math:`Z_{\rm CNO}` \                                 0.016 :math:`U[10^{-5},0.056]` :math:`0.019_{-0.005}^{+0.004}`
    :math:`Q_b`         MeV/nucleon                       0.15  :math:`U[10^{-6},5]`     :math:`0.16_{-0.09}^{+0.15}`   
    :math:`d`           kpc                               7.43  :math:`U[1,20]`          :math:`11_{-2}^{+3}`           
    :math:`\xi_b`       \                                 1.9   :math:`U[0.01,2]`        :math:`0.9_{-0.4}^{+0.6}`      
    :math:`\xi_p`       \                                 2.96  :math:`U[0.01,10]`       :math:`1.7_{-0.6}^{+0.9}`      
    :math:`M_{\rm NS}`  :math:`M_\odot`                   1.4   :math:`U[1.15,2.5]`      :math:`1.8\pm0.4`              
    :math:`R_{\rm NS}`  km                                11.2  :math:`U[9,17]`          :math:`12.2_{-1.9}^{+2.7}`     
    :math:`\cos i`      \                                 \     \                        :math:`0.18_{-0.05}^{+0.08}`   
    :math:`g`           :math:`10^{14}\ {\rm cm\,s^{-2}}` \     \                        :math:`2.1_{-0.6}^{+0.9}`      
    :math:`1+z`         \                                 \     \                        :math:`1.33_{-0.08}^{+0.11}`   
    =================== ================================= ===== ======================== ===============================


Bibliography
------------

.. _gal26:
* Galloway et al. (submitted to PASA, 2026) and accompanying data at `Monash bridges <https://doi.org/10.26180/33065816>`_
* `Galloway et al. (2024) <https://10.1093/mnras/stae2422>`_ and accompanying data at `Monash bridges <10.1093/mnras/stae2422>`_
* `Goodwin et al. (MNRAS 490, 2228, 2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.2228G>`_
* `Johnston et al. (MNRAS 477, 2112, 2018) <http://adsabs.harvard.edu/abs/2018MNRAS.477.2112J>`_
* `Galloway & Cumming (ApJ 652, 559, 2006) <http://adsabs.harvard.edu/abs/2006ApJ...652..559>`_

