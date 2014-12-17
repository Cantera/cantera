.. default-role:: math

****************
Reactor Networks
****************

Cantera's Reactor Network module is designed to simulate networks of
interconnected reactors. The contents of each reactor in the network are
assumed to be homogeneous, a model variously referred to as the Continuously
Stirred Tank Reactor (CSTR), Well-Stirred Reactor (WSR), or Perfectly Stirred
Reactor (PSR) model. Cantera solves the time-dependent governing equations
that describe the evolution of the chemical and thermodynamic state of the
reactors.

The contents of each reactor can undergo chemical reactions according to a
specified kinetic mechanism, and surface reactions may occur on the reactor
walls. Each reactor in a network may be connected so that the contents of one
reactor flow into another. Reactors may be also be in contact with one another
or the environment via walls which move or conduct heat.

The purpose of this document is to describe the governing equations of reactor
models as implemented in Cantera.

Wall Interactions
=================

At each wall where there are surface reactions, there is a net generation (or
destruction) of homogeneous phase species. The molar rate of production for
each species `k` on wall `w` is `\dot{s}_{k,w}`. The total (mass) production
rate for species `k` on all walls is:

.. math::

    \dot{m}_{k,wall} = W_k \sum_w A_w \dot{s}_{k,w}

where `W_k` is the molecular weight of species `k` and `A_w` is the area of
each wall. The net mass flux from all walls is then:

.. math::

    \dot{m}_{wall} = \sum_k \dot{m}_{k,wall}

The total rate of heat transfer through all walls is:

.. math::

    \dot{Q} = \sum_w f_w \dot{Q}_w

General Reactor
===============

The state variables for Cantera's general reactor model are

    - `m`, the mass of the reactor's contents
    - `V`, the reactor volume
    - `U`, the total internal energy of the reactors contents
    - `Y_k`, the mass fractions for each species

Reactor Volume
--------------

The reactor volume changes as a function of time due to the motion of one or
more walls:

.. math::

    \frac{dV}{dt} = \sum_w f_w A_w v_w(t)

where `f_w = \pm 1` indicates the facing of the wall, `A_w` is the surface
area of the wall, and `v_w(t)` is the velocity of the wall as a function of
time.

Mass Conservation
-----------------

The total mass of the reactor's contents changes as a result of flow through
the reactor's inlets and outlets, and production of homogeneous phase species
on the reactor walls:

.. math::

    \frac{dm}{dt} = \sum_{in} \dot{m}_{in} - \sum_{out} \dot{m}_{out} +
                    \dot{m}_{wall}

Species Conservation
--------------------

The rate at which species `k` is generated through homogeneous phase reactions
is `V \dot{\omega}_k W_k`, and the total rate at which species `k` is
generated is:

.. math::

    \dot{m}_{k,gen} = V \dot{\omega}_k W_k + \dot{m}_{k,wall}

The rate of change in the mass of each species is:

.. math::

    \frac{d(mY_k)}{dt} = \sum_{in} \dot{m}_{in} Y_{k,in} -
                         \sum_{out} \dot{m}_{out} Y_k +
                         \dot{m}_{k,gen}

Expanding the derivative on the left hand side and substituting the equation
for `dm/dt`, the equation for each homogeneous phase species is:

.. math::

    m \frac{dY}{dt} = \sum_{in} \dot{m}_{in} (Y_{k,in} - Y_k)+
                      \dot{m}_{k,gen} - Y_k \dot{m}_{wall}

Energy Conservation
-------------------

The equation for the total internal energy is found by writing the first law
for an open system:

.. math::

    \frac{dU}{dt} = - p \frac{dV}{dt} - \dot{Q} +
                    \sum_{in} \dot{m}_{in} h_{in} - h \sum_{out} \dot{m}_{out}

Ideal Gas Reactor
=================

The Ideal Gas Reactor model is similar to the General Reactor model, with the
reactor temperature `T` replacing the total internal energy `U` as a state
variable. For an ideal gas, we can rewrite the total internal energy in terms
of the mass fractions and temperature:

.. math::

    U = m \sum_k Y_k u_k(T)

    \frac{dU}{dt} = u \frac{dm}{dt}
                    + m c_v \frac{dT}{dt}
                    + m \sum_k u_k \frac{dY_k}{dt}

Substituting the corresponding derivatives yields an equation for the
temperature:

.. math::

    m c_v \frac{dT}{dt} = - p \frac{dV}{dt} - \dot{Q}
        + \sum_{in} \dot{m}_{in} \left( h_{in} - \sum_k u_k Y_{k,in} \right)
        - \frac{p V}{m} \sum_{out} \dot{m}_{out} - \sum_k \dot{m}_{k,gen} u_k

While this form of the energy equation is somewhat more complicated, it
significantly reduces the cost of evaluating the system Jacobian, since the
derivatives of the species equations are taken at constant temperature instead
of constant internal energy.

Constant Pressure Reactor
=========================

For this reactor model, the pressure is held constant. The volume is not a
state variable, but instead takes on whatever value is consistent with holding
the pressure constant. The total enthalpy replaces the total internal energy
as a state variable. Using the definition of the total enthalpy:

.. math::

    H = U + pV

    \frac{dH}{dt} = p \frac{dV}{dt} + V \frac{dp}{dt}

Noting that `dp/dt = 0` and substituting into the energy equation yields:

.. math::

   \frac{dH}{dt} = - \dot{Q} + \sum_{in} \dot{m}_{in} h_{in}
                   - h \sum_{out} \dot{m}_{out}

The species and continuity equations are the same as for the general reactor
model.

Ideal Gas Constant Pressure Reactor
===================================

As for the Ideal Gas Reactor, we replace the total enthalpy as a state
variable with the temperature by writing the total enthalpy in terms of the
mass fractions and temperature:

.. math::

    H = m \sum_k Y_k h_k(T)

    \frac{dH}{dt} = h \frac{dm}{dt} + m c_p \frac{dT}{dt}
                    + m \sum_k h_k \frac{dY_k}{dt}

Substituting the corresponding derivatives yields an equation for the
temperature:

.. math::

    m c_p \frac{dT}{dt} = - \dot{Q} - \sum_k h_k \dot{m}_{k,gen}
        + \sum_{in} \dot{m}_{in} \left(h_{in} - \sum_k h_k Y_{k,in} \right)
