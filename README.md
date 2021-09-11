# Symbol Timing Recovery

This repository contains MATLAB scripts focusing on symbol timing recovery
algorithms. The current implementation is based on the material from the book
*"Digital Communications: A Discrete-Time Approach"*, by Michael Rice.

The main script is a symbol timing recovery simulator, which can evaluate two
timing error detectors (TEDs): a maximum-likelihood (ML) TED (**ML-TED**) and a
zero-crossing TED (**ZC-TED**)*. The loop filter adopted by the simulator is a
proportional-plus-integrator (PI) controller, and the interpolator can be chosen
from linear, polyphase, quadratic, and cubic interpolator options. The
interpolator controller is a modulo-1 counter.

| File                        | Description                                           |
| --------------------------- |:------------------------------------------------------|
| `symbol_synchronizer.m`     | Main simulation.                                      |
| `symTimingLoop.m`           | Function that implements the timing recovery loop.    |
| `calcTedKp.m`               | Function to compute the timing error detector gain.   |
| `piLoopConstants.m`         | Function to compute the PI controller constants.      |
| `polyphaseFilterBank.m`     | Function to design the polyphase interpolator.        |
| `plotTedGain.m`             | Function to plot the TED gain vs. the rolloff factor. |
| `calcSCurve.m`              | Function to compute the TED's S-curve analytically.   |
| `simSCurve.m`               | Function to simulate the TED's S-curve.               |
| `sCurveDemo.mlx`            | A demonstration of TED S-curve and gain evaluations.  |
