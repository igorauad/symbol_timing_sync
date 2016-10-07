# Symbol Timing Synchronization Simulations

Several scripts used to study the material from the Book *" Digital Communications: A Discrete-Time Approach"*, by *Michael Rice*

Main simulation implements Symbol Timing Recovery using either a
**Maximum-likelihood (ML) Timing Error Detector (ML-TED)** or a **Zero-Crossing
TED (ZC-TED)**. The loop filter is a **Proportional-plus-integrator (PI)
Controller** and the interpolator can be chosen as a Linear Interpolator or a
**Polyphase Interpolator**. The **Interpolator Controller** is a **Modulo-1
Counter**.

| File        | Description         |
| ------------- |:--------------|
| `symbol_synchronizer.m`     | Main Simulation |
| `symTimingLoop.m`     | Function that implements the timing recovery loop. |
| `getTedKp.m`     | Function that computes the Timing Error Detector (TED) gain . |
| `timingLoopPIConstants.m`     | Function that computes the PI controller constants. |
| `polyphaseFilterBank.m`     | Function that computes the polyphase subfilters that are required when the polyphase interpolator is adopted. |
| `tedDesign.m`     | A short script to analyze TED design parameters. |