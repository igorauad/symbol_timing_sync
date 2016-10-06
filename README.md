# Symbol Timing Synchronization Simulations

Several scripts used to study the material from the Book *" Digital Communications: A Discrete-Time Approach"*, by *Michael Rice*

Main simulation implements Symbol Timing Recovery using a **Maximum-likelihood
(ML) Timing Error Detector (ML-TED)**, a **Proportional-plus-integrator (PI)
Controller**, an interpolator that can be chosen as the Linear Interpolator or
the **Polyphase Interpolator** and a **Modulo-1 Counter** as **Interpolator
Controller**.

| File        | Description         |
| ------------- |:--------------|
| `symbol_synchronizer.m`     | Main Simulation |
| `symTimingLoop.m`     | Function that implements the timing recovery loop. |
| `getTedKp.m`     | Function that computes the Timing Error Detector (TED) gain . |
| `timingLoopPIConstants.m`     | Function that computes the PI controller constants. |
| `polyphaseFilterBank.m`     | Function that computes the polyphase subfilters that are required when the polyphase interpolator is adopted. |
| `tedDesign.m`     | A short script to analyze TED design parameters. |