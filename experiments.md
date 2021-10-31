# Experiments

This page shall document relevant experiments using the symbol timing recovery loop simulator. The goal is to highlight some considerations in terms of parameters and the expected performance.

## Performance under Sampling Clock Frequency Offset

The interpolator choice can be critical when applying the symbol timing recovery scheme under significant sampling clock frequency offset with low oversampling ratios such as $L=2$. In this scenario, due to the clock frequency offset, the fractional timing offset mu keeps increasing or decreasing, depending on the sign of the offset. Hence, $\mu$ changes rapidly and frequently wraps around from 0 to 1 or vice-versa. In the end, the interpolator often has to compute interpolants far from the basepoint index.

For instance, consider a waveform sampled at a rate of 2 samples/symbol and that the fractional timing offset estimate reaches $\mu(k) = 1.0$. In this case, the interpolator would need to estimate the interpolant located after an interval of $0.5 T_s$ from the basepoint index, where $T_s$ is the symbol period, i.e., half a symbol period (or one sample period). This could be a challenging task for a linear interpolator, for example. In contrast, this problem does not arise when the oversampling ratio is high. For instance, with 8 samples/symbol, the interval resulting from $\mu(k) = 1.0$ is of only $0.125T_s$, which seems doable for a linear interpolator.

To evaluate this, configure `main.m` as follows:
- `L = 2`
- `timeOffset = 1`
- `fsOffsetPpm = 100`
- `EsN0 = 50`
- `intpl = 1`

> Note: The high $E_s/N_0$ is useful to focus on the degradation due to the
> symbol timing recovery scheme alone (instead of noise).

With this configuration, you can observe that:
- The performance achieved with this project's timing recovery loop is poor compared to the performance achieved by MATLAB's `comm.SymbolSynchronizer`. That is only because the latter does not offer a configurable interpolator. So, while this project's recovery loop uses a linear interpolator (`intpl = 1`), the MATLAB implementation uses its default quadratic interpolator. This performance discrepancy is already a hint about the importance of using a higher-order interpolator in this scenario.
- If you enable the real-time debugging option `debug_tl_runtime = 1`, you can see the constellation periodically loses accuracy. The rationale is that the linear interpolator can only perform well when the target interpolant falls on a line between the basepoint and the succeeding samples. Thus, the interpolation becomes poor periodically as the fractional timing offset ramps up from 0 to 1. In other words, the raised cosine pulse is a "rounded" waveform and not a triangular wave, so the interpolation over long intervals is not always accurate.

Next, observe what happens when a quadratic interpolator is adopted instead. Change the following interpolator parameter on `main.m`:
- `intpl = 2`

Now, the performance is significantly superior and similar to the performance achieved by the MATLAB implementation. After all, now both implementations are using quadratic interpolators. Nevertheless, note there is a remaining artifact of the imperfect interpolation. If you observe the real-time constellation scope (with option `debug_tl_runtime = 1`), you can see that the constellation points are periodically contracting (getting closer together) and expanding, as if they were sliding on diagonal axes around each of the four QPSK constellation points. This phenomenon is also visible in the static constellation plot obtained with `debug_tl_static = 1` on `main.m`. In contrast, with a cubic interpolator, the referred problem disappears. To verify that, set the following parameter on `main.m`:
- `intpl = 3`

In fact, you can observe the measured MER (printed on the console) achieved with the cubic interpolator is superior to the MER achieved with MATLAB's implementation using a quadratic interpolator.

The performance is even better when using a polyphase interpolator. To verify that, set the following parameter on `main.m`:
- `intpl = 0`

The rationale is that the polyphase interpolator has an interpolation factor of its own. The receiver may be running with an oversampling of two (`L=2`). Meanwhile, the polyphase interpolator can still rely on 32 polyphase branches, which is equivalent to using an interpolating factor of 32. Because of this property, the polyphase interpolator achieves the best performance with the parameters in this experiment. Its measured MER is nearly 10 dB better than the one achieved with the cubic interpolator, and 20 dB better than the linear interpolator.

Lastly, observe that these effects arise only because the oversampling ratio is low. For example, set the linear interpolator again but increase the oversampling ratio to $L=8$. That is, on `main.m`, set:
- `L = 8`
- `intpl = 1`

Note the constellation is clean now, even though the interpolator is a linear one.

