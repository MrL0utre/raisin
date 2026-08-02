# RaiSin test fixtures

`fixtures/chain.rzn2` is the hand-checkable scientific reference used by
`test_scientific.c`. It is a three-vertex path with delays 10, 20 and 5:

```text
red(10) -> black(20) -> red(5)
```

With all vertices on one FPGA, its longest-path cost is `10 + 20 + 5 = 35`
and its connectivity-minus-one cost is zero. With the parts `[0, 2, 0]` on
`arch_path.arch`, each of the two hyperedges traverses the two-link route
`0 -> 1 -> 2`. The route delay is 4 in each direction because link `0-1` has
zero delay and link `1-2` has delay 4. The expected longest-path cost is thus
`35 + 4 + 4 = 43`, the connectivity cost is 2, and both physical links carry
two signals.

`test_capacity.c` uses similarly small arrays to check resource and communication
overload detection and deterministic two-move repairs.

`test_algorithms.py` replaces the historical interactive C test programs. It
executes DBFS, DDFS, CCP, KFM, DKFM, DKFMFAST and representative multilevel
combinations through the supported CLI, validates every solution, and prevents
the objectives from exceeding the recorded final-v1.2 ceilings in
`baselines/v1.2-b14-arch0.json`. The same baseline bounds balance, elapsed time
and peak resident memory (when the host exposes process memory counters).
