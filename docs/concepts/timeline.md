# Observation timeline

## Hierarchy

The mission is a sequence of Cooler cycles. Each cycle contains $N_R$ Rings
followed by a Recycling window; a Ring contains $N_S$ Scans, a Scan $N_I$ IFGs
of duration $\Delta t$. The Recycling window lasts $N_\mathrm{rec}$ whole
Scans, so that every Ring starts at the same spin phase. With

$$
N_\mathrm{cyc} = N_S N_R + N_\mathrm{rec}
$$

Scans per cycle, the derived durations are

$$
T_\mathrm{scan} = N_I\,\Delta t, \qquad
T_\mathrm{ring} = N_S\,T_\mathrm{scan}, \qquad
T_\mathrm{cycle} = N_\mathrm{cyc}\,T_\mathrm{scan}.
$$

For the baseline ($\Delta t = 2$ s, $N_I = 720$, $N_S = 23$, $N_R = 3$,
$N_\mathrm{rec} = 18$): $T_\mathrm{scan} = 1440$ s, $T_\mathrm{ring} = 9.2$ h,
$T_\mathrm{cycle} = 34.8$ h, about 881 cycles in 3.5 years.

## Indices

Every quantity is a function of the global IFG index $g$ (ADR 0001). The IFG
spans $[g\,\Delta t, (g+1)\,\Delta t)$ and its pointing is evaluated at the
mid-time $t_g = (g + \tfrac12)\,\Delta t$. The Cooler cycle, the local index
and the position within the cycle are

$$
c = \left\lfloor \frac{g}{N_I N_\mathrm{cyc}} \right\rfloor, \qquad
\ell = g - c\,N_I N_\mathrm{cyc}, \qquad
s = \left\lfloor \frac{\ell}{N_I} \right\rfloor, \qquad
k = \ell - s\,N_I,
$$

with $s$ the Scan in the cycle and $k$ the IFG in the Scan. The Ring in the
cycle is $r = \lfloor s / N_S \rfloor$; the IFG is in the Recycling window if
$r \ge N_R$, in which case the pointing keeps the last Ring, $r = N_R - 1$.
The global Scan and Ring indices are $c\,N_\mathrm{cyc} + s$ and
$c\,N_R + r$ (Ring index $-1$ during recycling).

## Modes

The observing time (Recycling windows excluded) is split into Mode blocks of
duration $\tau_\mathrm{block}$ (`mode_switch_interval_s`, a multiple of
$\Delta t$, i.e. $L = \tau_\mathrm{block} / \Delta t$ IFGs). The block clock
stops during Recycling windows, so a block that does not end before one
resumes after it in the same Mode. A repeating pattern
$(m_0, \dots, m_{P-1})$ (`mode_pattern`) gives the Mode of each block:

$$
j = \left\lfloor \frac{c\,N_S N_R N_I + k_c}{L} \right\rfloor,
\qquad \mathrm{Mode} = m_{j \bmod P},
$$

with $k_c$ the IFG index in the cycle. Recycling IFGs take `recycling_mode`.
With the baseline Scan (1440 s) and Ring (23 Scans, 33120 s):

| Switch | `mode_switch_interval_s` |
|---|---|
| at the end of every Scan (spin phase 0, north ecliptic pole) | 1440 |
| Ring split into $n$ parts | 33120 / $n$ |
| one Mode per Ring | 33120 |
| one Mode per cycle | 99360 (3 Rings) |

Half Rings (16560 s, $n = 2$) switch at spin phase π of the middle Scan, the
south ecliptic pole, since $N_S = 23$ is odd. The baseline uses them with the
pattern (B10, B01, B01, B10): every Ring is observed in both Modes
(11.5 Scans each) and the Mode never changes between Rings, so each Mode
covers the full sky in six months. Blocks longer than a Ring leave each Mode a
subset of the Rings and a sparser coverage (one cycle per Mode: 79 % of the
sky per Mode at nside 64 after six months). Reversing the pattern swaps the
Modes.

After a Mode switch inside an observing period the first
$d = \lceil \tau / \Delta t \rceil$ IFGs are flagged, with $\tau$ =
`mode_switch_dead_time_s` (baseline 300 s, $d = 150$; 0 = instantaneous
switch). The start of an observing period (mission start, end of Recycling)
is not a switch: the BSMs are set during the Recycling window.

## Time budget

The fraction of time lost to the Recycling window and to mode switches is

$$
f_\mathrm{rec} = \frac{N_\mathrm{rec}}{N_\mathrm{cyc}}, \qquad
f_\mathrm{switch} = \frac{n_\mathrm{switch}\, d}{N_I N_\mathrm{cyc}},
$$

where $n_\mathrm{switch}$ counts the switches inside the observing time of a
cycle. For the baseline (half Rings, three switches per cycle):
$f_\mathrm{rec} = 20.69\,\%$, $f_\mathrm{switch} = 0.72\,\%$. Alternating Scan
by Scan with the same dead time would lose $68 \cdot 150 / 62640 = 16.3\,\%$,
always on the same arc after the north pole.

## BBIR schedule

The BBIR main-absorber set-point of cycle $c$ is $T_{c \bmod M}$, taken from
the list `setpoints_k` of length $M$. Set-points change only during Recycling
windows: the window of cycle $c$ already carries the set-point of cycle
$c + 1$ (the thermal transient is not modelled). Port $p$ sees
$T + \delta T_p$, with offsets `port_offsets_k`.
