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

Observing Scans follow a repeating pattern $(m_0, \dots, m_{P-1})$
(`ring_mode_pattern`). The pattern runs over the global index of *observing*
Scans,

$$
j = c\,N_S N_R + s, \qquad \mathrm{Mode} = m_{j \bmod P},
$$

so with the default pattern (B10, B01) the alternation is continuous across
Rings and cycles: since $N_S$ and $N_S N_R$ are odd, consecutive Rings and
cycles start with opposite Modes. Recycling Scans take `recycling_mode`.

The first $d$ IFGs of a Scan (`mode_switch_dead_ifgs`) are flagged when its
Mode differs from that of the previous Scan. Mode changes happen at spin phase
0, the point closest to the north ecliptic pole, where coverage is highest.

## Time budget

The fraction of time lost to the Recycling window and to mode switches is

$$
f_\mathrm{rec} = \frac{N_\mathrm{rec}}{N_\mathrm{cyc}}, \qquad
f_\mathrm{switch} = \frac{n_\mathrm{switch}\, d}{N_I N_\mathrm{cyc}},
$$

where $n_\mathrm{switch}$ counts the switching Scans per cycle. For the
baseline, alternating Modes and B00 during recycling,
$n_\mathrm{switch} = N_S N_R + 1$ (every observing Scan plus the entry into
recycling): $f_\mathrm{rec} = 20.69\,\%$, $f_\mathrm{switch} = 0.11\,\%$.

## BBIR schedule

The BBIR main-absorber set-point of cycle $c$ is $T_{c \bmod M}$, taken from
the list `setpoints_k` of length $M$. Set-points change only during Recycling
windows: the window of cycle $c$ already carries the set-point of cycle
$c + 1$ (the thermal transient is not modelled). Port $p$ sees
$T + \delta T_p$, with offsets `port_offsets_k`.
