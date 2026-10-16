# FOSSIL simulation

Simulation of the FOSSIL absolute FTS mission: observation timeline, pointing, sky signal and maps.

## Language

### Timeline

**IFG**:
One interferogram, the elementary time sample (~2 s); pointing is constant within it.
_Avoid_: sample, frame

**Scan**:
One full spin of the spacecraft: a (near) great circle traced by the boresight.
_Avoid_: rotation, sweep

**Ring**:
A set of consecutive Scans with the spin axis held fixed, so they repeat the same circle.
_Avoid_: pointing period

**Cooler cycle**:
One ADR period: the observing Rings followed by one Recycling window.
_Avoid_: chunk (chunk is an implementation term), ADR cycle

**Recycling window**:
The part of a Cooler cycle in which the ADR recycles and science data are flagged; the spacecraft keeps spinning.
_Avoid_: R, dead time

**Flag**:
A per-IFG marker of degraded or unusable data (Recycling window, Mode switch, SAA exceeded, slew).
_Avoid_: mask, cut

**Valid IFG**:
An IFG with no Flag set; only Valid IFGs enter maps.

### Pointing

**Spin axis**:
The spacecraft rotation axis, kept near the anti-Sun direction.

**Ring stepping**:
The rule that sets the Spin axis of each Ring. Baseline: per Cooler cycle the axis advances by the Sun's actual motion over the cycle, split into equal steps between Rings.
_Avoid_: repointing, slew (slew is the physical manoeuvre, not the rule)

**Precession**:
Slow motion of the Spin axis on a cone of fixed opening around the anti-Sun direction, frozen within each Ring.
_Avoid_: nutation, wobble

**Solar Aspect Angle (SAA)**:
Angle between the Spin axis and the anti-Sun direction.

### Instrument and modes

**Port**:
One of the two FTS inputs (1, 2). Port k is fed by BSMk, which routes it to Telescope k or to the BBIR.
_Avoid_: input, arm (arm is the interferometer path)

**BSM**:
Beam-switching mechanism; selects sky or BBIR for its Port.

**BBIR**:
Blackbody Internal Reference. Each Port may see its own effective BBIR spectrum.
_Avoid_: calibrator, load

**Mode**:
The pair of BSM states, written Bxy: x = BSM1, y = BSM2, 1 = sky, 0 = BBIR. The L output modulates +(I₁ − I₂), the R output −(I₁ − I₂).
- B00: both on BBIR (cross-calibration).
- B10: Port 1 on sky, Port 2 on BBIR (science).
- B01: Port 1 on BBIR, Port 2 on sky (science; opposite sky sign to B10).
- B11: both on the same sky direction.
_Avoid_: configuration, state

**Detector unit (DU)**:
One of LF-L, HF-L, LF-R, HF-R: an FTS output (L/R) split by the 300 GHz dichroic.
_Avoid_: detector, channel
