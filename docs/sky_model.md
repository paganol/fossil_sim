# Sky model

The sky is the sum of three parts ({doc}`adr/0002-isotropic-backgrounds-in-time-domain`):

- the **Isotropic backgrounds**, isotropic in the CMB rest frame, evaluated
  for every IFG with the exact Boost to the spacecraft frame;
- the **Static sky**, the time-independent anisotropies (Galactic
  foregrounds, CMB anisotropies), defined in the barycentric frame;
- the **zodiacal light**, emission of the interplanetary dust, which depends
  on the position of the observer and therefore on time.

```{toctree}
:maxdepth: 2

concepts/boost
concepts/isotropic_backgrounds
concepts/static_sky
concepts/zodiacal_light
```
