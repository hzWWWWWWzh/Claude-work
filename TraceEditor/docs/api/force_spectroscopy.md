# Table of Contents

* [trace\_editor.core.force\_spectroscopy](#trace_editor.core.force_spectroscopy)
  * [kT](#trace_editor.core.force_spectroscopy.kT)
  * [nm\_bp](#trace_editor.core.force_spectroscopy.nm_bp)
  * [nm\_base](#trace_editor.core.force_spectroscopy.nm_base)
  * [omega0](#trace_editor.core.force_spectroscopy.omega0)
  * [bilinear\_interpolation](#trace_editor.core.force_spectroscopy.bilinear_interpolation)
  * [intersection](#trace_editor.core.force_spectroscopy.intersection)
  * [median\_filter](#trace_editor.core.force_spectroscopy.median_filter)
  * [MagnetForce](#trace_editor.core.force_spectroscopy.MagnetForce)
  * [zMagnetForce](#trace_editor.core.force_spectroscopy.zMagnetForce)
  * [correct\_phase\_jump](#trace_editor.core.force_spectroscopy.correct_phase_jump)
  * [zFJC](#trace_editor.core.force_spectroscopy.zFJC)
  * [zWLC](#trace_editor.core.force_spectroscopy.zWLC)
  * [zDNA\_unzip](#trace_editor.core.force_spectroscopy.zDNA_unzip)
  * [zFiber\_old](#trace_editor.core.force_spectroscopy.zFiber_old)
  * [zFiber](#trace_editor.core.force_spectroscopy.zFiber)
  * [z\_nucleosome](#trace_editor.core.force_spectroscopy.z_nucleosome)
  * [zHMf](#trace_editor.core.force_spectroscopy.zHMf)
  * [zHMf\_new](#trace_editor.core.force_spectroscopy.zHMf_new)
  * [zHMf2](#trace_editor.core.force_spectroscopy.zHMf2)
  * [zMedianFilter](#trace_editor.core.force_spectroscopy.zMedianFilter)
  * [ExponentialDrift](#trace_editor.core.force_spectroscopy.ExponentialDrift)
  * [power\_spectrum](#trace_editor.core.force_spectroscopy.power_spectrum)

<a id="trace_editor.core.force_spectroscopy"></a>

# trace\_editor.core.force\_spectroscopy

<a id="trace_editor.core.force_spectroscopy.kT"></a>

#### kT

Thermal energy [pN*nm]

<a id="trace_editor.core.force_spectroscopy.nm_bp"></a>

#### nm\_bp

length dsDNA [nm/bp]

<a id="trace_editor.core.force_spectroscopy.nm_base"></a>

#### nm\_base

length ssDNA = [nm/base];

<a id="trace_editor.core.force_spectroscopy.omega0"></a>

#### omega0

inverted pitch of dsDNA [nm-1]

<a id="trace_editor.core.force_spectroscopy.bilinear_interpolation"></a>

#### bilinear\_interpolation

```python
def bilinear_interpolation(x, y, points)
```

Interpolate (x,y) from values associated with four points.

The four points are a list of four triplets:  (x, y, value).
The four points can be in any order.  They should form a rectangle.

    >>> bilinear_interpolation(12, 5.5,
    ...                        [(10, 4, 100),
    ...                         (20, 4, 200),
    ...                         (10, 6, 150),
    ...                         (20, 6, 300)])
    165.0

<a id="trace_editor.core.force_spectroscopy.intersection"></a>

#### intersection

```python
def intersection(curve1, curve2, xrange=[-np.inf, np.inf])
```

Find the intercept of two curves, given by the different x data

<a id="trace_editor.core.force_spectroscopy.median_filter"></a>

#### median\_filter

```python
def median_filter(x, k)
```

Apply a length-k median filter to a 1D array x.
Boundaries are extended by repeating endpoints.

<a id="trace_editor.core.force_spectroscopy.MagnetForce"></a>

#### MagnetForce

```python
def MagnetForce(x, L1=1.4, L2=0.8, f0=0.02, fmax=85.0, alpha=0.7, dt=0)
```

Calculate the force from the magnet position.

**Arguments**:

- `x`: Shift (mm) #
- `L1`: L1 (mm) #
- `L2`: L2 (mm) #
- `f0`: Fmin (pN) #
- `fmax`: Fmax (pN) #
- `alpha`: alpha (-) #
- `dt`: offset (pts) #

**Returns**:

F (pN) #

<a id="trace_editor.core.force_spectroscopy.zMagnetForce"></a>

#### zMagnetForce

```python
def zMagnetForce(x,
                 L1=1.4,
                 L2=0.8,
                 f0=0.02,
                 fmax=85.0,
                 alpha=0.7,
                 L=1,
                 dt=0,
                 z_bead=0)
```

Calculate the force from the magnet position.

**Arguments**:

- `x`: Shift (mm) #
- `L1`: L1 (mm) #
- `L2`: L2 (mm) #
- `f0`: Fmin (pN) #
- `fmax`: Fmax (pN) #
- `alpha`: alpha (-) #
- `L`: L (bp) #
- `dt`: offset (pts) #
- `z_bead`: Zbead (um) #

**Returns**:

Z (um) #

<a id="trace_editor.core.force_spectroscopy.correct_phase_jump"></a>

#### correct\_phase\_jump

```python
def correct_phase_jump(x, y, z, xy_margin=np.inf, z_margin=5)
```

Correct the z-jumps caused by limited range of phase calibration

**Arguments**:

- `x`: x timetrace
- `y`: y timetrace
- `z`: z timetrace
- `xy_margin`: maximal x/y jump (um)
- `z_margin`: minimal z jump (um)

**Returns**:

correced z timetrace

<a id="trace_editor.core.force_spectroscopy.zFJC"></a>

#### zFJC

```python
def zFJC(f, k=0.3, z0=0.0, L=25.0, b=None)
```

Compute force induced extension of a Freely Jointed Chain

**Arguments**:

- `f`: force (pN)
- `k`: stiffness (pN/nm), computed from b if not supplied
- `z0`: rest extension (nm)
- `L`: contour length (nm)
- `b`: Kuhn length (nm)

**Returns**:

extension (nm)

<a id="trace_editor.core.force_spectroscopy.zWLC"></a>

#### zWLC

```python
def zWLC(x, L=1, A=50, S=1000, z_bead=0)
```

extendable Worm-Like Chain model

**Arguments**:

- `x`: F (pN) # Force
- `A`: A (nm) # persistence length
- `L`: L (bp) # Contour length of DNA
- `S`: S (pN) # Stretch modulus
- `tw`: twist (turns) # twist

**Returns**:

Z (um) # extension of DNA

<a id="trace_editor.core.force_spectroscopy.zDNA_unzip"></a>

#### zDNA\_unzip

```python
def zDNA_unzip(f, L=1, Lds_bp=1, Lss0_bp=0)
```

two state unzipping model

**Arguments**:

- `f`: force (pN)
- `L`: contourlength (bp)
- `Lds_bp`: 
- `Lss0_bp`: 

**Returns**:

DNA extension (nm)

<a id="trace_editor.core.force_spectroscopy.zFiber_old"></a>

#### zFiber\_old

```python
def zFiber_old(x,
               L=6000,
               A=50,
               S=1000,
               n4=0,
               n=0,
               k=0.1,
               NLD=1.2,
               dg1=17,
               dg2=5,
               dg3=50,
               NRL=197,
               L1=89,
               L2=80,
               d=1,
               z_bead=0)
```

Compute force-dependent extension of chromatin fiber according to Meng2019NAR

**Arguments**:

- `x`: F (pN) # Force
- `L`: L (bp) `contour` length (bp)
- `A`: A (nm) # persistence length
- `S`: S (pN) # Stretch modulus
- `n4`: n4 # number of tetrasomes
- `n`: n # number of nucleosomes (incl. tetrasomes)
- `k`: k (pN_nm) # stiffness (pN/nm)
- `NLD`: NLD (nm) `Nucleosome` Line Density (nm)
- `dg1`: G1 (kT) `s` tacking energy (kT)
- `dg2`: G2 (kT) # partial unwrapping energy (kT) of (l2-L1) basepairs
- `dg3`: G3 (kT) # full unwrapping energy (kT) of last L2 basepairs
- `NRL`: NRL (bp) # Nucleosome Repeat Length (bp)
- `L1`: L1 (bp) # wrapped amount of DNA (bp) in unstacked conformation
- `L2`: L2 (bp) # wrapped amount of DNA (bp) in singly wrapped conformation
- `d`: degeneracy # degeneracy (0 = cooperative, 1 = degenerate)
- `z_bead`: Zbead (um) # Height offset

**Returns**:

Z (um) # fiber extension

<a id="trace_editor.core.force_spectroscopy.zFiber"></a>

#### zFiber

```python
def zFiber(x,
           L=6000,
           A=50,
           S=1000,
           n0=0,
           n4=0,
           n=0,
           k=0.1,
           NLD=1.2,
           dg_unwrap=12,
           dg1=17,
           dg2=5,
           dg3=50,
           NRL=197,
           L0=147,
           L1=89,
           L2=80,
           d=1,
           z_bead=0)
```

Compute force-dependent extension of chromatin fiber according to Meng2019NAR

**Arguments**:

- `x`: F (pN) # Force
- `L`: L (bp) `contour` length (bp)
- `A`: A (nm) # persistence length
- `S`: S (pN) # Stretch modulus
- `n0`: n0 # number of nucleosomes that don not stack into a fiber
- `n4`: n4 # number of tetrasomes
- `n`: n # number of nucleosomes (incl. tetrasomes)
- `k`: k (pN_nm) # stiffness (pN/nm)
- `NLD`: NLD (nm) `Nucleosome` Line Density (nm)
- `dg_unwrap`: Gwrap (kT) # wrapping energy for a free (non-stacked) nucleosome (kT)
- `dg1`: G1 (kT) `s` tacking energy (kT)
- `dg2`: G2 (kT) # partial unwrapping energy (kT) of (l2-L1) basepairs
- `dg3`: G3 (kT) # full unwrapping energy (kT) of last L2 basepairs
- `NRL`: NRL (bp) # Nucleosome Repeat Length (bp)
- `L0`: L0 (bp) # wrapped amount of DNA (bp) in a wrapped, non-stacked nucleosome
- `L1`: L1 (bp) # wrapped amount of DNA (bp) in unstacked conformation
- `L2`: L2 (bp) # wrapped amount of DNA (bp) in singly wrapped conformation
- `d`: degeneracy # degeneracy (0 = cooperative, 1 = degenerate)
- `z_bead`: Zbead (um) # Height offset

**Returns**:

Z (um) # fiber extension

<a id="trace_editor.core.force_spectroscopy.z_nucleosome"></a>

#### z\_nucleosome

```python
def z_nucleosome(x,
                 L=1000,
                 A=50,
                 S=1000,
                 NRL=197,
                 dg_unwrap=12,
                 dg2=5,
                 dg3=50,
                 L0=147,
                 L1=89,
                 L2=80)
```

Compute force-dependent extension of a single nucleosome with multiple states (fully wrapped, partially unwrapped, single wrap, fully unwrapped)

**Arguments**:

- `x`: F (pN) # Force
- `L`: L (bp) `contour` length (bp)
- `A`: A (nm) # persistence length
- `S`: S (pN) # Stretch modulus
- `n0`: n0 # number of nucleosomes that don not stack into a fiber
- `n4`: n4 # number of tetrasomes
- `n`: n # number of nucleosomes (incl. tetrasomes)
- `k`: k (pN_nm) # stiffness (pN/nm)
- `NLD`: NLD (nm) `Nucleosome` Line Density (nm)
- `dg_unstacked`: Gwrap (kT) # wrapping energy for a free (non-stacked) nucleosome (kT)
- `dg1`: G1 (kT) `s` tacking energy (kT)
- `dg2`: G2 (kT) # partial unwrapping energy (kT) of (l2-L1) basepairs
- `dg3`: G3 (kT) # full unwrapping energy (kT) of last L2 basepairs
- `NRL`: NRL (bp) # Nucleosome Repeat Length (bp)
- `L0`: L0 (bp) # wrapped amount of DNA (bp) in a wrapped, non-stacked nucleosome
- `L1`: L1 (bp) # wrapped amount of DNA (bp) in unstacked conformation
- `L2`: L2 (bp) # wrapped amount of DNA (bp) in singly wrapped conformation
- `d`: degeneracy # degeneracy (0 = cooperative, 1 = degenerate)
- `z_bead`: Zbead (um) # Height offset

**Returns**:

Z (um) # fiber extension

<a id="trace_editor.core.force_spectroscopy.zHMf"></a>

#### zHMf

```python
def zHMf(x,
         L=3646,
         A=50,
         S=1000,
         k=0.5,
         z_dimer=4,
         L_dimer=30,
         angle=10,
         n_dimer=100,
         n_gap=0,
         g1=1.9,
         g2=4.11,
         Z0=0)
```

Calculate extension of HMf fiber according to Henneman2020NAR

**Arguments**:

- `x`: F (pN) # Force
- `L`: L (bp) # [1, inf]  contour length
- `A`: A (nm) # [1, 500] persistence length
- `S`: S (pN) # [1, 5000]Stretch modulus
- `k`: k (pN_nm) [0.01, inf] # Stiffness of stacked fiber
- `z_dimer`: Zdimer (nm) # [0.01, 10] maximal extension of FJC of stacked fiber
- `L_dimer`: Ldimer (bp) # [1, 50] Footprint of dimer
- `angle`: angle (deg) # [0, 90]  deflection angle of DNA induced by dimer
- `n_dimer`: Ndimer (-) # [0, 300] number of bound proteins
- `n_gap`: Ngap (-) # [0, 100] Number of defects in hypernucleosome
- `g1`: G1 (kT) # [0, 10] Stacking energy
- `g2`: G2 (kT) # [0, 10] Wrapping energy
- `Z0`: Zbead (um) # [-2, 0] Offset due to bead atachement

**Returns**:

Z (um) # extension of the fiber

<a id="trace_editor.core.force_spectroscopy.zHMf_new"></a>

#### zHMf\_new

```python
def zHMf_new(x,
             L=3646,
             A=50,
             S=1000,
             k=0.5,
             z_dimer=4,
             L_dimer=30,
             angle=10,
             n_dimer=100,
             n_gap=0,
             g1=1.9,
             g2=4.11,
             Z0=0)
```

Calculate extension of HMf fiber according to Henneman2020NAR

**Arguments**:

- `x`: F (pN) # Force
- `L`: L (bp) # [1, inf]  contour length
- `A`: A (nm) # [1, 500] persistence length
- `S`: S (pN) # [1, 5000]Stretch modulus
- `k`: k (pN_nm) [0.01, inf] # Stiffness of stacked fiber
- `z_dimer`: Zdimer (nm) # [0.01, 10] maximal extension of FJC of stacked fiber
- `L_dimer`: Ldimer (bp) # [1, 50] Footprint of dimer
- `angle`: angle (deg) # [0, 90]  deflection angle of DNA induced by dimer
- `n_dimer`: Ndimer (-) # [0, 300] number of bound proteins
- `n_gap`: Ngap (-) # [0, 100] Number of defects in hypernucleosome
- `g1`: G1 (kT) # [0, 10] Stacking energy
- `g2`: G2 (kT) # [0, 10] Wrapping energy
- `Z0`: Zbead (um) # [-2, 0] Offset due to bead atachement

**Returns**:

Z (um) # extension of the fiber

<a id="trace_editor.core.force_spectroscopy.zHMf2"></a>

#### zHMf2

```python
def zHMf2(x,
          L=3646,
          A=50,
          S=1000,
          k=0.5,
          z_dimer=4,
          L_dimer=30,
          angle=10,
          n_dimer=100,
          n_gap=0,
          g1=1.9,
          g2=4.11,
          Z0=0)
```

Calculate extension of HMf fiber according to Henneman2020NAR.

This function is a wrapper to allow for non-integer n_dimer and n_gap values.

**Arguments**:

- `x`: F (pN) # Force
- `L`: L (bp) # [1, inf]  contour length
- `A`: A (nm) # [1, 500] persistence length
- `S`: S (pN) # [1, 5000] Stretch modulus
- `k`: k (pN_nm) # [0.01, inf] Stiffness of stacked fiber
- `z_dimer`: Zdimer (nm) # [0.01, 10] maximal extension of FJC of stacked fiber
- `L_dimer`: Ldimer (bp) # [1, 50] Footprint of dimer
- `angle`: angle (deg) # [0, 90]  deflection angle of DNA induced by dimer
- `n_dimer`: Ndimer (-) # [0, 300] number of bound proteins
- `n_gap`: Ngap (-) # [0, 100] Number of defects in hypernucleosome
- `g1`: G1 (kT) # [0, 10] Stacking energy
- `g2`: G2 (kT) # [0, 10] Wrapping energy
- `Z0`: Zbead (um) # [-2, 0] Offset due to bead attachment

**Returns**:

Z (um) # extension of the fiber

<a id="trace_editor.core.force_spectroscopy.zMedianFilter"></a>

#### zMedianFilter

```python
def zMedianFilter(x, w=1.0, t=None)
```

Calculate median filter.

**Arguments**:

- `x`: z (um) # Heigth of tether
- `w`: w (s) # Width of median filter
- `t`: t (s) # Time trace to determine filter window

**Returns**:

Z (um) # Median filtered height of tether

<a id="trace_editor.core.force_spectroscopy.ExponentialDrift"></a>

#### ExponentialDrift

```python
def ExponentialDrift(x, Z0=0, dZ_dt=0, tau=np.inf)
```

Calculate drift, assuming an exponential decay

**Arguments**:

- `x`: Time (s)  # Time
- `Z0`: Z0 (um) # offset
- `dZ_dt`: dZ_dt (nm_s) # drift
- `tau`: tau (s) # decay time

**Returns**:

Z (um) #

<a id="trace_editor.core.force_spectroscopy.power_spectrum"></a>

#### power\_spectrum

```python
def power_spectrum(x, r_um=1.4)
```

Placeholder to prepare GUI for asking radius of bead

**Arguments**:

- `x`: Z (um) # Trace to calculate power spectrum
- `r_um`: r (um) # Radius of the bead

**Returns**:

None

