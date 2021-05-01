- [Musical Systems](#org46f1ddc)
  - [Musical Systems](#org1d8a203)
  - [Set Classes](#org2b41158)
  - [Forte Table](#org5a8e007)
  - [Class Elements](#orgaa87e6a)
  - [Webern Matrices](#org4fb2f0f)



<a id="org46f1ddc"></a>

# Musical Systems

We provide here the routines to work with new musical systems as described in FIXME. Note that these routines run in [Sage](https://www.sagemath.org/), not straight [Python](https://www.python.org/), as it makes it easier to deal with permutation groups.

Sage is free and open source, and can be used in any platform. If one wants to use it without having to install it, [Cocalc](https://cocalc.com/) allows you to run it online, but it is necessary to create an account.

The file `MusicalSystem.sage` provides all necessary routines. It provides FIXME classes:

-   `MusicalSystem`: allows you to create new musical systems and manipulate some of its elements.
-   `SetClass`: to deal with set classes in new musical systems.
-   `ForteTable:` to produce Forte tables for new musical systems.
-   `ClassElements:` to see the set of all set classes and their elements of a fixed size in a musical system.
-   `Webern Matrix:` to produce examples of Webern matrices in new systems.

We describe the functionality of each below.


<a id="org1d8a203"></a>

## Musical Systems

First we need to load the file:

```sage
load('MusicalSystem.sage')
```

One can create the standard (equal temperament) system by calling `MusicalSystem` with nor arguments:

```sage
SMS = MusicalSystem()
SMS
```

    MusicalSystem(rho=(0,1,2,3,4,5,6,7,8,9,10,11), phi=(1,11)(2,10)(3,9)(4,8)(5,7)): Musical System given by rotation (0,1,2,3,4,5,6,7,8,9,10,11) and reflection (1,11)(2,10)(3,9)(4,8)(5,7)

To create a new system, let's first create the symmetric group on integers modulo 12 (or the set {0, 1, 2, 3, &#x2026; , 11}) in Sage:

```sage
S12 = SymmetricGroup(range(12))
```

**IMPORTANT:** Note that Sage, like some group theorists, compose functions from left to right, instead of the usual right to left used in all mathematics. So, for Sage, `(sigma * tau)(x)` means `tau(sigma(x))`, and not the usual `sigma(tau(x))`. On the other hand, */I will always refer to composition the usual way*, i.e., from right to left.

Let's also create the standard rotation/transposition and reflection/inversion:

```sage
rho = S12([(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)])
phi = S12([(1, 11), (2, 10), (3, 9), (4, 8), (5, 7)])
```

We can obtain the corresponding rotation and reflection from a system with `.rho` and `.phi`:

```sage
SMS.rho == rho, SMS.phi == phi
```

    (True, True)

Let's create a new system given by the rotation `(0, 4, 11, 3, 7, 2, 6, 10, 5, 9, 1, 8)`:

```sage
rho1 = S12([(0, 4, 11, 3, 7, 2, 6, 10, 5, 9, 1, 8)])
phi1 = phi  # same as the standard system
MS1 = MusicalSystem(rho1, phi1)
MS1
```

    MusicalSystem(rho=(0,4,11,3,7,2,6,10,5,9,1,8), phi=(1,11)(2,10)(3,9)(4,8)(5,7)): Musical System given by rotation (0,4,11,3,7,2,6,10,5,9,1,8) and reflection (1,11)(2,10)(3,9)(4,8)(5,7)

As observed in the paper, every new system can be given by a conjugation of the original dihedral group by some permutation. We can use a permutation to create a new system with `MusicalSystem.from_conj`:

```sage
sigma = S12([(0,1,8,7,6,9), (2,5,3,11,4)])  # some permutation
MS2 = MusicalSystem.from_conj(sigma)
MS2
```

    MusicalSystem(rho=(0,10,4,1,8,5,11,2,3,9,6,7), phi=(1,9)(2,5)(3,8)(4,6)(7,10)): Musical System given by rotation (0,10,4,1,8,5,11,2,3,9,6,7) and reflection (1,9)(2,5)(3,8)(4,6)(7,10)

Let's save the corresponding rotation and reflection:

```sage
rho2, phi2 = MS2.rho, MS2.phi
```

We can measure intervals in any of the systems. Let's see what are the intervals between the pitch classes 3 and 10 in these systems:

```sage
SMS.interval(3, 10), MS1.interval(3, 10), MS2.interval(3, 10)
```

    (7, 4, 5)

Let's check:

```sage
(rho^7)(3) == 10, (rho1^4)(3) == 10, (rho2^5)(3) == 10
```

    (True, True, True)

We can also ask for the interval class:

```sage
SMS.interval_class(3, 10), MS1.interval_class(3, 10), MS2.interval_class(3, 10)
```

    (5, 4, 5)

Of course, we can change the order and obtain the same result:

```sage
SMS.interval_class(10, 3), MS1.interval_class(10, 3), MS2.interval_class(10, 3)
```

    (5, 4, 5)

Given a pitch class set (given as a list), we can ask for the normal form in the corresponding system:

```sage
pc_set = [3, 5, 10, 11]
SMS.normal_form(pc_set), MS1.normal_form(pc_set), MS2.normal_form(pc_set)
```

    ([10, 11, 3, 5], [11, 3, 10, 5], [10, 5, 11, 3])

We can also find what are the symmetries of a pitch class set for each system:

```sage
pc_set = [0, 4, 8]
SMS.symmetries(pc_set)
```

    [(0, 0), (4, 0), (8, 0), (0, 1), (4, 1), (8, 1)]

The output tells us that the pitch class set `[0, 4, 8]` is mapped onto itself by `rho^0 * phi^0` (the identity), `rho^4 * phi^0 = \rho^4`, `rho^8 phi^0 = rho^8`, `rho^0 * phi^1 = \phi`, `rho^4 * phi^1 = rho^4 * phi`, and `rho^8 * phi^1 = rho^8 * phi`.

Let's check if this same set has any symmetries in the other systems:

```sage
MS1.symmetries(pc_set)
```

    [(0, 0), (0, 1)]

So, in the systems `MS1`, only the identity and `phi1` preserve the set.

```sage
MS2.symmetries(pc_set)
```

    [(0, 0), (4, 1)]

In the system `MS2`, only the identity and `rho_2^4 * phi_2` preserve the set.

We can also compute class sums of pitch class sets in any system:

```sage
pc_set = [0, 3, 8, 10]
SMS.class_sum(pc_set), MS1.class_sum(pc_set), MS2.class_sum(pc_set)
```

    (9, 9, 1)

We also have functions to compute symmetries of pitch class sets. To apply the rotation of the system to all pitch classes in the set:

```sage
pc_set = [0, 1, 6]
SMS.vrho(pc_set), MS1.vrho(pc_set), MS2.vrho(pc_set)
```

    ([1, 2, 7], [10, 8, 4], [7, 10, 8])

By default, the result is in *normal form*. Given the option `sort='none'`, it give the result with the outputs ordered according to input:

```sage
pc_set = [0, 1, 6]
SMS.vrho(pc_set, sort='none'), MS1.vrho(pc_set, sort='none'), MS2.vrho(pc_set, sort='none')
```

    ([1, 2, 7], [4, 8, 10], [10, 8, 7])

We can also ask to order in increasing numerical order with the option `sort='sort'`:

```sage
pc_set = [0, 1, 6]
SMS.vrho(pc_set, sort='sort'), MS1.vrho(pc_set, sort='sort'), MS2.vrho(pc_set, sort='sort')
```

    ([1, 2, 7], [4, 8, 10], [7, 8, 10])

Finally, we can sort according to the interval to the pitch class 0 in the system with `sort='rho`:

```sage
pc_set = [0, 1, 6]
SMS.vrho(pc_set, sort='rho'), MS1.vrho(pc_set, sort='rho'), MS2.vrho(pc_set, sort='rho')
```

    ([1, 2, 7], [4, 10, 8], [10, 8, 7])

We can also compute powers of the rotation with optional argument `power`, for instance, let's take the 4th power:

```sage
pc_set = [0, 1, 6]
SMS.vrho(pc_set, power=4), MS1.vrho(pc_set, power=4), MS2.vrho(pc_set, power=4)
```

    ([4, 5, 10], [1, 11, 7], [4, 8, 2])

We also have the corresponding method `.vphi` for the reflection:

```sage
pc_set = [1, 2, 6]
SMS.vphi(pc_set, sort='none'), MS1.vphi(pc_set, sort='none'), MS2.vphi(pc_set, sort='none')
```

    ([11, 10, 6], [11, 10, 6], [9, 5, 4])

(Remember that `phi_1 = phi`.)

To mix those, we can call `.vmap`. The first argument is the power of the rotation and the second is the power of the reflection, while the third is the pitch class set. Note that the power of the reflection is computed **first** and the power of the rotation **second**. It has the same sorting options as `.vhro` and `.vphi`:

```sage
pc_set = [1, 2, 6]
SMS.vmap(4, 1, pc_set, sort='none'), MS1.vmap(4, 1, pc_set, sort='none'), MS2.vmap(4  , 1, pc_set, sort='none')
```

    ([3, 2, 10], [6, 8, 1], [10, 9, 11])

Finally, given two pitch class sets, we can ask what symmetries of the system maps one to the other:

```sage
pc1 = [0, 3, 5]
pc2 = [2, 4, 7]
MS2.symmetry_maps(pc1, pc2)
```

    [(7, 1)]

This means that `rho2^7 * phi2` maps `[0, 3, 5]` to `[2, 4, 7]` in the system `MS2`:

```sage
MS2.vmap(7, 1, pc1)
```

    [2, 7, 4]

(Note that order does not matter.)


<a id="org2b41158"></a>

## Set Classes

We can also create set classes in different systems. Still with the systems `SMS`, `MS1`, and `MS2` above, we can create set classes with:

```sage
pc_set = [0, 4, 8]
sc, sc1, sc2 = SetClass(pc_set, MS=SMS), SetClass(pc_set, MS=MS1),  SetClass(pc_set, MS=MS2)
sc, sc1, sc2
```

    (SetClass([0, 4, 8]), with rho = (0,1,2,3,4,5,6,7,8,9,10,11), phi = (1,11)(2,10)(3,9)(4,8)(5,7).  Prime form: (048),
     SetClass([0, 4, 8]), with rho = (0,4,11,3,7,2,6,10,5,9,1,8), phi = (1,11)(2,10)(3,9)(4,8)(5,7).  Prime form: (04e),
     SetClass([0, 4, 8]), with rho = (0,10,4,1,8,5,11,2,3,9,6,7), phi = (1,9)(2,5)(3,8)(4,6)(7,10).  Prime form: (048))

We can ask for the number of internal symmetries:

```sage
sc.nsym(), sc1.nsym(), sc2.nsym()
```

    ([3, 3], [1, 1], [1, 1])

The first element is the number of tranpositional symmetries (including the identity) and the second is the number of reflexive symmetries. We can actually see what the symmetries are with `.symmetries`:

```sage
sc.symmetries(), sc1.symmetries(), sc2.symmetries()
```

    ([(0, 0), (4, 0), (8, 0), (0, 1), (4, 1), (8, 1)],
     [(0, 0), (2, 1)],
     [(0, 0), (4, 1)])

We can also ask for the complement of a set class:

```sage
sc2, sc2.complement()
```

    (SetClass([0, 4, 8]), with rho = (0,10,4,1,8,5,11,2,3,9,6,7), phi = (1,9)(2,5)(3,8)(4,6)(7,10).  Prime form: (048),
     SetClass([1, 2, 3, 5, 6, 7, 9, 10, 11]), with rho = (0,10,4,1,8,5,11,2,3,9,6,7), phi = (1,9)(2,5)(3,8)(4,6)(7,10).  Prime form: (0t4185e36))

And we can ask for class sums:

```sage
sc1.class_sum()
```

    3

Finally, we can ask for all sets in a set class. For instance:

```sage
sc1.elements()
```

    [(0, 4, 11),
     (1, 8, 0),
     (2, 6, 10),
     (3, 7, 2),
     (4, 11, 3),
     (5, 9, 1),
     (6, 10, 5),
     (7, 2, 6),
     (8, 0, 4),
     (9, 1, 8),
     (10, 5, 9),
     (11, 3, 7)]

Note that the results are in normal form.

As another example, if I want to know all tetrachors that can be mapped to `[10, 4, 5]` in `MS2`, we can do

```sage
SetClass([10, 4, 5], MS=MS2).elements()
```

    [(0, 1, 8),
     (0, 10, 8),
     (1, 8, 2),
     (1, 11, 2),
     (2, 3, 7),
     (2, 6, 7),
     (3, 7, 0),
     (3, 9, 0),
     (4, 1, 11),
     (4, 5, 11),
     (5, 3, 9),
     (5, 11, 9),
     (6, 7, 4),
     (6, 10, 4),
     (7, 0, 1),
     (7, 4, 1),
     (8, 2, 3),
     (8, 5, 3),
     (9, 0, 10),
     (9, 6, 10),
     (10, 4, 5),
     (10, 8, 5),
     (11, 2, 6),
     (11, 9, 6)]


<a id="org5a8e007"></a>

## Forte Table

We can also ask for the *Forte Table* for a system. For instance, here is the Forte Table for trichords and nonachors in the standard system:

```sage
ft = ForteTable(3)
print(ft)
```

    (012)  [ 2 , 1 , 0 , 0 , 0 , 0]   1,  1  [ 8 , 7 , 6 , 6 , 6 , 3]  (012345678)
    (013)  [ 1 , 1 , 1 , 0 , 0 , 0]   1,  0  [ 7 , 7 , 7 , 6 , 6 , 3]  (012345679)
    (014)  [ 1 , 0 , 1 , 1 , 0 , 0]   1,  0  [ 7 , 6 , 7 , 7 , 6 , 3]  (012345689)
    (015)  [ 1 , 0 , 0 , 1 , 1 , 0]   1,  0  [ 7 , 6 , 6 , 7 , 7 , 3]  (012345789)
    (016)  [ 1 , 0 , 0 , 0 , 1 , 1]   1,  0  [ 7 , 6 , 6 , 6 , 7 , 4]  (012346789)
    (024)  [ 0 , 2 , 0 , 1 , 0 , 0]   1,  1  [ 6 , 8 , 6 , 7 , 6 , 3]  (01234568t)
    (025)  [ 0 , 1 , 1 , 0 , 1 , 0]   1,  0  [ 6 , 7 , 7 , 6 , 7 , 3]  (01234578t)
    (026)  [ 0 , 1 , 0 , 1 , 0 , 1]   1,  0  [ 6 , 7 , 6 , 7 , 6 , 4]  (01234678t)
    (027)  [ 0 , 1 , 0 , 0 , 2 , 0]   1,  1  [ 6 , 7 , 6 , 6 , 8 , 3]  (01235678t)
    (036)  [ 0 , 0 , 2 , 0 , 0 , 1]   1,  1  [ 6 , 6 , 8 , 6 , 6 , 4]  (01234679t)
    (037)  [ 0 , 0 , 1 , 1 , 1 , 0]   1,  0  [ 6 , 6 , 7 , 7 , 7 , 3]  (01235679t)
    (048)  [ 0 , 0 , 0 , 3 , 0 , 0]   3,  3  [ 6 , 6 , 6 , 9 , 6 , 3]  (01245689t)

The first and last column have the set classes, the second and second to last have interval vectors, and the two middle columns have the number of transpositional and inversive symmetries, respectively. Note that we do not give the traditional names associated to the rows.

Let's see it for a different system, say `MS2`, now with tetrachords and octachords:

```sage
ft2 = ForteTable(4, MS=MS2)
print(ft2)
```

    (0153)  [ 0 , 1 , 2 , 1 , 2 , 0]   1,  1  [ 4 , 5 , 6 , 5 , 6 , 2]  (0t185236)
    (0182)  [ 1 , 0 , 2 , 2 , 1 , 0]   1,  1  [ 5 , 4 , 6 , 6 , 5 , 2]  (0t185e39)
    (01e9)  [ 0 , 0 , 4 , 0 , 0 , 2]   4,  4  [ 4 , 4 , 8 , 4 , 4 , 4]  (0t18e296)
    (0412)  [ 1 , 1 , 1 , 1 , 2 , 0]   1,  0  [ 5 , 5 , 5 , 5 , 6 , 2]  (0t485e29)
    (0415)  [ 1 , 2 , 2 , 0 , 1 , 0]   1,  1  [ 5 , 6 , 6 , 4 , 5 , 2]  (04185e29)
    (041e)  [ 1 , 1 , 2 , 1 , 0 , 1]   1,  0  [ 5 , 5 , 6 , 5 , 4 , 3]  (0t185e29)
    (0452)  [ 0 , 2 , 1 , 0 , 3 , 0]   1,  1  [ 4 , 6 , 5 , 4 , 7 , 2]  (0t415236)
    (0453)  [ 0 , 1 , 2 , 1 , 1 , 1]   1,  0  [ 4 , 5 , 6 , 5 , 5 , 3]  (0t485236)
    (0482)  [ 0 , 2 , 1 , 1 , 2 , 0]   1,  0  [ 4 , 6 , 5 , 5 , 6 , 2]  (0t415e36)
    (0483)  [ 0 , 2 , 0 , 3 , 0 , 1]   1,  1  [ 4 , 6 , 4 , 7 , 4 , 3]  (0t485e36)
    (048e)  [ 0 , 3 , 0 , 2 , 0 , 1]   1,  1  [ 4 , 7 , 4 , 6 , 4 , 3]  (0t418e36)
    (04e3)  [ 0 , 2 , 0 , 2 , 0 , 2]   2,  2  [ 4 , 6 , 4 , 6 , 4 , 4]  (0t48e236)
    (0t12)  [ 1 , 1 , 1 , 1 , 1 , 1]   1,  0  [ 5 , 5 , 5 , 5 , 5 , 3]  (0t415e29)
    (0t15)  [ 1 , 2 , 1 , 1 , 1 , 0]   1,  0  [ 5 , 6 , 5 , 5 , 5 , 2]  (0t418529)
    (0t18)  [ 2 , 1 , 2 , 1 , 0 , 0]   1,  1  [ 6 , 5 , 6 , 5 , 4 , 2]  (0t4185e9)
    (0t1e)  [ 1 , 1 , 2 , 0 , 1 , 1]   1,  0  [ 5 , 5 , 6 , 4 , 5 , 3]  (0t418e29)
    (0t41)  [ 3 , 2 , 1 , 0 , 0 , 0]   1,  1  [ 7 , 6 , 5 , 4 , 4 , 2]  (0t4185e2)
    (0t42)  [ 2 , 1 , 0 , 0 , 2 , 1]   1,  1  [ 6 , 5 , 4 , 4 , 6 , 3]  (0t415e23)
    (0t45)  [ 2 , 1 , 1 , 1 , 1 , 0]   1,  0  [ 6 , 5 , 5 , 5 , 5 , 2]  (0t418523)
    (0t48)  [ 2 , 2 , 1 , 1 , 0 , 0]   1,  0  [ 6 , 6 , 5 , 5 , 4 , 2]  (0t4185e3)
    (0t4e)  [ 2 , 1 , 0 , 1 , 1 , 1]   1,  0  [ 6 , 5 , 4 , 5 , 5 , 3]  (0t418e23)
    (0t52)  [ 1 , 1 , 0 , 1 , 2 , 1]   1,  0  [ 5 , 5 , 4 , 5 , 6 , 3]  (0t415239)
    (0t53)  [ 1 , 0 , 1 , 2 , 2 , 0]   1,  1  [ 5 , 4 , 5 , 6 , 6 , 2]  (0t485239)
    (0t5e)  [ 2 , 0 , 0 , 1 , 2 , 1]   1,  1  [ 6 , 4 , 4 , 5 , 6 , 3]  (0t418239)
    (0t82)  [ 1 , 0 , 2 , 1 , 1 , 1]   1,  0  [ 5 , 4 , 6 , 5 , 5 , 3]  (0t415e39)
    (0t83)  [ 1 , 0 , 1 , 3 , 1 , 0]   1,  0  [ 5 , 4 , 5 , 7 , 5 , 2]  (0t485e39)
    (0t85)  [ 2 , 0 , 1 , 2 , 1 , 0]   1,  1  [ 6 , 4 , 5 , 6 , 5 , 2]  (0t418539)
    (0t8e)  [ 1 , 1 , 1 , 1 , 1 , 1]   1,  0  [ 5 , 5 , 5 , 5 , 5 , 3]  (0t418e39)
    (0te2)  [ 2 , 0 , 0 , 0 , 2 , 2]   2,  2  [ 6 , 4 , 4 , 4 , 6 , 4]  (0t41e239)

We can also ask for the output in LaTeX:

```sage
print(ft2.latex())
```

    \begin{tabular}{llrll}
      \multicolumn{2}{c}{\textbf{Tetrachords}} &  & \multicolumn{2}{c}{\textbf{Octachords}} \\
      \midrule
      $(0153)$ & $012120$ & $1$, $1$ & $456562$ & $(0t185236)$\\
      $(0182)$ & $102210$ & $1$, $1$ & $546652$ & $(0t185e39)$\\
      $(01e9)$ & $004002$ & $4$, $4$ & $448444$ & $(0t18e296)$\\
      $(0412)$ & $111120$ & $1$, $0$ & $555562$ & $(0t485e29)$\\
      $(0415)$ & $122010$ & $1$, $1$ & $566452$ & $(04185e29)$\\
      $(041e)$ & $112101$ & $1$, $0$ & $556543$ & $(0t185e29)$\\
      $(0452)$ & $021030$ & $1$, $1$ & $465472$ & $(0t415236)$\\
      $(0453)$ & $012111$ & $1$, $0$ & $456553$ & $(0t485236)$\\
      $(0482)$ & $021120$ & $1$, $0$ & $465562$ & $(0t415e36)$\\
      $(0483)$ & $020301$ & $1$, $1$ & $464743$ & $(0t485e36)$\\
      $(048e)$ & $030201$ & $1$, $1$ & $474643$ & $(0t418e36)$\\
      $(04e3)$ & $020202$ & $2$, $2$ & $464644$ & $(0t48e236)$\\
      $(0t12)$ & $111111$ & $1$, $0$ & $555553$ & $(0t415e29)$\\
      $(0t15)$ & $121110$ & $1$, $0$ & $565552$ & $(0t418529)$\\
      $(0t18)$ & $212100$ & $1$, $1$ & $656542$ & $(0t4185e9)$\\
      $(0t1e)$ & $112011$ & $1$, $0$ & $556453$ & $(0t418e29)$\\
      $(0t41)$ & $321000$ & $1$, $1$ & $765442$ & $(0t4185e2)$\\
      $(0t42)$ & $210021$ & $1$, $1$ & $654463$ & $(0t415e23)$\\
      $(0t45)$ & $211110$ & $1$, $0$ & $655552$ & $(0t418523)$\\
      $(0t48)$ & $221100$ & $1$, $0$ & $665542$ & $(0t4185e3)$\\
      $(0t4e)$ & $210111$ & $1$, $0$ & $654553$ & $(0t418e23)$\\
      $(0t52)$ & $110121$ & $1$, $0$ & $554563$ & $(0t415239)$\\
      $(0t53)$ & $101220$ & $1$, $1$ & $545662$ & $(0t485239)$\\
      $(0t5e)$ & $200121$ & $1$, $1$ & $644563$ & $(0t418239)$\\
      $(0t82)$ & $102111$ & $1$, $0$ & $546553$ & $(0t415e39)$\\
      $(0t83)$ & $101310$ & $1$, $0$ & $545752$ & $(0t485e39)$\\
      $(0t85)$ & $201210$ & $1$, $1$ & $645652$ & $(0t418539)$\\
      $(0t8e)$ & $111111$ & $1$, $0$ & $555553$ & $(0t418e39)$\\
      $(0te2)$ & $200022$ & $2$, $2$ & $644464$ & $(0t41e239)$
    \end{tabular}


<a id="orgaa87e6a"></a>

## Class Elements

We can also print all set classes while listing every element in each set class. (Long output!)

For instance, let's look at out set classes of size four and their elements in the traditional system:

```sage
ce = ClassElements(4)
print(ce)
```

The output is too long to be displayed here, but it gives a series of columns, with the prime form of the set class on top, and the elements in the class below it.

We can also break the result in smaller number of colums:

```sage
ce.print_col(ncol=6)
```

Again, the output is too long, but is it more suitable for printing.


<a id="org4fb2f0f"></a>

## Webern Matrices

We can also construct Webern matrices. To see the original matrix, we can do:

```sage
wm = WebernMatrix()
print(wm)
```

        |   0  11   3   4   8   7   9   5   6   1   2  10  |
    ----|--------------------------------------------------|----
     0  |   0  11   3   4   8   7   9   5   6   1   2  10  |  10
     1  |   1   0   4   5   9   8  10   6   7   2   3  11  |  11
     9  |   9   8   0   1   5   4   6   2   3  10  11   7  |   7
     8  |   8   7  11   0   4   3   5   1   2   9  10   6  |   6
     4  |   4   3   7   8   0  11   1   9  10   5   6   2  |   2
     5  |   5   4   8   9   1   0   2  10  11   6   7   3  |   3
     3  |   3   2   6   7  11  10   0   8   9   4   5   1  |   1
     7  |   7   6  10  11   3   2   4   0   1   8   9   5  |   5
     6  |   6   5   9  10   2   1   3  11   0   7   8   4  |   4
    11  |  11  10   2   3   7   6   8   4   5   0   1   9  |   9
    10  |  10   9   1   2   6   5   7   3   4  11   0   8  |   8
     2  |   2   1   5   6  10   9  11   7   8   3   4   0  |   0
    ----|--------------------------------------------------|----
        |   2   1   5   6  10   9  11   7   8   3   4   0  |

The numbers on the left are the powers of the rotation that take the first row into the corresponding row. The numbers on top are the powers of the rotation that when composed with the reflection take the first row into the corresponding *column*. The right numbers and bottom numbers are similar, but with the *retrogrades*. See FIXME.

We can also get the whole matrix (without the labels) with `.matrix`:

```sage
wm.matrix
```

    [(0, 11, 3, 4, 8, 7, 9, 5, 6, 1, 2, 10),
     [1, 0, 4, 5, 9, 8, 10, 6, 7, 2, 3, 11],
     [9, 8, 0, 1, 5, 4, 6, 2, 3, 10, 11, 7],
     [8, 7, 11, 0, 4, 3, 5, 1, 2, 9, 10, 6],
     [4, 3, 7, 8, 0, 11, 1, 9, 10, 5, 6, 2],
     [5, 4, 8, 9, 1, 0, 2, 10, 11, 6, 7, 3],
     [3, 2, 6, 7, 11, 10, 0, 8, 9, 4, 5, 1],
     [7, 6, 10, 11, 3, 2, 4, 0, 1, 8, 9, 5],
     [6, 5, 9, 10, 2, 1, 3, 11, 0, 7, 8, 4],
     [11, 10, 2, 3, 7, 6, 8, 4, 5, 0, 1, 9],
     [10, 9, 1, 2, 6, 5, 7, 3, 4, 11, 0, 8],
     [2, 1, 5, 6, 10, 9, 11, 7, 8, 3, 4, 0]]

If you want just the labels, we can get them with `.left`, `.top`, `.right`, `.bottom`:

```sage
wm.left, wm.top, wm.right, wm.bottom
```

    ([0, 1, 9, 8, 4, 5, 3, 7, 6, 11, 10, 2],
     [0, 11, 3, 4, 8, 7, 9, 5, 6, 1, 2, 10],
     [10, 11, 7, 6, 2, 3, 1, 5, 4, 9, 8, 0],
     [2, 1, 5, 6, 10, 9, 11, 7, 8, 3, 4, 0])

Or, we can extract rows and columns (indexing starting at 0, as usual in Python/Sage):

```sage
wm.row(3), wm.column(8)
```

    ([8, 7, 11, 0, 4, 3, 5, 1, 2, 9, 10, 6],
     [6, 7, 3, 2, 10, 11, 9, 1, 0, 5, 4, 8])

We can also print it with LaTeX:

```sage
print(wm.latex())
```

    \begin{tabular}{c|cccccccccccc|c}
         &  $\phi$ &  $\rho^{11} \phi$ &  $\rho^{3} \phi$ &  $\rho^{4} \phi$ &  $\rho^{8} \phi$ &  $\rho^{7} \phi$ &  $\rho^{9} \phi$ &  $\rho^{5} \phi$ &  $\rho^{6} \phi$ &  $\rho \phi$ &  $\rho^{2} \phi$ &  $\rho^{10} \phi$ &    \\
      \hline
      $1$ & $0$ & $11$ & $3$ & $4$ & $8$ & $7$ & $9$ & $5$ & $6$ & $1$ & $2$ & $10$ & $\rho^{10}$ \\
      $\rho$ & $1$ & $0$ & $4$ & $5$ & $9$ & $8$ & $10$ & $6$ & $7$ & $2$ & $3$ & $11$ & $\rho^{11}$ \\
      $\rho^{9}$ & $9$ & $8$ & $0$ & $1$ & $5$ & $4$ & $6$ & $2$ & $3$ & $10$ & $11$ & $7$ & $\rho^{7}$ \\
      $\rho^{8}$ & $8$ & $7$ & $11$ & $0$ & $4$ & $3$ & $5$ & $1$ & $2$ & $9$ & $10$ & $6$ & $\rho^{6}$ \\
      $\rho^{4}$ & $4$ & $3$ & $7$ & $8$ & $0$ & $11$ & $1$ & $9$ & $10$ & $5$ & $6$ & $2$ & $\rho^{2}$ \\
      $\rho^{5}$ & $5$ & $4$ & $8$ & $9$ & $1$ & $0$ & $2$ & $10$ & $11$ & $6$ & $7$ & $3$ & $\rho^{3}$ \\
      $\rho^{3}$ & $3$ & $2$ & $6$ & $7$ & $11$ & $10$ & $0$ & $8$ & $9$ & $4$ & $5$ & $1$ & $\rho$ \\
      $\rho^{7}$ & $7$ & $6$ & $10$ & $11$ & $3$ & $2$ & $4$ & $0$ & $1$ & $8$ & $9$ & $5$ & $\rho^{5}$ \\
      $\rho^{6}$ & $6$ & $5$ & $9$ & $10$ & $2$ & $1$ & $3$ & $11$ & $0$ & $7$ & $8$ & $4$ & $\rho^{4}$ \\
      $\rho^{11}$ & $11$ & $10$ & $2$ & $3$ & $7$ & $6$ & $8$ & $4$ & $5$ & $0$ & $1$ & $9$ & $\rho^{9}$ \\
      $\rho^{10}$ & $10$ & $9$ & $1$ & $2$ & $6$ & $5$ & $7$ & $3$ & $4$ & $11$ & $0$ & $8$ & $\rho^{8}$ \\
      $\rho^{2}$ & $2$ & $1$ & $5$ & $6$ & $10$ & $9$ & $11$ & $7$ & $8$ & $3$ & $4$ & $0$ & $1$ \\
      \hline
         &  $\rho^{2} \phi$ &  $\rho \phi$ &  $\rho^{5} \phi$ &  $\rho^{6} \phi$ &  $\rho^{10} \phi$ &  $\rho^{9} \phi$ &  $\rho^{11} \phi$ &  $\rho^{7} \phi$ &  $\rho^{8} \phi$ &  $\rho^{3} \phi$ &  $\rho^{4} \phi$ &  $\phi$ &
    \end{tabular}

We can also create new Webern matrices using other systems, but first we need a first row. Let's use `MS2`. We can get possible first rows with:

```sage
first_webern_rows = MS2.find_all_first_webern_rows()
first_webern_rows
```

    [[(0, 1, 3), (2, 7, 8), (4, 5, 6), (9, 10, 11), [3, 3]],
     [(0, 1, 3), (2, 7, 8), (4, 5, 9), (6, 10, 11), [3, 3]],
     [(0, 1, 3), (2, 7, 8), (5, 6, 10), (4, 9, 11), [3, 3]],
     [(0, 1, 3), (5, 6, 10), (4, 9, 11), (2, 7, 8), [1, 1]],
     [(0, 2, 8), (1, 3, 7), (4, 5, 6), (9, 10, 11), [3, 3]],
     [(0, 2, 8), (1, 3, 7), (4, 5, 9), (6, 10, 11), [3, 3]],
     [(0, 2, 8), (1, 3, 7), (5, 6, 10), (4, 9, 11), [3, 3]],
     [(0, 2, 8), (4, 5, 9), (1, 3, 7), (6, 10, 11), [1, 1]],
     [(1, 6, 11), (2, 4, 7), (3, 5, 10), (0, 8, 9), [3, 3]],
     [(1, 6, 11), (2, 4, 7), (8, 9, 10), (0, 3, 5), [3, 3]],
     [(1, 6, 11), (3, 5, 10), (2, 4, 7), (0, 8, 9), [1, 1]],
     [(3, 5, 10), (0, 8, 9), (1, 7, 11), (2, 4, 6), [3, 3]],
     [(3, 8, 10), (0, 5, 9), (1, 6, 11), (2, 4, 7), [3, 3]],
     [(3, 8, 10), (0, 5, 9), (1, 7, 11), (2, 4, 6), [3, 3]],
     [(3, 8, 10), (0, 5, 9), (4, 7, 11), (1, 2, 6), [3, 3]],
     [(3, 8, 10), (1, 7, 11), (0, 5, 9), (2, 4, 6), [1, 1]],
     [(4, 5, 6), (0, 1, 2), (3, 7, 8), (9, 10, 11), [1, 1]],
     [(4, 5, 6), (9, 10, 11), (0, 1, 2), (3, 7, 8), [3, 3]],
     [(4, 5, 9), (6, 10, 11), (0, 1, 2), (3, 7, 8), [3, 3]],
     [(4, 7, 11), (1, 2, 6), (3, 5, 10), (0, 8, 9), [3, 3]],
     [(5, 6, 10), (4, 9, 11), (0, 1, 2), (3, 7, 8), [3, 3]],
     [(8, 9, 10), (0, 3, 5), (1, 7, 11), (2, 4, 6), [3, 3]],
     [(8, 9, 10), (0, 3, 5), (4, 7, 11), (1, 2, 6), [3, 3]],
     [(8, 9, 10), (4, 7, 11), (0, 3, 5), (1, 2, 6), [1, 1]],
     [(0, 4, 7), (1, 8, 10), (3, 5, 11), (2, 6, 9), [1, 1]],
     [(0, 4, 7), (2, 6, 9), (3, 5, 11), (1, 8, 10), [1, 1]],
     [(1, 8, 11), (2, 3, 5), (0, 6, 9), (4, 7, 10), [1, 1]],
     [(1, 8, 11), (4, 7, 10), (2, 3, 5), (0, 6, 9), [1, 1]],
     [(2, 3, 6), (0, 7, 9), (4, 8, 10), (1, 5, 11), [1, 1]],
     [(2, 3, 6), (1, 5, 11), (0, 7, 9), (4, 8, 10), [1, 1]],
     [(2, 9, 11), (3, 6, 7), (4, 5, 8), (0, 1, 10), [1, 1]],
     [(2, 9, 11), (4, 5, 8), (3, 6, 7), (0, 1, 10), [1, 1]],
     [(3, 7, 9), (0, 6, 10), (1, 4, 5), (2, 8, 11), [1, 1]],
     [(3, 7, 9), (2, 8, 11), (1, 4, 5), (0, 6, 10), [1, 1]],
     [(3, 9, 11), (2, 5, 8), (6, 7, 10), (0, 1, 4), [1, 1]],
     [(3, 9, 11), (6, 7, 10), (0, 1, 4), (2, 5, 8), [1, 1]],
     [(2, 8, 9), (0, 4, 5), (3, 7, 11), (1, 6, 10), [1, 1]],
     [(2, 8, 9), (3, 7, 11), (1, 6, 10), (0, 4, 5), [1, 1]],
     [(3, 5, 6), (0, 2, 9), (1, 10, 11), (4, 7, 8), [1, 1]],
     [(3, 5, 6), (1, 10, 11), (4, 7, 8), (0, 2, 9), [1, 1]],
     [(3, 6, 10), (0, 1, 5), (2, 4, 8), (7, 9, 11), [1, 1]],
     [(3, 6, 10), (7, 9, 11), (2, 4, 8), (0, 1, 5), [1, 1]],
     [(4, 7, 9), (0, 2, 6), (1, 3, 5), (8, 10, 11), [1, 1]],
     [(4, 7, 9), (8, 10, 11), (1, 3, 5), (0, 2, 6), [1, 1]],
     [(7, 8, 10), (0, 4, 9), (2, 5, 6), (1, 3, 11), [1, 1]],
     [(7, 8, 10), (1, 3, 11), (2, 5, 6), (0, 4, 9), [1, 1]],
     [(8, 9, 11), (2, 4, 5), (0, 1, 6), (3, 7, 10), [1, 1]],
     [(8, 9, 11), (3, 7, 10), (0, 1, 6), (2, 4, 5), [1, 1]],
     [(0, 1, 8), (2, 3, 7), (4, 5, 11), (6, 9, 10), [3, 3]],
     [(0, 3, 9), (2, 6, 7), (1, 4, 11), (5, 8, 10), [1, 1]],
     [(0, 3, 9), (5, 8, 10), (1, 2, 11), (4, 6, 7), [3, 3]],
     [(0, 3, 9), (5, 8, 10), (1, 4, 11), (2, 6, 7), [3, 3]],
     [(0, 3, 9), (5, 8, 10), (2, 6, 11), (1, 4, 7), [3, 3]],
     [(1, 2, 8), (0, 3, 7), (4, 5, 11), (6, 9, 10), [3, 3]],
     [(1, 2, 8), (0, 3, 7), (4, 6, 10), (5, 9, 11), [3, 3]],
     [(1, 2, 8), (0, 3, 7), (6, 9, 11), (4, 5, 10), [3, 3]],
     [(1, 2, 8), (4, 5, 11), (6, 9, 10), (0, 3, 7), [1, 1]],
     [(1, 2, 11), (3, 5, 8), (4, 6, 7), (0, 9, 10), [1, 1]],
     [(1, 2, 11), (4, 6, 7), (0, 9, 10), (3, 5, 8), [3, 3]],
     [(1, 4, 11), (2, 6, 7), (0, 9, 10), (3, 5, 8), [3, 3]],
     [(2, 6, 11), (1, 4, 7), (0, 9, 10), (3, 5, 8), [3, 3]],
     [(3, 5, 9), (0, 8, 10), (1, 2, 11), (4, 6, 7), [3, 3]],
     [(3, 5, 9), (0, 8, 10), (1, 4, 11), (2, 6, 7), [3, 3]],
     [(3, 5, 9), (0, 8, 10), (2, 6, 11), (1, 4, 7), [3, 3]],
     [(3, 5, 9), (2, 6, 11), (1, 4, 7), (0, 8, 10), [1, 1]],
     [(4, 5, 11), (6, 9, 10), (2, 3, 8), (0, 1, 7), [3, 3]],
     [(4, 6, 10), (0, 1, 7), (2, 3, 8), (5, 9, 11), [1, 1]],
     [(4, 6, 10), (5, 9, 11), (0, 1, 8), (2, 3, 7), [3, 3]],
     [(4, 6, 10), (5, 9, 11), (2, 3, 8), (0, 1, 7), [3, 3]],
     [(6, 9, 11), (2, 3, 7), (0, 1, 8), (4, 5, 10), [1, 1]],
     [(6, 9, 11), (4, 5, 10), (0, 1, 8), (2, 3, 7), [3, 3]],
     [(6, 9, 11), (4, 5, 10), (2, 3, 8), (0, 1, 7), [3, 3]],
     [(0, 2, 7), (1, 3, 8), (5, 10, 11), (4, 6, 9), [3, 3]],
     [(0, 2, 7), (4, 6, 9), (5, 10, 11), (1, 3, 8), [1, 1]],
     [(0, 2, 7), (5, 10, 11), (3, 8, 9), (1, 4, 6), [2, 2]],
     [(0, 2, 7), (5, 10, 11), (4, 6, 9), (1, 3, 8), [2, 2]],
     [(0, 5, 10), (1, 4, 6), (3, 8, 9), (2, 7, 11), [1, 1]],
     [(0, 5, 10), (2, 7, 11), (3, 8, 9), (1, 4, 6), [2, 2]],
     [(0, 5, 10), (2, 7, 11), (4, 6, 9), (1, 3, 8), [2, 2]],
     [(0, 5, 10), (3, 8, 9), (2, 4, 11), (1, 6, 7), [3, 3]],
     [(0, 5, 10), (3, 8, 9), (2, 7, 11), (1, 4, 6), [3, 3]],
     [(0, 5, 10), (3, 8, 9), (6, 7, 11), (1, 2, 4), [3, 3]],
     [(2, 4, 11), (0, 3, 10), (1, 6, 7), (5, 8, 9), [2, 2]],
     [(2, 4, 11), (0, 3, 10), (1, 7, 8), (5, 6, 9), [2, 2]],
     [(2, 4, 11), (1, 6, 7), (5, 8, 9), (0, 3, 10), [3, 3]],
     [(2, 4, 11), (5, 8, 9), (1, 6, 7), (0, 3, 10), [1, 1]],
     [(3, 9, 10), (0, 5, 8), (2, 4, 11), (1, 6, 7), [3, 3]],
     [(3, 9, 10), (0, 5, 8), (2, 7, 11), (1, 4, 6), [3, 3]],
     [(3, 9, 10), (0, 5, 8), (6, 7, 11), (1, 2, 4), [3, 3]],
     [(3, 9, 10), (1, 2, 4), (5, 6, 11), (0, 7, 8), [2, 2]],
     [(3, 9, 10), (1, 2, 4), (6, 7, 11), (0, 5, 8), [2, 2]],
     [(3, 9, 10), (6, 7, 11), (1, 2, 4), (0, 5, 8), [1, 1]],
     [(4, 9, 10), (0, 7, 8), (5, 6, 11), (1, 2, 3), [1, 1]],
     [(4, 9, 10), (1, 2, 3), (5, 6, 11), (0, 7, 8), [2, 2]],
     [(4, 9, 10), (1, 2, 3), (6, 7, 11), (0, 5, 8), [2, 2]],
     [(4, 9, 10), (5, 6, 11), (0, 2, 3), (1, 7, 8), [3, 3]],
     [(4, 9, 10), (5, 6, 11), (0, 2, 7), (1, 3, 8), [3, 3]],
     [(4, 9, 10), (5, 6, 11), (0, 7, 8), (1, 2, 3), [3, 3]],
     [(4, 10, 11), (0, 2, 3), (1, 6, 7), (5, 8, 9), [2, 2]],
     [(4, 10, 11), (0, 2, 3), (1, 7, 8), (5, 6, 9), [2, 2]],
     [(4, 10, 11), (1, 7, 8), (0, 2, 3), (5, 6, 9), [1, 1]],
     [(4, 10, 11), (5, 6, 9), (0, 2, 3), (1, 7, 8), [3, 3]],
     [(4, 10, 11), (5, 6, 9), (0, 2, 7), (1, 3, 8), [3, 3]],
     [(4, 10, 11), (5, 6, 9), (0, 7, 8), (1, 2, 3), [3, 3]],
     [(5, 8, 9), (0, 3, 10), (2, 7, 11), (1, 4, 6), [3, 3]],
     [(5, 10, 11), (4, 6, 9), (0, 2, 3), (1, 7, 8), [3, 3]],
     [(5, 10, 11), (4, 6, 9), (0, 7, 8), (1, 2, 3), [3, 3]],
     [(6, 7, 11), (1, 2, 4), (5, 8, 9), (0, 3, 10), [3, 3]],
     [(0, 2, 4), (1, 5, 6), (7, 8, 9), (3, 10, 11), [1, 1]],
     [(0, 2, 4), (3, 10, 11), (7, 8, 9), (1, 5, 6), [2, 2]],
     [(0, 2, 4), (7, 8, 9), (1, 5, 6), (3, 10, 11), [1, 1]],
     [(7, 8, 11), (0, 5, 6), (1, 3, 10), (2, 4, 9), [2, 2]],
     [(7, 8, 11), (1, 3, 10), (0, 5, 6), (2, 4, 9), [1, 1]],
     [(7, 8, 11), (2, 4, 9), (0, 5, 6), (1, 3, 10), [1, 1]],
     [(7, 10, 11), (0, 2, 5), (1, 3, 6), (4, 8, 9), [2, 2]],
     [(7, 10, 11), (1, 3, 6), (0, 2, 5), (4, 8, 9), [1, 1]],
     [(7, 10, 11), (4, 8, 9), (0, 2, 5), (1, 3, 6), [1, 1]],
     [(0, 2, 10), (5, 7, 11), (1, 6, 8), (3, 4, 9), [2, 2]],
     [(0, 2, 10), (5, 7, 11), (6, 8, 9), (1, 3, 4), [2, 2]],
     [(0, 5, 7), (2, 10, 11), (1, 6, 8), (3, 4, 9), [2, 2]],
     [(0, 5, 7), (2, 10, 11), (6, 8, 9), (1, 3, 4), [2, 2]],
     [(0, 10, 11), (2, 3, 4), (5, 6, 7), (1, 8, 9), [2, 2]],
     [(1, 3, 9), (2, 4, 10), (6, 7, 8), (0, 5, 11), [2, 2]],
     [(1, 6, 9), (5, 7, 8), (0, 10, 11), (2, 3, 4), [2, 2]],
     [(1, 6, 9), (5, 7, 8), (3, 4, 10), (0, 2, 11), [2, 2]],
     [(3, 4, 10), (0, 2, 11), (5, 6, 7), (1, 8, 9), [2, 2]],
     [(5, 6, 8), (0, 7, 11), (1, 3, 9), (2, 4, 10), [2, 2]],
     [(5, 6, 8), (0, 7, 11), (1, 4, 9), (2, 3, 10), [2, 2]],
     [(6, 7, 8), (0, 5, 11), (1, 4, 9), (2, 3, 10), [2, 2]],
     [(0, 3, 11), (4, 6, 8), (1, 2, 9), (5, 7, 10), [6, 6]],
     [(0, 4, 11), (3, 6, 8), (1, 2, 9), (5, 7, 10), [6, 6]],
     [(0, 4, 11), (3, 6, 8), (1, 5, 7), (2, 9, 10), [6, 6]],
     [(0, 4, 11), (3, 6, 8), (2, 5, 7), (1, 9, 10), [6, 6]],
     [(0, 4, 11), (3, 6, 8), (2, 5, 10), (1, 7, 9), [6, 6]],
     [(0, 4, 11), (3, 6, 8), (2, 7, 10), (1, 5, 9), [6, 6]],
     [(0, 4, 11), (3, 6, 8), (5, 7, 9), (1, 2, 10), [6, 6]],
     [(0, 6, 8), (3, 4, 11), (1, 2, 9), (5, 7, 10), [6, 6]],
     [(0, 6, 8), (3, 4, 11), (1, 5, 7), (2, 9, 10), [6, 6]],
     [(0, 6, 8), (3, 4, 11), (2, 5, 7), (1, 9, 10), [6, 6]],
     [(0, 6, 8), (3, 4, 11), (2, 5, 10), (1, 7, 9), [6, 6]],
     [(0, 6, 8), (3, 4, 11), (2, 7, 10), (1, 5, 9), [6, 6]],
     [(0, 6, 8), (3, 4, 11), (5, 7, 9), (1, 2, 10), [6, 6]],
     [(0, 6, 11), (3, 4, 8), (1, 2, 9), (5, 7, 10), [6, 6]],
     [(0, 6, 11), (3, 4, 8), (2, 5, 7), (1, 9, 10), [6, 6]],
     [(0, 6, 11), (3, 4, 8), (2, 5, 10), (1, 7, 9), [6, 6]],
     [(0, 6, 11), (3, 4, 8), (2, 7, 10), (1, 5, 9), [6, 6]],
     [(0, 6, 11), (3, 4, 8), (5, 7, 9), (1, 2, 10), [6, 6]],
     [(0, 8, 11), (3, 4, 6), (1, 2, 9), (5, 7, 10), [6, 6]],
     [(0, 8, 11), (3, 4, 6), (5, 7, 9), (1, 2, 10), [6, 6]],
     [(1, 5, 7), (2, 9, 10), (0, 3, 11), (4, 6, 8), [6, 6]],
     [(1, 5, 7), (2, 9, 10), (0, 6, 11), (3, 4, 8), [6, 6]],
     [(1, 5, 7), (2, 9, 10), (0, 8, 11), (3, 4, 6), [6, 6]],
     [(1, 5, 7), (2, 9, 10), (6, 8, 11), (0, 3, 4), [6, 6]],
     [(2, 5, 7), (1, 9, 10), (0, 3, 11), (4, 6, 8), [6, 6]],
     [(2, 5, 7), (1, 9, 10), (0, 8, 11), (3, 4, 6), [6, 6]],
     [(2, 5, 7), (1, 9, 10), (6, 8, 11), (0, 3, 4), [6, 6]],
     [(2, 5, 10), (1, 7, 9), (0, 3, 11), (4, 6, 8), [6, 6]],
     [(2, 5, 10), (1, 7, 9), (0, 8, 11), (3, 4, 6), [6, 6]],
     [(2, 5, 10), (1, 7, 9), (6, 8, 11), (0, 3, 4), [6, 6]],
     [(2, 7, 10), (1, 5, 9), (0, 3, 11), (4, 6, 8), [6, 6]],
     [(2, 7, 10), (1, 5, 9), (0, 8, 11), (3, 4, 6), [6, 6]],
     [(2, 7, 10), (1, 5, 9), (6, 8, 11), (0, 3, 4), [6, 6]],
     [(5, 7, 9), (1, 2, 10), (0, 3, 11), (4, 6, 8), [6, 6]],
     [(5, 7, 9), (1, 2, 10), (6, 8, 11), (0, 3, 4), [6, 6]],
     [(6, 8, 11), (0, 3, 4), (1, 2, 9), (5, 7, 10), [6, 6]],
     [(0, 4, 8), (1, 2, 5), (7, 9, 10), (3, 6, 11), [1, 1]],
     [(0, 4, 8), (1, 7, 10), (3, 6, 11), (2, 5, 9), [1, 1]],
     [(0, 4, 8), (2, 5, 9), (1, 7, 10), (3, 6, 11), [1, 1]],
     [(0, 4, 8), (3, 6, 11), (1, 7, 10), (2, 5, 9), [6, 6]],
     [(0, 4, 8), (3, 6, 11), (7, 9, 10), (1, 2, 5), [6, 6]],
     [(0, 4, 8), (7, 9, 10), (1, 2, 5), (3, 6, 11), [1, 1]],
     [(1, 7, 10), (0, 4, 6), (2, 5, 9), (3, 8, 11), [1, 1]],
     [(1, 7, 10), (2, 5, 9), (0, 4, 6), (3, 8, 11), [6, 6]],
     [(1, 7, 10), (3, 8, 11), (0, 4, 6), (2, 5, 9), [1, 1]],
     [(2, 7, 9), (0, 3, 6), (4, 8, 11), (1, 5, 10), [1, 1]],
     [(2, 7, 9), (0, 4, 6), (1, 5, 10), (3, 8, 11), [1, 1]],
     [(2, 7, 9), (0, 4, 8), (1, 5, 10), (3, 6, 11), [1, 1]],
     [(2, 7, 9), (1, 5, 10), (0, 4, 6), (3, 8, 11), [6, 6]],
     [(2, 7, 9), (1, 5, 10), (0, 4, 8), (3, 6, 11), [6, 6]],
     [(2, 7, 9), (1, 5, 10), (4, 8, 11), (0, 3, 6), [6, 6]],
     [(2, 7, 9), (3, 6, 11), (0, 4, 8), (1, 5, 10), [1, 1]],
     [(2, 7, 9), (3, 8, 11), (1, 5, 10), (0, 4, 6), [1, 1]],
     [(2, 7, 9), (4, 8, 11), (1, 5, 10), (0, 3, 6), [1, 1]],
     [(4, 8, 11), (0, 3, 6), (1, 7, 10), (2, 5, 9), [6, 6]],
     [(4, 8, 11), (0, 3, 6), (7, 9, 10), (1, 2, 5), [6, 6]],
     [(4, 8, 11), (1, 2, 5), (7, 9, 10), (0, 3, 6), [1, 1]],
     [(4, 8, 11), (1, 7, 10), (0, 3, 6), (2, 5, 9), [1, 1]],
     [(4, 8, 11), (2, 5, 9), (1, 7, 10), (0, 3, 6), [1, 1]],
     [(4, 8, 11), (7, 9, 10), (0, 3, 6), (1, 2, 5), [1, 1]],
     [(7, 9, 10), (0, 4, 6), (1, 2, 5), (3, 8, 11), [1, 1]],
     [(7, 9, 10), (1, 2, 5), (0, 4, 6), (3, 8, 11), [6, 6]],
     [(7, 9, 10), (3, 8, 11), (1, 2, 5), (0, 4, 6), [1, 1]],
     [(1, 2, 7), (0, 3, 8), (5, 9, 10), (4, 6, 11), [3, 3]],
     [(1, 2, 7), (4, 6, 11), (5, 9, 10), (0, 3, 8), [3, 3]],
     [(1, 2, 7), (5, 9, 10), (0, 3, 8), (4, 6, 11), [6, 6]],
     [(0, 7, 10), (1, 4, 8), (2, 5, 11), (3, 6, 9), [1, 1]],
     [(0, 7, 10), (2, 5, 11), (3, 6, 9), (1, 4, 8), [2, 2]],
     [(0, 7, 10), (3, 6, 9), (2, 5, 11), (1, 4, 8), [1, 1]],
     [(2, 3, 9), (0, 6, 7), (1, 4, 10), (5, 8, 11), [1, 1]],
     [(2, 3, 9), (1, 4, 10), (0, 6, 7), (5, 8, 11), [2, 2]],
     [(2, 3, 9), (5, 8, 11), (1, 4, 10), (0, 6, 7), [1, 1]],
     [(6, 7, 9), (0, 4, 10), (2, 3, 11), (1, 5, 8), [1, 1]],
     [(6, 7, 9), (1, 5, 8), (0, 4, 10), (2, 3, 11), [2, 2]],
     [(6, 7, 9), (2, 3, 11), (0, 4, 10), (1, 5, 8), [1, 1]]]

The result is divided in four trichords, all in the same class, and the number of symmetries of the corresponding class.

For instance, let's choose the thirteenth row:

```sage
first_webern_rows[12]
```

    [(3, 8, 10), (0, 5, 9), (1, 6, 11), (2, 4, 7), [3, 3]]

We can see that the trichords are indeed related by symmetry, by asking what are the symmetries that map one to another. For instance:

```sage
wrow = first_webern_rows[12]
t1, t2, t3, t4, _ = wrow
MS2.symmetry_maps(t1, t2), MS2.symmetry_maps(t2, t4)
```

    ([(1, 1)], [(2, 0)])

This means that `rho_2 * phi2` maps `[3, 8, 10]` to `[0, 5, 9]` and `rho2^2` maps `[0, 5, 9]` to `[2, 4, 7]`.

We can see if this first row has hexachords that are also related by symmetries:

```sage
MS2.symmetry_maps(t1 + t2, t3 + t4)
```

    [(11, 1), (7, 1), (3, 1), (2, 0), (10, 0), (6, 0)]

Indeed, the hexachords are related by three rotations (transpositions) and three reflection (inversions)!

So, now we can use this first row to create a Webern matrix in this system. We can scramble the order of the trichors, and the other of the pitch classes inside each trichord:

```sage
first_row = [0, 9, 5, 7, 4, 2, 1, 6, 11, 10, 3, 8]
wm2 = WebernMatrix(first_row, MS=MS2)
print(wm2)
```

        |   0   9   5  11   2   7   3  10   6   1   8   4  |
    ----|--------------------------------------------------|----
     0  |   0   9   5   7   4   2   1   6  11  10   3   8  |   4
     3  |   1   0   3   4   5   6  11  10   9   8   7   2  |   7
     7  |   2   8   0  11   9   4   6   5  10   3   1   7  |  11
     1  |  10   6  11   0   1   3   8   7   2   4   9   5  |   5
    10  |   6   2   1   9   0   5  10   3   8   7  11   4  |   2
     5  |   5   4   6   8   2   0   3   1   7  11  10   9  |   9
     9  |   9  11   4   3   7   8   0   2   1   6   5  10  |   1
     2  |   4   7   2  10   8   9   5   0   3   1   6  11  |   6
     6  |  11   1   7   5   3  10   9   8   0   2   4   6  |  10
    11  |   7   3   8   6  10  11   4   9   5   0   2   1  |   3
     4  |   8  10   9   1  11   7   2   4   6   5   0   3  |   8
     8  |   3   5  10   2   6   1   7  11   4   9   8   0  |   0
    ----|--------------------------------------------------|----
        |   8   5   1   7  10   3  11   6   2   9   4   0  |

To print in LaTeX, we can also give the names to the maps. In this case, we can use `rho2` and `phi2`:

```sage
print(wm2.latex('\\rho_2', '\\phi_2'))
```

    \begin{tabular}{c|cccccccccccc|c}
         &  $\phi_2$ &  $\rho_2^{9} \phi_2$ &  $\rho_2^{5} \phi_2$ &  $\rho_2^{11} \phi_2$ &  $\rho_2^{2} \phi_2$ &  $\rho_2^{7} \phi_2$ &  $\rho_2^{3} \phi_2$ &  $\rho_2^{10} \phi_2$ &  $\rho_2^{6} \phi_2$ &  $\rho_2 \phi_2$ &  $\rho_2^{8} \phi_2$ &  $\rho_2^{4} \phi_2$ &    \\
      \hline
      $1$ & $0$ & $9$ & $5$ & $7$ & $4$ & $2$ & $1$ & $6$ & $11$ & $10$ & $3$ & $8$ & $\rho_2^{4}$ \\
      $\rho_2^{3}$ & $1$ & $0$ & $3$ & $4$ & $5$ & $6$ & $11$ & $10$ & $9$ & $8$ & $7$ & $2$ & $\rho_2^{7}$ \\
      $\rho_2^{7}$ & $2$ & $8$ & $0$ & $11$ & $9$ & $4$ & $6$ & $5$ & $10$ & $3$ & $1$ & $7$ & $\rho_2^{11}$ \\
      $\rho_2$ & $10$ & $6$ & $11$ & $0$ & $1$ & $3$ & $8$ & $7$ & $2$ & $4$ & $9$ & $5$ & $\rho_2^{5}$ \\
      $\rho_2^{10}$ & $6$ & $2$ & $1$ & $9$ & $0$ & $5$ & $10$ & $3$ & $8$ & $7$ & $11$ & $4$ & $\rho_2^{2}$ \\
      $\rho_2^{5}$ & $5$ & $4$ & $6$ & $8$ & $2$ & $0$ & $3$ & $1$ & $7$ & $11$ & $10$ & $9$ & $\rho_2^{9}$ \\
      $\rho_2^{9}$ & $9$ & $11$ & $4$ & $3$ & $7$ & $8$ & $0$ & $2$ & $1$ & $6$ & $5$ & $10$ & $\rho_2$ \\
      $\rho_2^{2}$ & $4$ & $7$ & $2$ & $10$ & $8$ & $9$ & $5$ & $0$ & $3$ & $1$ & $6$ & $11$ & $\rho_2^{6}$ \\
      $\rho_2^{6}$ & $11$ & $1$ & $7$ & $5$ & $3$ & $10$ & $9$ & $8$ & $0$ & $2$ & $4$ & $6$ & $\rho_2^{10}$ \\
      $\rho_2^{11}$ & $7$ & $3$ & $8$ & $6$ & $10$ & $11$ & $4$ & $9$ & $5$ & $0$ & $2$ & $1$ & $\rho_2^{3}$ \\
      $\rho_2^{4}$ & $8$ & $10$ & $9$ & $1$ & $11$ & $7$ & $2$ & $4$ & $6$ & $5$ & $0$ & $3$ & $\rho_2^{8}$ \\
      $\rho_2^{8}$ & $3$ & $5$ & $10$ & $2$ & $6$ & $1$ & $7$ & $11$ & $4$ & $9$ & $8$ & $0$ & $1$ \\
      \hline
         &  $\rho_2^{8} \phi_2$ &  $\rho_2^{5} \phi_2$ &  $\rho_2 \phi_2$ &  $\rho_2^{7} \phi_2$ &  $\rho_2^{10} \phi_2$ &  $\rho_2^{3} \phi_2$ &  $\rho_2^{11} \phi_2$ &  $\rho_2^{6} \phi_2$ &  $\rho_2^{2} \phi_2$ &  $\rho_2^{9} \phi_2$ &  $\rho_2^{4} \phi_2$ &  $\phi_2$ &
    \end{tabular}

Note that we need the double `\` for the LaTeX names, as in `'\\rho_2`.
