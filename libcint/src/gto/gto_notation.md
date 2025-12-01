# Notation definitions

## Common functions

- Cartesian GTO

  $$
  | l_x, l_y, l_z, \alpha, \bm r \rangle = | \bm l, \alpha, \bm r \rangle = x^{l_x} y^{l_y} z^{l_z} \exp (-\alpha |\bm r|^2)
  $$

  We may use $r^2$ to denote $|\bm r|^2$ for simplicity. Please note that the coordinates above are centered by the atom, not the grid points centered by the origin $(0, 0, 0)$.

  We will denote the $x^{l_x} y^{l_y} z^{l_z}$ part as the polynomial part, and $\exp (-\alpha |\bm r|^2)$ part as the exponent part.

  This cartesian GTO function is not normalized. We expect the normalization factor to be multiplied by contraction coefficients $C_{k p}$, which is specified in struct fields `bas` and `env` of [`CInt`].

- Primitive GTO (`prim`)

  For certain angular momentum $\bm l$, the primitive GTO is denoted as

  $$
  | \bm l, \alpha_p, \bm r_g \rangle
  $$

  Note that subscript $p$ denotes primitive index, not component of coordinates ($x, y, z$). Also, in codes, we usually evaluate the exponent part and polynomial part separately. We denote the exponent part as $| \bm 0, \alpha_p, \bm r_g \rangle = \exp (- \alpha |\bm r|^2)$, but please note that this is not the same to the full primitive GTO $| \bm l, \alpha_p, \bm r_g \rangle$.

- Contracted GTO (`ctr`, `gto`)

  We usually denote atomic orbital index as $\mu, \nu$, and the atomic orbital as

  $$
  \phi_{\mu g} = \phi_\mu (\bm r_g) = \sum_{p} C_{k p} | \bm l, \alpha_p, \bm r_g \rangle
  $$

  where $k$ is the contracted function index, $C_{k p}$ is the contraction coefficient (defined by basis set, normalization factor included additionally).

- Atomic orbital (`ao`)

  Atomic orbital by definition is the contracted GTO. In program, we evaluated variable `gto` shell-by-shell, and assemble them into variable `ao`. They are eventually the same thing, so we may mix the notations of `ao` and `gto` in function APIs.

## Notes on contraction coefficients $C_{k p}$

  - $k$ is the contraction index. The size of $k$ is usually 1 (e.g., Pople, Ahlrichs), but may also be larger than 1 (e.g., Dunning, Jensen).
  - For example of Jensen basis set pc-0. An $s$ shell of carbon in Gaussian94 format is
    ```log
    S    3   1.00
          0.446510D+03           0.197880D-01
          0.674620D+02           0.145340D+00
          0.150670D+02           0.185290D+00
    S    3   1.00
          0.674620D+02          -0.999100D-02
          0.150670D+02           0.285190D+00
          0.393880D+01           0.566870D+00
    ```
    This can actually be considered as only one shell, with two contracted functions (i.e. $n_\mathrm{ctr} = 2$ at index $k$), each contracted function having four primitive functions (i.e. $n_\mathrm{prim} = 4$ at index $p$). The un-normalized contraction coefficients $\tilde C_{k p}$ and corresponding exponents $\alpha_p$ are:

    |                                     | $p = 0$         | $p = 1$         | $p = 2$         | $p = 3$         |
    |-------------------------------------|----------------:|----------------:|----------------:|----------------:|
    | coefficient $\tilde C_{k=0, p}$     | 0.019788        |  0.145340       | 0.185290        | 0.000000        |
    | coefficient $\tilde C_{k=1, p}$     | 0.000000        | -0.009991       | 0.285190        | 0.566870        |
    | exponent    $\alpha_p$              | 446.510         | 67.4620         | 15.0670         | 3.93880         |

    Please note that for libcint, the coefficient $\tilde C_{kp}$ does not include normalization factor, and you need to merge the normalization factor into the coefficient $C_{k p}$ to struct field `env` of [`CInt`] manually.

## Notes on relation of atomic orbital index $\mu$ and contracted GTO index $k$

  - We are considering a single shell here.
  - $\mu$ represents the cartesian atomic orbital basis. The size of $\mu$ is the multiplication of $n_\mathrm{cart} = (l + 1) (l + 2) / 2$, and $n_\mathrm{ctr}$ (the number of contracted functions in the shell). We can say that index $\mu$ can be mapped to a tuple of angular momentum $\bm l$ and contracted function index $k$:

    $$
    \mu := (\bm l, k) = (l_x, l_y, l_z, k)
    $$

  - $l_x$, $l_y$, $l_z$ are the cartesian exponents of the basis function, satisfying $l_x + l_y + l_z = l$. Note that for a shell of basis, the angular momentum $l$ is fixed for all contracted functions.
  - The ordering of $(l_x, l_y, l_z)$ follows the convention in PySCF, i.e., the most rapidly changing index is $l_z$, then $l_y$, and finally $l_x$. Counts of all orderings is $n_\mathrm{cart}$.
  - $k$ is the index of contracted function. The size of $k$ is $n_\mathrm{ctr}$.

  For example of angular momentum $l=2$ (d shell) with $n_\mathrm{ctr} = 2$, the ordering of $\mu$ (without offset that given by [`make_loc`](CInt::make_loc_with_type)) is:

  | $\mu$ | $(l_x, l_y, l_z)$ | $k$ || $\mu$ | $(l_x, l_y, l_z)$ | $k$ |
  |--|--|--|--|--|--|--|
  |  0 | (2, 0, 0) | 0 ||  6 | (2, 0, 0) | 1 |
  |  1 | (1, 1, 0) | 0 ||  7 | (1, 1, 0) | 1 |
  |  2 | (1, 0, 1) | 0 ||  8 | (1, 0, 1) | 1 |
  |  3 | (0, 2, 0) | 0 ||  9 | (0, 2, 0) | 1 |
  |  4 | (0, 1, 1) | 0 || 10 | (0, 1, 1) | 1 |
  |  5 | (0, 0, 2) | 0 || 11 | (0, 0, 2) | 1 |

## Common indices

| index | size | notes |
|-------|------|-------|
| $p$ | `nprim` | primitive GTO |
| $k$ | `nctr` | GTO contraction |
| $g$ | `ngrids` | grid point |
|     | [`BLKSIZE`] | grid block size for AO evaluation in iterations |
|     | [`BLKSIMDD`] | grid block size in SIMD lanes for AO evaluation in iterations |
| $l$ | `l` | angular momentum |
| $t$ | `3` | coordinate index $x, y, z$ |