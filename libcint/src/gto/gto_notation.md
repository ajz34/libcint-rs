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
  - For example of Jensen basis set pc-0, carbon atom, the $s$ shell. There is a 
    ```
    ```

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