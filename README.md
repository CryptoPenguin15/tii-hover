# HOV/HOVER key recovery

This package implements the structured HOV/HOVER attack from Saarinen's
**ePrint 2026/1778**.

See also https://github.com/mjosaarinen/tii-solved/tree/main for verification.

## Results

| Challenge | m | r | n | ratio | method | wall time | system |
|-----------|---|---|------|--------|--------|-----------|--------|
| **TII-252** | 10 | 11 | 1008 | 89.1% | HOVER | 1 min 55 s | m1 |
| **TII-246** | 10 | 11 | 1009 | 89.1% | HOVER | 1 min 54 s | m1 |
| **TII-240** | 10 | 11 | 1010 | 89.1% | HOVER | 1 min 55 s | m1 |
| **TII-248** | 9 | 7 | 482 | 86.9% | HOVER | 11 s | m1 |
| **TII-83** | 8 | 5 | 253 | 84.2% | HOVER | 5.5 s | m1 |

## Run

Activate SageMath environment. Sage/Cython automatically compiles
`hov_native.pyx` on first use; a C compiler and pkg-config must be available
in that environment. It links Sage's M4RI library and contains no
architecture-specific intrinsics. M4GB and CUDA are not required.

```bash
./recover_key.sh <TII instance number>
```

By default the program also checks a matched random code has zero HOV kernel.
`--skip-random-control` skips that stage for diagnostic reruns.
Final exact key verification always runs.

Successful output is JSON plus `<output>.txt` in official two-line TII format.
A failed acceptance check returns exit code 2 without writing a new key.
An existing output is left in place on failure: check exit code and input hash.

## Resource profiles

| Input | Degree | Shortening | Variables kappa | Matrix rows x columns | Packed matrix alone |
|---|---:|---:|---:|---:|---:|
| TII-100 | 4 | 16 | 48 | 214,510 x 194,580 | 4.860 GiB |
| TII-252 | 3 | 36 | 74 | 66,452 x 64,824 | 0.502 GiB |

Sage, M4RI workspaces and completion matrices require additional memory.
Kernel elimination is in place, avoiding a second full matrix copy. This uses
dense bit-packed M4RI elimination, so its scaling differs from the original
implementation's iterative sparse kernel computation. It is intended for
reachable degree-3/4 challenge profiles, not full Classic McEliece parameters.

```bash
python hover.py data/pk_McEliece_100.txt --profile-only
```

`--degree auto` is the default: it selects the first of 3,4 meeting a necessary
matrix-dimension screen, without guaranteeing a clean kernel. `--trials`
defaults to 3 public shortenings, `--seed` to 0, and `--max-matrix-gib` to 6.
The memory screen measures packed matrix storage, not total process memory.

## Implemented recovery chain

1. Read the original TII matrix and defining field polynomial. Perform public
   permutation, systematic reduction and shortening using equation (3).
2. Construct the modified HOV map: first derivatives for p=3; values and second
   square-free derivatives for p=4. Compute its binary kernel and independently
   replay all basis vectors through the streaming operator.
3. Construct the contraction tensor and equation (11). Check the compatible
   algebra's dimension m, identity, multiplication closure and an irreducible
   degree-m field generator.
4. Split directly over the input's GF(2^m). Check one-dimensional eigenlines,
   one Frobenius cycle, direct-image ranks, separated derivative kernels and
   every cyclic-window rank bound.
5. When a direction has rank p^2, recover its rational-normal-scroll rulings
   from tangent quadrics. Their Pluecker coordinates give a rational-normal
   curve of degree p(p-1). Cross-ratio GRS recovery is accepted only after exact
   reconstruction of the GRS code. This is a new tangent/Pluecker completion
   implementation of the geometry in Proposition 2, rather than the author's
   quadratic-ideal/linear-syzygy code. The same construction supports p=4 when
   all its exact rank gates pass.
6. Otherwise, accumulate block-omitting windows across shortenings and enumerate
   cyclic phases. Require the full expected Frobenius-intersection dimension
   trace before GRS extraction. This fallback does not implement every
   partial-space/Schur-product completion described by Hemmert; it refuses
   profiles where full block sums cannot be reached.
7. Recover g from even-weight codeword-product derivatives, determine the
   projective pole using an odd-weight word when available, and complete missing
   locators using the public row transformation.
8. Require monic irreducible g, distinct support, correct parameters and exact
   equality of reconstructed and public binary parity-check row spaces.

## Sources

Markku-Juhani O. Saarinen, HOVER: Higher-Order Vanishing Endomorphism Recovery,
ePrint 2026/1778, especially sections 3-5 and equations (3), (11), (20):
https://eprint.iacr.org/2026/1778

Original public inputs and official format:
https://github.com/ElenaKirshanova/tii_decoding_challenge/tree/main/public_keyRec

## AI attribution

This implementation was generated entirely by ChatGPT (OpenAI),
with prompts and guidance from Anders Nilson.
