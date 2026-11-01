#!/usr/bin/env python3
"""Portable TII Track-2 input, finite fields, key verification and preflight.

This is infrastructure for a recovery implementation, NOT a recovery solver.
No dependencies; also importable with Sage's Python on Apple Silicon.
"""
import argparse
import ast
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path


def binary_rref(rows):
    """Unique row-space representation, using Python integers as packed rows."""
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    for pivot in sorted(basis):
        for other in basis:
            if other > pivot and (basis[other] >> pivot) & 1:
                basis[other] ^= basis[pivot]
    return tuple(basis[p] for p in sorted(basis, reverse=True))


class BinaryField:
    def __init__(self, modulus):
        self.modulus = modulus
        self.m = modulus.bit_length() - 1
        self.q = 1 << self.m
        if self.m < 1 or not modulus & 1:
            raise ValueError('Invalid binary field modulus')

    def mul(self, x, y):
        out = 0
        while y:
            if y & 1:
                out ^= x
            y >>= 1
            x <<= 1
            if x & self.q:
                x ^= self.modulus
        return out

    def power(self, x, e):
        out = 1
        while e:
            if e & 1:
                out = self.mul(out, x)
            x = self.mul(x, x)
            e >>= 1
        return out

    def inverse(self, x):
        if not x:
            raise ZeroDivisionError('Zero field element')
        return self.power(x, self.q - 2)

    def evaluate(self, coeffs, x):
        out = 0
        for c in reversed(coeffs):
            out = self.mul(out, x) ^ c
        return out


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def remainder(a, b, field):
    a, b = trim(a), trim(b)
    if not b:
        raise ZeroDivisionError('Zero polynomial')
    inv = field.inverse(b[-1])
    while len(a) >= len(b):
        c = field.mul(a[-1], inv)
        offset = len(a) - len(b)
        for j, v in enumerate(b):
            a[offset+j] ^= field.mul(c, v)
        a = trim(a)
    return a


def poly_gcd(a, b, field):
    while b:
        a, b = b, remainder(a, b, field)
    if not a:
        return []
    inv = field.inverse(a[-1])
    return [field.mul(c, inv) for c in a]


def poly_square_mod(a, g, field):
    out = [0] * (2 * len(a) - 1)
    for j, c in enumerate(a):
        out[2*j] = field.mul(c, c)
    return remainder(out, g, field)


def irreducible(g, field):
    """Rabin's criterion over GF(2^m), no probabilistic assumptions."""
    g = trim(g)
    degree = len(g) - 1
    if degree < 1:
        return False
    x = remainder([0, 1], g, field)
    h = x
    for i in range(1, degree + 1):
        for _ in range(field.m):
            h = poly_square_mod(h, g, field)
        diff = list(h) + [0] * max(0, len(x)-len(h))
        for j, c in enumerate(x):
            diff[j] ^= c
        diff = trim(diff)
        if i <= degree // 2 and len(poly_gcd(g, diff, field)) != 1:
            return False
    return not diff


@dataclass
class PublicKey:
    rows: list
    n: int
    modulus_coefficients: list
    sha256: str

    @property
    def m(self):
        return len(self.modulus_coefficients) - 1

    @property
    def t(self):
        if len(self.rows) % self.m:
            raise ValueError('Row count is not divisible by extension degree')
        return len(self.rows) // self.m

    @property
    def k(self):
        return self.n - len(self.rows)

    @property
    def field(self):
        return BinaryField(sum(c << j for j, c in enumerate(self.modulus_coefficients)))

    def sage_objects(self):
        from sage.all import GF, PolynomialRing, matrix
        base = GF(2)
        ring = PolynomialRing(base, 'z')
        ff = GF(2**self.m, name='a', modulus=ring(self.modulus_coefficients))
        H = matrix(base, [[(r >> j) & 1 for j in range(self.n)] for r in self.rows])
        return H, ff


def read_public_key(path):
    raw = Path(path).read_bytes()
    lines = [s.strip() for s in raw.decode('ascii').splitlines() if s.strip()]
    if len(lines) < 2:
        raise ValueError('Expected matrix rows and final field polynomial')
    try:
        coeffs = ast.literal_eval(lines[-1])
    except (SyntaxError, ValueError) as e:
        raise ValueError('Invalid trailing field polynomial') from e
    if (not isinstance(coeffs, list) or len(coeffs) < 2
            or any(type(c) is not int or c not in (0, 1) for c in coeffs)
            or coeffs[-1] != 1):
        raise ValueError('Expected low-to-high monic binary modulus coefficients')
    rows, n = [], None
    for line in lines[:-1]:
        if not re.fullmatch(r'\[\s*[01](?:\s+[01])*\s*\]', line):
            raise ValueError('Invalid matrix row: expected bracketed, space-separated bits')
        bits = list(map(int, line[1:-1].split()))
        if n is not None and len(bits) != n:
            raise ValueError('Inconsistent matrix row lengths')
        n = len(bits)
        rows.append(sum(b << j for j, b in enumerate(bits)))
    result = PublicKey(rows, n, coeffs, hashlib.sha256(raw).hexdigest())
    if n > 2**result.m:
        raise ValueError('More support coordinates than field elements')
    if not irreducible(coeffs, BinaryField(3)):
        raise ValueError('Field modulus is reducible over GF(2)')
    if len(binary_rref(rows)) != len(rows):
        raise ValueError('Public matrix does not have full row rank')
    if result.k <= 0:
        raise ValueError('Invalid code dimension')
    _ = result.t
    return result


def parse_field_expression(text, field):
    """Parse TII's polynomial-in-a notation without eval or Sage preparsing."""
    text = text.replace(' ', '')
    if not text:
        raise ValueError('Empty field element')
    result = 0
    for term in text.split('+'):
        if term in ('0', '1'):
            value = int(term)
        elif term == 'a':
            value = field.power(2, 1)
        elif re.fullmatch(r'a\^[0-9]+', term):
            exponent = int(term[2:])
            if exponent > 100000:
                raise ValueError('Field exponent exceeds parser limit')
            value = field.power(2, exponent)
        else:
            raise ValueError('Unsupported field expression: ' + term)
        result ^= value
    return result


def read_secret_key(path, field):
    path = Path(path)
    if path.suffix == '.json':
        obj = json.loads(path.read_text())
        return obj['g'], obj['support']
    lines = [s.strip() for s in path.read_text().splitlines() if s.strip()]
    if len(lines) != 2:
        raise ValueError('TII secret key must have exactly two lines: g, support')
    out = []
    for s in lines:
        if not s.startswith('[') or not s.endswith(']'):
            raise ValueError('Expected bracketed coefficient/support list')
        out.append([parse_field_expression(x, field) for x in s[1:-1].split(',')])
    return tuple(out)


def parity_rows(support, g, field):
    rows = [0] * ((len(g)-1) * field.m)
    for j, a in enumerate(support):
        value = field.inverse(field.evaluate(g, a))
        for i in range(len(g)-1):
            for bit in range(field.m):
                if (value >> bit) & 1:
                    rows[i*field.m+bit] |= 1 << j
            value = field.mul(value, a)
    return rows


def verify_key(pk, g, support):
    field = pk.field
    if len(g) != pk.t + 1 or g[-1] != 1:
        raise ValueError('Goppa polynomial must be monic of degree t')
    if len(support) != pk.n or len(set(support)) != pk.n:
        raise ValueError('Support must contain n distinct elements in public-column order')
    if any(type(x) is not int or not 0 <= x < field.q for x in g + support):
        raise ValueError('Field coefficients must be polynomial-basis integers')
    if not irreducible(g, field):
        raise ValueError('Goppa polynomial is reducible')
    rebuilt = parity_rows(support, g, field)
    if binary_rref(rebuilt) != binary_rref(pk.rows):
        raise ValueError('Reconstructed code differs from the public code')
    return True


def generic_kernel_dimension(k, d, multiplicities):
    """Weis eq. (4): positive-part Hilbert-series heuristic, NOT a rank proof."""
    coeffs = [math.comb(k, i) if i <= k else 0 for i in range(d+1)]
    for mu in multiplicities:
        shift = d - mu + 1
        for i in range(d, shift-1, -1):
            coeffs[i] -= coeffs[i-shift]
    return coeffs[d] if all(c > 0 for c in coeffs) else 0


def preflight(pk, max_degree=8):
    """Gradient-route screening for uniform s=d-1, one value-only holdout.

    Enumerates shortening dimensions and the infinity-point multiplicity.
    Uses conservative degree-zero forcing sum(mu)>d*(n_short-2*t-1).
    It deliberately does not assert the sharper probabilistic Theorem-3 bound.
    """
    choices = []
    nonunit = len(pk.rows)
    for d in range(3, max_degree+1):
        found = None
        for k in range(d, pk.k+1):
            for infinity in range(d):
                multiplicities = [d-1] * (nonunit-1) + [1]
                if infinity:
                    multiplicities.append(infinity)
                kappa = generic_kernel_dimension(k, d, multiplicities)
                zeros = k*(d-1) + (nonunit-1)*(d-1) + 1 + infinity
                curve_degree = k+nonunit-2*pk.t-1
                if kappa and zeros > d*curve_degree:
                    N = math.comb(k, d)
                    # Full Hasse derivative matrix, before any reductions.
                    conditions = sum(sum(math.comb(k,u) for u in range(mu))
                                     for mu in multiplicities)
                    found = dict(degree=d, shortened_dimension=k,
                                 shortened_length=k+nonunit, infinity_multiplicity=infinity,
                                 predicted_kernel=kappa, monomials=N,
                                 full_hasse_conditions=conditions,
                                 dense_matrix_GiB=N*conditions/8/2**30,
                                 forced_zeros=zeros, polynomial_degree_bound=d*curve_degree)
                    break
            if found:
                break
        if found:
            choices.append(found)
    return choices


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('inspect', 'preflight', 'verify'):
        p = sub.add_parser(name)
        p.add_argument('public_key', type=Path)
        if name == 'verify':
            p.add_argument('secret_key', type=Path)
    args = parser.parse_args()
    pk = read_public_key(args.public_key)
    if args.command == 'inspect':
        out = dict(n=pk.n, k=pk.k, t=pk.t, m=pk.m, rank=len(pk.rows),
                   modulus=pk.modulus_coefficients, sha256=pk.sha256)
    elif args.command == 'preflight':
        out = dict(method='conservative gradient-route screening',
                   warning='Heuristic kernel dimensions; dense sizes describe the unreduced matrix, not a memory lower bound.',
                   candidates=preflight(pk))
    else:
        g, support = read_secret_key(args.secret_key, pk.field)
        verify_key(pk, g, support)
        out = dict(verified=True, n=pk.n, m=pk.m, t=pk.t,
                   public_key_sha256=pk.sha256)
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
