#!/usr/bin/env python3
"""Shared field-conversion and verified key-export helpers for hover.py.
"""
import json
import os
from pathlib import Path

from tii_cpu import verify_key

__all__ = ["integer_element", "int_coefficients", "atomic_json", "write_key"]


def integer_element(x):
    return sum(int(c) << j for j, c in enumerate(x.polynomial()))


def int_coefficients(g):
    return [integer_element(c) for c in g.list()]


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(path)


def write_key(path, pk, result, method):
    gg, support = result
    verify_key(pk, gg, support)
    obj = dict(
        g=gg, support=support, m=pk.m, t=pk.t, n=pk.n,
        modulus=pk.modulus_coefficients,
        public_key_sha256=pk.sha256, method=method, verified=True,
    )
    atomic_json(path, obj)

    def text_element(x):
        terms = [
            '1' if i == 0 else 'a' if i == 1 else 'a^{}'.format(i)
            for i in range(pk.m) if (x >> i) & 1
        ]
        return ' + '.join(reversed(terms)) if terms else '0'

    text = '\n'.join(
        '[' + ', '.join(text_element(v) for v in array) + ']'
        for array in (gg, support)
    ) + '\n'
    Path(str(path) + '.txt').write_text(text)
    print('VERIFIED key written to {}'.format(path), flush=True)
