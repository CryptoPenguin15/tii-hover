import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tii_cpu import (BinaryField, PublicKey, binary_rref, irreducible,
                     parity_rows, parse_field_expression, read_public_key,
                     verify_key)


class Tests(unittest.TestCase):
    def test_real_input(self):
        p = read_public_key(Path(__file__).resolve().parents[1] / 'data/pk_McEliece_100.txt')
        self.assertEqual((p.m,p.t,p.n,p.k), (8,8,254,190))
        self.assertEqual(p.sha256, '276ec7ccfcb69885537fe8745dd676c028939ff6ff619d662ff7111d3c77bd0d')

    def test_field(self):
        f = BinaryField(0x11d)
        for x in range(1,256):
            self.assertEqual(f.mul(x,f.inverse(x)), 1)
        self.assertTrue(irreducible([1,0,1,1,1,0,0,0,1], BinaryField(3)))
        self.assertFalse(irreducible([1,0,1], BinaryField(3)))

    def test_row_space(self):
        self.assertEqual(binary_rref([3,5,6]), binary_rref([6,3]))
        self.assertNotEqual(binary_rref([3,5]), binary_rref([3]))

    def test_generated_key_and_corruption(self):
        f = BinaryField(0b1011)
        g = next([a,b,1] for a in range(8) for b in range(8)
                 if irreducible([a,b,1],f))
        support = list(range(8))
        rows = list(binary_rref(parity_rows(support,g,f)))
        self.assertEqual(len(rows),6)
        pk = PublicKey(rows,8,[1,1,0,1],'synthetic')
        self.assertTrue(verify_key(pk,g,support))
        with self.assertRaises(ValueError):
            verify_key(pk,g,[0]*8)
        corrupt = PublicKey([rows[0]^1]+rows[1:],8,[1,1,0,1],'synthetic')
        with self.assertRaises(ValueError):
            verify_key(corrupt,g,support)

    def test_strict_parser(self):
        for text in ('[1 0]\n[0 1]\n[1, 0, 1]\n',
                     '[1 2]\n[1, 1]\n',
                     '[1 0]\n[0]\n[1, 1]\n'):
            with tempfile.TemporaryDirectory() as td:
                p = Path(td)/'bad.txt'
                p.write_text(text)
                with self.assertRaises(ValueError):
                    read_public_key(p)
        with self.assertRaises(ValueError):
            parse_field_expression('__import__("os")',BinaryField(0x11d))


if __name__ == '__main__':
    unittest.main()
