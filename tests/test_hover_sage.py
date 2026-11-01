"""HOVER operator, geometry, checkpoints and public-only recovery tests."""
import itertools
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
try:
    import sage.all
    HAVE_SAGE=True
except ImportError:
    HAVE_SAGE=False


@unittest.skipUnless(HAVE_SAGE,'requires SageMath and a C compiler')
class HoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from sage.all import GF,PolynomialRing,set_random_seed
        from tii_cpu import BinaryField,binary_rref,parity_rows,read_public_key
        from recover_key import integer_element,int_coefficients
        cls.temp=tempfile.TemporaryDirectory()
        cls.root=Path(cls.temp.name)
        cls.F=GF(128,'a',modulus=PolynomialRing(GF(2),'z')([1,1,0,0,0,0,0,1]))
        R=PolynomialRing(cls.F,'X')
        set_random_seed(120)
        cls.secret_polynomial=R.irreducible_element(4)
        cls.secret_support=list(cls.F)
        random.Random(32).shuffle(cls.secret_support)
        rows=list(binary_rref(parity_rows([integer_element(x) for x in cls.secret_support],
                  int_coefficients(cls.secret_polynomial),BinaryField(0x83))))
        path=cls.root/'public.txt'
        path.write_text('\n'.join('['+' '.join(str((r>>j)&1) for j in range(128))+']'
                                  for r in rows)+'\n[1, 1, 0, 0, 0, 0, 0, 1]\n')
        cls.pk=read_public_key(path)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_native_operators_against_independent_definition(self):
        from sage.all import GF,matrix,random_matrix,set_random_seed
        from hover import native
        set_random_seed(42)
        A=random_matrix(GF(2),7,3)
        for p in (3,4):
            monomials=list(itertools.combinations(range(7),p))
            orders=[1] if p==3 else [0,2]
            jets=[s for u in orders for s in itertools.combinations(range(7),u)]
            expected=matrix(GF(2),len(jets)*A.ncols(),len(monomials))
            for j in range(A.ncols()):
                for s,S in enumerate(jets):
                    for t,T in enumerate(monomials):
                        if set(S)<=set(T):
                            expected[j*len(jets)+s,t]=all(A[x,j] for x in T if x not in S)
            actual=native().construct(A,p)
            self.assertEqual(actual,expected)
            standard=expected.right_kernel_matrix()
            K=native().kernel_inplace(actual,standard.nrows())
            self.assertEqual(K.row_space(),standard.row_space())
            self.assertTrue(native().replay(A,K,p))
            col=next(i for i in range(expected.ncols()) if not expected.column(i).is_zero())
            K[0,col]+=1
            self.assertFalse(native().replay(A,K,p))

    def test_scroll_and_grs_recovery(self):
        from sage.all import matrix,random_matrix,set_random_seed
        from hover import scroll_support,grs_support
        F=self.F
        set_random_seed(333)
        alpha=list(F)
        multipliers=[F.random_element() or F.one() for _ in alpha]
        C=matrix(F,[[u*x**i for u,x in zip(multipliers,alpha)] for i in range(7)])
        result=grs_support(C)
        self.assertEqual(len(set(result)),128)
        w=[[F.random_element() for _ in range(3)] for x in alpha]
        V=matrix(F,[[multipliers[j]*alpha[j]**i*w[j][h]
                     for j in range(128)] for h in range(3) for i in range(3)])
        S=random_matrix(F,9,9)
        while not S.is_invertible(): S=random_matrix(F,9,9)
        visible,result=scroll_support(S*V)
        self.assertEqual(visible,list(range(128)))
        self.assertEqual(len(set(result)),128)

    def test_partial_support_with_nonaffine_pole(self):
        from hover import finish_partial
        from tii_cpu import verify_key
        H,_=self.pk.sage_objects()
        visible=list(range(5,128))
        pole=self.secret_support[0]
        beta=[1/(self.secret_support[j]+pole) for j in visible]
        result=finish_partial(H,self.pk,visible,beta)
        self.assertTrue(verify_key(self.pk,*result))

    def test_degree_four_scroll_geometry(self):
        from sage.all import GF,matrix,random_matrix,set_random_seed
        from hover import scroll_support
        F=GF(256,'a')
        alpha=list(F)[:150]
        set_random_seed(92)
        w=[[F.random_element() for _ in range(4)] for x in alpha]
        V=matrix(F,[[alpha[j]**i*w[j][h] for j in range(len(alpha))]
                     for h in range(4) for i in range(4)])
        S=random_matrix(F,16,16)
        while not S.is_invertible(): S=random_matrix(F,16,16)
        visible,result=scroll_support(S*V,4)
        self.assertEqual(len(visible),150)
        self.assertEqual(len(set(result)),150)

    def test_public_only_hover_and_portable_checkpoint(self):
        from hover import recover,profile,read_kernel
        from tii_cpu import verify_key
        cache=self.root/'cache'
        result=recover(self.pk,3,3,0,cache,0.1)
        self.assertTrue(verify_key(self.pk,*result))
        files=list(cache.glob('hov_*.json'))
        path=next(p for p in files if '.certificate.' not in p.name)
        K=read_kernel(path,self.pk,profile(self.pk,3))
        self.assertEqual(K.nrows(),self.pk.m)
        obj=json.loads(path.read_text())
        obj['public_key_sha256']='wrong'
        bad=self.root/'bad.json'
        bad.write_text(json.dumps(obj))
        with self.assertRaises(ValueError): read_kernel(bad,self.pk,profile(self.pk,3))


if __name__=='__main__':
    unittest.main()
