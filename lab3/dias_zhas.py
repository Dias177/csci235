from matrix import *
from vector import *
from rational import *

def tests():
    a11 = Rational(1, 2)
    a12 = Rational(1, 3)
    a21 = Rational(-2, 7)
    a22 = Rational(2, 8)

    b11 = Rational(-1, 3)
    b12 = Rational(2, 7)
    b21 = Rational(2, 5)
    b22 = Rational(-1, 7)

    m1 = Matrix(a11, a12, a21, a22)
    m2 = Matrix(b11, b12, b21, b22)

    print("matrix m1\n{}".format(m1))
    print("matrix m2\n{}".format(m2))

    print("m1*m2\n{}".format(m1@m2))

    m3 = m1.inverse()
    print("inverse of m1, matrix m3\n{}".format(m3))

    mm1 = (m1@m2)@m3 - m1@(m2@m3)
    print("(m1*m2)*m3 - m1*(m2*m3)\n{}".format(mm1))

    mm2 = m1@(m2+m3) - m1@m2 - m1@m3
    print("m1*(m2+m3) - m1*m2 - m1*m3\n{}".format(mm2))

    mm3 = (m1+m2)@m3 - m1@m3 - m2@m3
    print("(m1+m2)*m3 - m1*m3 - m2*m3\n{}".format(mm3))

    v = Vector(1,1)
    mm4 = m1(m2(v)) - (m1@m2)(v)
    print("m1(m2(v)) - (m1*m2)(v)\n{}".format(mm4))

    mm5 = m1.determinant()*m2.determinant() - (m1@m2).determinant()
    print("det(m1)*det(m2) - det(m1*m2)\n{}".format(mm5))

    print("m1*inv(m1)\n{}".format(m1@m3))
    print("inv(m1)*m1\n{}".format(m3@m1))
