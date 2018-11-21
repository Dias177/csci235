from matrix import *
from vector import *
from rational import *

def tests():
    a11 = Rational(1,2)
    a12 = Rational(1,3)
    a21 = Rational(-2,7)
    a22 = Rational(2,8)

    b11 = Rational(-1,3)
    b12 = Rational(2, 7)
    b21 = Rational(2,5)
    b22 = Rational(-1,7)

    m1 = Matrix(a11, a12,a21,a22)
    m2 = Matrix(b11,b12,b21,b22)

    print(m1*m2)
