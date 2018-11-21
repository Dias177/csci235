from numbers import *

def gcd(n1, n2):
    if n1 == 0 and n2 == 0:
        raise ArithmeticError( "gcd(0,0) does not exist" )
    n1 = abs(n1)
    n2 = abs(n2)
    while n1 != 0 and n2 != 0:
        if n1 > n2:
            n1 = n1%n2
        else:
            n2 = n2 % n1
    return n1 + n2

class Rational(Number):
    def __init__(self, num, denum = 1) :
        self.num = num
        self.denum = denum

        self.normalize()

    def normalize(self):
        common = gcd(self.denum, self.num)

        if self.denum == 0:
            raise ZeroDivisionError( "cannot divide by 0")
        self.num = self.num // common
        self.denum = self.denum // common

        if self.denum < 0:
            self.denum *= -1
            self.num *= -1

    def __repr__(self):
        if self.denum == 1:
            return "{}".format(self.num)

        return "{}/{}".format(self.num,self.denum)

    def __neg__(self):
        self.num *= -1
        return self

    def __add__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = self.num*new.denum + self.denum*new.num
        denum = self.denum*new.denum

        return Rational(num, denum)

    def __sub__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = self.num*new.denum - self.denum*new.num
        denum = self.denum*new.denum

        return Rational(num, denum)

    def __radd__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = self.num*new.denum + self.denum*new.num
        denum = self.denum*new.denum

        return Rational(num, denum)

    def __rsub__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = -self.num*new.denum + self.denum*new.num
        denum = self.denum*new.denum

        return Rational(num, denum)

    def __mul__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, int):
            return Rational(self.num*other, self.denum)

        new = Rational(other.num, other.denum)
        return Rational(self.num*new.num, self.denum*new.denum)

    def __rmul__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = self.num*new.num
        denum = self.denum*new.denum

        return Rational(num, denum)

    def __truediv__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = self.num*new.denum
        denum = self.denum*new.num

        return Rational(num, denum)

    def __rtruediv__(self, other):
        if isinstance(other, str):
            raise NotImplementedError

        if isinstance(other, Rational):
            new = Rational(other.num, other.denum)

        if isinstance(other, int):
            new = Rational(other, 1)

        num = self.denum*new.num
        denum = self.num*new.denum

        return Rational(num, denum)
