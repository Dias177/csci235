
from scipy import *
from scipy.linalg import *
import timeit

# runge_kutta1 is Euler's method:

def runge_kutta1( F, x, h ) :
   k1 = F( x )

   return x + h * k1

def runge_kutta21( F, x, h ) :
   k1 = F( x )
   k2 = F( x + h*k1 )

   return x + h * ( 0.5 * k1 + 0.5 * k2 )

def runge_kutta22( F, x, h ) :
   k1 = F( x )
   k2 = F( x + 0.5 * h * k1 )

   return x + h * k2

def runge_kutta31( F, x, h ) :
   k1 = F( x )
   k2 = F( x + h * (2.0/3.0) * k1 )
   k3 = F( x + h * ( (1.0/3.0) * k1 + (1.0/3.0) * k2 ) )

   return x + h * ( (1.0/4.0) * k1 + (3.0/4.0) * k3 )

def runge_kutta41( F, x, h ) :
   k1 = F( x )
   k2 = F( x + 0.5 * h * k1 )
   k3 = F( x + 0.5 * h * k2 )
   k4 = F( x + h * k3 )

   return x + (1.0/6.0) * h * (k1 + 2.0*k2 + 2.0*k3 + k4 )

def runge_kutta4_kuntzmann( F, x, h) :
   k1 = F( x )
   k2 = F( x + h * (2.0/5.0) * k1 )
   k3 = F( x + h * ( ( -3.0/20.0) * k1 + (3.0/4.0) * k2 ) )
   k4 = F( x + h * ( (19.0/44.0) * k1 + (-15.0/44.0) * k2 + (40.0/44.0) * k3 ) )
   return x + h * ( (55.0/360.0) * k1 + (125.0/360.0) * k2 + (125.0/360.0) * k3 + (55.0/360.0) * k4 )

def runge_kutta5( F, x, h) :
   k1 = F( x )
   k2 = F( x + h * (1.0/4.0) * k1 )
   k3 = F( x + h * ( (1.0/8.0) * k1 + (1.0/8.0) * k2 ) )
   k4 = F( x + h * (1.0/2.0) * k3 )
   k5 = F( x + h * ( (3.0/16.0) * k1 + (-3.0/8.0) * k2 + (3.0/8.0) * k3 + (9.0/16.0) * k4 ) )
   k6 = F( x + h * ( (-3.0/7.0) * k1 + (8.0/7.0) * k2 + (6.0/7.0) * k3 + (-12.0/7.0) * k4 + (8.0/7.0) * k5 ) )

   return x + h * (   ( 7.0 / 90.0 ) * k1 + ( 32.0 / 90.0 ) * k3 + ( 12.0 / 90.0 ) * k4 + ( 32.0 / 90.0 ) * k5 + ( 7.0 / 90.0 ) * k6 )

def approx( h = 4.0E-3 ) :
   print( "testing Runge-Kutta methods on the catenary" )

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta41( cat, p, h )
      x += h

   p = runge_kutta41( cat, p, x1 - x )
   x = x1

   print( "h = ", h )
   print( "final value = ", p )

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   print( "error = ", error )

def approx1( h = 1.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta1( cat, p, h )
      x += h

   p = runge_kutta1( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx12( h = 2.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta1( cat, p, h )
      x += h

   p = runge_kutta1( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx13( h = 4.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta1( cat, p, h )
      x += h

   p = runge_kutta1( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx2( h = 1.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta21( cat, p, h )
      x += h

   p = runge_kutta21( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx22( h = 2.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta21( cat, p, h )
      x += h

   p = runge_kutta21( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx23( h = 4.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta21( cat, p, h )
      x += h

   p = runge_kutta21( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx3( h = 1.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta41( cat, p, h )
      x += h

   p = runge_kutta41( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx32( h = 2.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta41( cat, p, h )
      x += h

   p = runge_kutta41( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def approx33( h = 4.0E-3 ) :

   x0 = 0.0
   x1 = 1.0

   mu = 2.0

   s0 = array( [ 1.0 / mu, 0.0 ] )

   p = s0
   x = x0

   def cat( p ) :
      return array( [ p[1], mu * sqrt( 1.0 + p[1] * p[1] ) ] )

   while x + h < x1 :
      p = runge_kutta41( cat, p, h )
      x += h

   p = runge_kutta41( cat, p, x1 - x )
   x = x1

   expected = array( [ cosh( mu * x1 ) / mu, sinh( mu * x1 ) ] )
   error = p - expected

   return error

def table():
    print(" Value of h  | 0.001", " " * 27, " | 0.002", " " * 27, " | 0.004")
    print(" Euler ", " "*4,"|", approx1(), " "*8, "|", approx12(), " "*8, "|", approx13())
    print(" Heun ", " " * 5, "|", approx2(), " |", approx22(), " |", approx23())
    print(" Runge-Kutta", "|", approx3(), " |", approx32(), " |", approx33())
# 1) Time was tested on the Euler method with h = 0.001:
# C++
# real 0m0.005s
# user 0m0.002s
# sys 0m0.002s
#
# python
# n = 100, t = 0.01277 s

# 2) Time was tested on the Heun method with h = 0.002:
# C++
# real 0m0.005s
# user 0m0.002s
# sys 0m0.002s
#
# python
# n = 100, t = 0.01398 s

# 3) Time was tested on the Runge-Kutta method with h = 0.004:
# C++
# real 0m0.005s
# user 0m0.002s
# sys 0m0.002s

# python
# n = 100, t = 0.01297 s
