from scipy import *
from scipy.linalg import *

def solution( ) :
    print( "here are the answers of Dias Zhassanbay" )

    m1 = array( [ [ 1/2, 1/3 ], [ -2/7, 2/ 8 ]] )
    m2 = array( [ [ -1/3, 2/7 ], [ 2/5, -1/7 ]] )
    m3 = array( [ [ -1/5,2/3], [1/8, 3/11] ] )

    print( "This is part 1\n" )

    print( "the dot product of m1 and m2 is " )
    print( dot(m1,m2) )

    print( "\nThis is part 2\n" )

    print( "the inverse of m1 is ")
    print( inv(m1))

    print( "\nThis is part 3\n" )

    print( "the dot product of m1.m2 and m3 is " )
    print( dot(dot(m1,m2), m3) )

    print( "the dot product of m1 and m2.m3 is " )
    print( dot(m1, dot(m2,m3)) )

    print( "The dot products are equal due to associative property of matrix multiplication" )

    print( "\nThis is part 4\n" )

    print( "the dot product of m1 and m2 + m3 is " )
    print( dot(m1, m2 + m3))

    print( "The sum of two dot products m1.m2 and m1.m3 is " )
    print( dot(m1,m2) + dot(m1,m3) )

    print( "The dot product of m1 and m2 + m3 equals m1.m2 + m1.m3" )
    print( "Thus, matrix multiplication is distributive on the right side" )

    print( "the dot product of m1+m2 and m3 is " )
    print( dot(m1 + m2,m3) )

    print( "the sum of two dot products m1.m3 and m2.m3 is " )
    print( dot(m1,m3) + dot(m2,m3) )

    print( "The dot product of m1 + m2 and m3 equals m1.m3 + m2.m3" )
    print( "Thus, matrix multiplication is distributive on the left side" )
    print( "Therefore, matrix multiplication is distributive on both sides" )

    print( "\nThis is part 5\n" )

    v = array( [ 3.0, -1 ] )

    print( "\nThis is part 6\n" )

    print( "the product of the determinant of m1 and determinant of m2 is " )
    print( det(m1) * det(m2) )

    print( "the determinant of m1.m2 is " )
    print( det(dot(m1,m2)) )

    print( "the product of two determinants of m1 and m2 equals to the determinant of m1.m2" )
    print( "Therefore, determinant commutes over multiplication" )

    print( "\nThis is part 7\n" )

    print( "the dot product of m1 and its inverse is " )
    print( dot(m1, inv(m1)) )

    print( "the dot product of the inverse of m1 and m1 is " )
    print( dot(inv(m1), m1) )

    print( "Thus, the inverse of m1 is indeed inverse" )
