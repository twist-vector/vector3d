#
#
#            Nimrod's Runtime Library
#        (c) Copyright 2011 Tom Krauss
#
#    See the file "copying.txt", included in this
#    distribution, for details about the copyright.
#
# The vector3D package provides simple 3-dimensional
# vectors and matrices.  The vectors represent <x,y,z>
# coordinates and the package provides typical 3D operations
# on them.  These include add, subtract, etc.  Also included
# are 3D matrices which generally represent rotations.
#

import
  math

const
  EPS = 5.0e-12
  ## Epsilon used for approximate float comparisons

type
  TVector3D* = tuple[x, y, z: float] 
    ## a 3D vector, consisting of components x, y, and z


type 
  TMatrix3D* = object
    data: array[0..2, array[0..2, float]]
      ## a 3D matrix (3x3), which often represents a rotation matrix




proc `=~` *(a, b: TVector3D): bool =
  ## Compare two 3D vectors `a` and `b` for equality.
  result = abs(a.x-b.x)<EPS and abs(a.y-b.y)<EPS and abs(a.z-b.z)<EPS

                    
proc dot*(a: TVector3D, b: TVector3D): float {.inline.} = 
  ## Returns the inner (dot) product of two 3D vectors `a` and `b`.
  result = a.x * b.x + a.y * b.y + a.z * b.z

proc cross*(a: TVector3D, b: TVector3D): TVector3D {.inline.} = 
  ## Returns the outer (cross) product of two 3D vectors `a` and `b`.
  result.x = a.y * b.z - a.z * b.y
  result.y = a.z * b.x - a.x * b.z
  result.z = a.x * b.y - a.y * b.x

proc `*`*(a: float, b: TVector3D): TVector3D {.inline.} =
  ## Returns product of a float `a` and a 3D vector `b`.
  result.x = a * b.x
  result.y = a * b.y
  result.z = a * b.z

proc `*`*(a: TVector3D, b: float): TVector3D {.inline.} = 
  ## Returns product of a 3D vector `a` and a float `b`.
  result = b * a


proc `+`*(a: TVector3D, b: TVector3D): TVector3D {.inline.} = 
  ## Returns the sum of two 3D vectors `a` and `b`.
  result.x = a.x + b.x
  result.y = a.y + b.y
  result.z = a.z + b.z

proc `+`*(a: float, b: TVector3D): TVector3D {.inline.} = 
  ## Returns the sum of a float `a` and a 3D vectors `b` where 
  ## the float is added to each element of `b`.
  result.x = a + b.x
  result.y = a + b.y
  result.z = a + b.z

proc `+`*(a: TVector3D, b: float): TVector3D {.inline.} = 
  ## Returns the sum of a 3D vector `a` and a float `b` where 
  ## the float is added to each element of `b`.
  result = b + a


proc `-`*(a: TVector3D): TVector3D {.inline.} = 
  ## Returns the negation of the 3D vector `a`.
  result.x = -a.x
  result.y = -a.y
  result.z = -a.z

proc `-`*(a: TVector3D, b: TVector3D): TVector3D {.inline.} = 
  ## Subtracts the 3D vector `b` from the 3D vector `a` returning 
  ## the result vector.
  result = a + (-b)

proc `-`*(a: float, b: TVector3D): TVector3D {.inline.} = 
  ## Subtracts the 3D vector `b` from the float `a` returning 
  ## the resulting vector.  The float `a` is extended to be a 3D 
  ## vector with all elements equal to `a`.
  result = a + (-b)

proc `-`*(a: TVector3D, b: float): TVector3D {.inline.} = 
  ## Subtracts the float `b` from the 3D vector `a` returning 
  ## the resulting vector.  The float `b` is extended to be a 3D 
  ## vector with all elements equal to `b`.
  result = a + (-b)


proc `/`*(a: TVector3D, b: float): TVector3D {.inline.} = 
  ## Divides (scales) the vector `a` by the float `b`
  result.x = a.x/b
  result.y = a.y/b
  result.z = a.z/b




proc matrix3D*(d: array[0..8, float]): TMatrix3D = 
  var k: int = 0
  for i in countup(0,2):
    for j in countup(0,2):
      result.data[i][j] = d[k]
      k = k + 1

proc `*`*(a: TMatrix3D, b: TVector3D): TVector3D =
  ## Returns product of a 3D matrix `a` and a 3D vector `b`.
  result.x = a.data[0][0]*b[0] + a.data[0][1]*b[1] + a.data[0][2]*b[2]
  result.y = a.data[1][0]*b[0] + a.data[1][1]*b[1] + a.data[1][2]*b[2]
  result.z = a.data[2][0]*b[0] + a.data[2][1]*b[1] + a.data[2][2]*b[2]


proc `*`*(a: TMatrix3D, b: TMatrix3D): TMatrix3D = 
  ## 3D Matrix multiply
  for i in countup(0,2):
    for j in countup(0,2):
      result.data[i][j] = 0.0
      for k in countup(0,2):
        result.data[i][j] = result.data[i][j] + a.data[i][k]*b.data[k][j]



proc rotationMatrix3D*(theta: float): TMatrix3D = 
  ## Return a rotation matrix which represents a rotation about
  ## the z-axis through the angle theta. 
  result = matrix3D( [ cos(theta), sin(theta), 0.0, 
                      -sin(theta), cos(theta), 0.0, 
                          0.0,       0.0,      1.0] )

proc eulerRotationMatrix3D*(phi: float, theta: float, psi: float): TMatrix3D = 
  ## Return a rotation matrix which represents a rotation about
  ## through the three angles phi, theta, and psi.  The first rotation
  ## is through phi about the z-axis.  The second is through theta about
  ## the original x-axis.  The last is through psi about the new z-axis.
  var D = matrix3D( [ cos(phi), sin(phi), 0.0, 
                     -sin(phi), cos(phi), 0.0, 
                         0.0,     0.0,    1.0] )
                         
  var C = matrix3D( [    1.0,     0.0,       0.0,
                         0.0,  cos(theta), sin(theta), 
                         0.0, -sin(theta), cos(theta)] )
                      
  var B = matrix3D( [ cos(psi), sin(psi), 0.0, 
                     -sin(psi), cos(psi), 0.0, 
                      0.0,       0.0,     1.0] )
  result = B*(C*D)


proc norm*(a: TVector3D): float {.inline.} = 
  ## Returns the norm (length) of the 3D vector `a`.
  result = sqrt(a.x*a.x + a.y*a.y + a.z*a.z)

proc normalize*(a: TVector3D): TVector3D {.inline.} = 
  ## Returns the unit vector in the same direction as 
  ## the 3D vector `a`.
  result = a/norm(a)



proc wgs_84_norm*(ECEF_pos: TVector3D): TVector3D =
  var a = 6378137.0            # Semi-major (equatorial) axis of WGS_84 model
  var b = 6356752.314245179    # Semi-minor (polar) axis of WGS_84 model

  var norm_vector = (2.0/(a*a) * ECEF_pos.x,
                     2.0/(a*a) * ECEF_pos.y,
                     2.0/(b*b) * ECEF_pos.z)
  result = normalize(norm_vector)




proc `$`*(a: TVector3D): string =
  result = "<" & $a.x & ", " & 
                 $a.y & ", " & 
                 $a.z & ">" 


proc `$`*(a: TMatrix3D): string =
  result = "|" & $a.data[0][0] & ", " & $a.data[0][1] & ", " & $a.data[0][2] & "|\n" &
           "|" & $a.data[1][0] & ", " & $a.data[1][1] & ", " & $a.data[1][2] & "|\n" & 
           "|" & $a.data[2][0] & ", " & $a.data[2][1] & ", " & $a.data[2][2] & "|\n"



when isMainModule:
  var a     : TVector3D = (1.0, 2.0, 3.0)
  var b     : TVector3D = (-1.0, 2.0, -3.0)
  
  assert( a =~ (1.0,2.0,3.0) )
  assert( b =~ (-1.0,2.0,-3.0) )
  
  assert( dot(a,b) == -6.0 )
  assert( dot(b,a) == -6.0 )
  assert( cross(a,b) =~ (-12.0, 0.0, 4.0) )
  assert( cross(b,a) =~ (12.0, 0.0, -4.0) )
  assert( (10.0*a) =~ (10.0, 20.0, 30.0) )
  assert( (10.0*a) =~ (a*10.0) )

  assert( (-a) =~ (-1.0, -2.0, -3.0) )
  assert( (a+10.0) =~ (11.0, 12.0, 13.0) )
  assert( (10.0+a) =~ (11.0, 12.0, 13.0) )
  assert( (a+b) =~ (0.0, 4.0, 0.0) )
  assert( (a-b) =~ (2.0, 0.0, 6.0) )
  assert( (a-1.0) =~ (0.0, 1.0, 2.0) )
  assert( (1.0-a) =~ (0.0, -1.0, -2.0) )
  
  var R = rotationMatrix3D( 90.0*pi/180.0 )
  var uv: TVector3D = (1.0, 0.0, 0.0)
  assert( (R*uv) =~ (0.0, -1.0, 0.0) )
