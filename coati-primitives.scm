#|
Copyright (c) 2014 Richard van Roy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Large parts of this egg are ported from Chipmunk2D's cpVect.h (c) 2007 - Scott Lembcke and Howling Moon Software.
|#

(module coati-primitives
  *
  (import chicken scheme foreign)
  (use data-structures srfi-1 srfi-4)

#>
#include <float.h>
#include "triangulate.h"
<#


;-------------------------------------------------------
; %
;-------------------------------------------------------

(define (%wrap-degree v)
  (if (negative? v) (+ 360 v) v))

(define (%f64vector-part v size)
  (assert (zero? (modulo (f64vector-length v) size)))
  (let loop ((r (list)) (n (f64vector-length v)))
    (if (= n 0) r
	(loop (cons (subf64vector v (- n size) n) r)
	      (- n size)))))

;-------------------------------------------------------
; Float
;-------------------------------------------------------

(define fmod (foreign-lambda double "fmod" double double))

(define (clamp f mmin mmax)
  (min (max f mmin) mmax))

(define (sqr x) (* x x))

(define double-min (foreign-value "DBL_MIN" double))

(define infinity (foreign-value "INFINITY" double))

;-------------------------------------------------------
; Constants
;-------------------------------------------------------

(define epsilon 1e-6)

(define pi (acos -1.0))

(define pi/2 (/ pi 2.0))

(define 2pi (* 2.0 pi))

(define -pi (- pi))

(define 360/2pi (/ 360.0 2pi))

(define 2pi/360 (/ 2pi 360.0))

;-------------------------------------------------------
; Vectors
;-------------------------------------------------------

; Returs a new vector
(define (create-vect x y)
  (f64vector x y))

(define (vect-x v)
  (f64vector-ref v 0))

(define (vect-y v)
  (f64vector-ref v 1))

; Constant for the zero vector.
(define (vect-zero)
  (create-vect 0 0))

; Check if two vectors are equal. (Be careful when comparing floating point numbers!)
(define (vect=? a b)
  (and (equal? (vect-x a) (vect-x b))
       (equal? (vect-x a) (vect-y b))))

; Add two vectors.
(define (vect+ a b)
  (create-vect (+ (vect-x a) (vect-x b))
	     (+ (vect-y a) (vect-y b))))

; Subtract two vectors or negate a vector.
(define (vect- a #!optional b)
  (if b (create-vect (- (vect-x a) (vect-x b))
		   (- (vect-y a) (vect-y b)))
      (create-vect (- (vect-x a))
		 (- (vect-y a)))))

; Scalar multiplication.
(define (vect* v s)
  (create-vect (* (vect-x v) s)
	     (* (vect-y v) s)))

; Vector dot product.
(define (vect-dot a b)
  (+ (* (vect-x a)
	(vect-x b))
     (* (vect-y a)
	(vect-y b))))

; 2D vector cross product analog.
; The cross product of 2D vectors results in a 3D vector with only a z component.
; This function returns the magnitude of the z value.
(define (vect-cross a b)
  (- (* (vect-x a)
	(vect-y b))
     (* (vect-y a)
	(vect-x b))))
; Returns a perpendicular vector. (90 degree rotation)
(define (vect-perp v)
  (create-vect (- (vect-y v)) (vect-x v)))

; Returns a perpendicular vector. (-90 degree rotation)
(define (vect-vperp v)
  (create-vect (vect-y v) (- (vect-x v))))

; Returns the vector projection of /a/ onto /b/.
(define (vect-project a b)
  (vect* a (/ (vect-dot a b)
	      (vect-dot b b))))

; Returns the unit length vector for the given angle (in radians).
(define (angle->vect a)
  (create-vect (cos a) (sin a)))

; Returns the angular direction v is pointing in (in radians).
(define (vect->angle v)
  (atan (vect-y v) (vect-x v)))

; Uses complex number multiplication to rotate /a/ by /b/. Scaling will occur if /a/ is not a unit vector.
(define (vect-rotate a b)
  (create-vect (+ (* (vect-x a) (vect-x b))
		(* (vect-y a) (vect-y b)))
	     (- (* (vect-x a) (vect-y b)
		   (vect-y a) (vect-x b)))))

; Inverse of vect-rotate
(define (vect-unrotate a b)
  (create-vect (+ (* (vect-x a) (vect-x b))
		(* (vect-y a) (vect-y b)))
	     (- (* (vect-y a) (vect-x b)
		   (vect-x a) (vect-y b)))))

; Returns the squared length of v. Faster than cpvlength() when you only need to compare lengths.
(define (vect-length-squared v)
  (vect-dot v v))

; Returns the length of v.
(define (vect-length v)
  (sqrt (vect-dot v v)))

; Linearly interpolate between /a/ and /b/.
(define (vect-lerp v1 v2 t)
  (vect+ (vect* v1 (- 1.0 t))
	 (vect* v2 t)))

; Returns a normalized copy of v.
(define (vect-normalize v)
  (vect* v (/ 1.0 (+ (vect-length v) double-min))))

; Clamp v to length len.
(define (vect-clamp v len)
  (if (> (vect-dot v v) (sqr len))
      (vect* (vect-normalize v) len)
      v))

; Linearly interpolate between v1 towards v2 by distance d.
(define (vect-lerp-const v1 v2 dist)
  (vect+ v1 (vect+ (vect-clamp v2 v1) dist)))

; Returns the distance between v1 and v2.
(define (vect-dist v1 v2)
  (vect-length (vect- v1 v2)))

; Returns the squared distance between v1 and v2. Faster than cpvdist() when you only need to compare distances.
(define (vect-dist-squared v1 v2)
  (vect-length-squared (vect- v1 v2)))

; Returns true if the distance between v1 and v2 is less than dist.
(define (vect-near? a b dist)
  (< (vect-dist-squared a b) (sqr dist)))

; Spherical linearly interpolate between /a/ and /b/.
(define (vect-spherical-lerp a b t)
  (let* ((dot (vect-dot (vect-normalize a) (vect-normalize b)))
	 (omega (clamp dot -1.0 1.0)))
    (if (< omega 0.001)
	(vect-lerp a b t)
	(let ((denom (/ 1.0 (sin omega))))
	  (vect+ (vect* a (* (sin (* (- 1.0 t) omega)) denom))
		 (vect* b (* (sin (* (* t omega) denom)))))))))

; Spherical linearly interpolate between /a/ towards /b/ by no more than angle /angle/ in radians.
(define (vect-spherical-lerp-const a b angle)
  (let* ((dot (vect-dot (vect-normalize a) (vect-normalize b)))
	 (omega (clamp dot -1.0 1.0)))
    (vect-spherical-lerp a b (/ (min angle omega) omega))))

;-------------------------------------------------------
; Bounding Boxes
;-------------------------------------------------------

; Returs a new bounding box.
(define (create-bb l b r t)
  (f64vector l b r t))

(define (bb-l bb)
  (f64vector-ref bb 0))

(define (bb-b bb)
  (f64vector-ref bb 1))

(define (bb-r bb)
  (f64vector-ref bb 2))

(define (bb-t bb)
  (f64vector-ref bb 3))

; Constructs a /bb/ for a circle with the given position and radius.
(define (create-bb-for-circle p r)
  (create-bb (- (vect-x p) r)
	   (- (vect-y p) r)
	   (+ (vect-x p) r)
	   (+ (vect-y p) r)))

; Returns true if /a/ and /b/ intersect.
(define (bb-intersects? a b)
  (and (<= (bb-l a) (bb-r b))
       (<= (bb-l b) (bb-r a))
       (<= (bb-b a) (bb-t b))
       (<= (bb-b b) (bb-t a))))

; Returns true if /other/ lies completely within /bb/.
(define (bb-contains? bb other)
  (and (<= (bb-l bb) (bb-l other))
       (>= (bb-r bb) (bb-r other))
       (<= (bb-b bb) (bb-b other))
       (>= (bb-t bb) (bb-t other))))

; Returns true if /bb/ contains /v/.
(define (bb-constains-vect? bb v)
  (and (<= (bb-l bb) (vect-x v))
       (>= (bb-r bb) (vect-x v))
       (<= (bb-b bb) (vect-y v))
       (>= (bb-t bb) (vect-y v))))

; Returns a bounding box that holds both bounding boxes.
(define (bb-merge a b)
  (create-bb (min (bb-l a) (bb-l b))
	   (min (bb-b a) (bb-b b))
	   (max (bb-r a) (bb-r b))
	   (max (bb-t a) (bb-t b))))

; Returns a bounding box that holds both /bb/ and /v/.
(define (bb-expand bb v)
  (create-bb (min (bb-l bb) (vect-x v))
	   (min (bb-b bb) (vect-x v))
	   (max (bb-r bb) (vect-y v))
	   (max (bb-t bb) (vect-y v))))

; Returns the center of a bounding box.
(define (bb-center bb)
  (vect-lerp (create-vect (bb-l bb) (bb-b bb))
	     (create-vect (bb-r bb) (bb-t bb))
	     0.5))

; Returns the area of the bounding box.
(define (bb-area bb)
  (* (- (bb-r bb) (bb-l bb))
     (- (bb-t bb) (bb-b bb))))

; Merges /a/ and /b/ and returns the area of the merged bounding box.
(define (bb-merged-area a b)
  (* (- (max (bb-r a) (bb-r b))
	(min (bb-l a) (bb-l b)))
     (- (max (bb-t a) (bb-t b))
	(min (bb-b a) (bb-b b)))))

; Returns the fraction along the segment query the bounding box is hit. Returns /infinity/ if it doesn't hit.
(define (bb-segment-query bb a b)

  (let* ((idx (/ 1 (- (vect-x b) (vect-x a))))
	 (tx1 (if (= (bb-l bb) (vect-x a))
		  (- infinity)
		  (* (- (bb-l bb) (vect-x a)) idx)))
	 (tx2 (if (= (bb-r bb) (vect-x a))
		  (- infinity)
		  (* (- (bb-r bb) (vect-x a)) idx)))
	 (txmin (min tx1 tx2))
	 (txmax (max tx1 tx2))
	 ;
	 (idy (/ 1 (- (vect-y b) (vect-y a))))
	 (ty1 (if (= (bb-b bb) (vect-y a))
		  (- infinity)
		  (* (- (bb-b bb) (vect-y a)) idy)))
	 (ty2 (if (= (bb-t bb) (vect-y a))
		  (- infinity)
		  (* (- (bb-t bb) (vect-y a)) idy)))
	 (tymin (min ty1 ty2))
	 (tymax (max ty1 ty2)))
    (if (and (<= tymin txmax)
	     (<= txmin tymax))
	(let ((mmin (max txmin tymin))
	      (mmax (min txmax tymax)))
	  (if (and (<= 0.0 mmax)
		   (<= mmin 1.0))
	      (max mmin 0.0)
	      infinity))
	infinity)))

; Return true if the bounding box intersects the line segment with ends /a/ and /b/.
(define (bb-intersects-segment? bb a b)
  (not (= (bb-segment-query bb a b) infinity)))

(define (bb-clamp-vect bb v)
  (create-vect (clamp (vect-x v) (bb-l bb) (bb-r bb))
	     (clamp (vect-y v) (bb-b bb) (bb-t bb))))

;-------------------------------------------------------
; Lines
;-------------------------------------------------------

; Makes a line from two vectors
(define (create-line a b)
  (f64vector (vect-x a) (vect-y a)
	     (vect-x b) (vect-y b)))

;-------------------------------------------------------
; Polygon
;-------------------------------------------------------

; Creates a new polygon from a list of vectors.
(define-syntax create-polygon
  (syntax-rules ()
    ((_  vects)
     (list->f64vector (append-map f64vector->list vects)))))

; Triangulates the given polygon and returns a list of triangles.
(define (polygon-triangulate polygon)
  (let* ((return-size (- (* (f64vector-length polygon) 3) 12))
         (res (create-f64vector return-size)))
    ((foreign-lambda* void ((f64vector polygon)
                            (f64vector res)
                            (integer polygonSize)
                            (integer returnSize)) "
	std::vector<Vector> vec;
	vec.resize(polygonSize/2);
	memcpy(vec.data(), polygon, sizeof(Vector)*(polygonSize/2));
	std::vector<Triangle> tmp = triangulate(vec);
	memcpy(res, tmp.data(), returnSize*sizeof(Vector));")
     polygon res (f64vector-length polygon) return-size)
    ; Chop the result into a list of triangles.
    (%f64vector-part res 3)))

; Return #t if the given polygon is convex.
(define (polygon-convex? polygon)
  ((foreign-lambda* bool ((f64vector polygon)
                          (unsigned-integer length)) "
	C_return( isConvex((Vector*)polygon, length) );")
   polygon (f64vector-length polygon)))

; Converts a polygon to a list of vertices.
(define (polygon->vects polygon)
  (%f64vector-part polygon 2))

(define (%sort-vects vects)
  (sort vects
	(lambda (a b) 
	    (or (< (vect-x a) (vect-x b))
		(= (vect-x a) (vect-y b))))))

(define (%cross o a b)
  (- (* (- (vect-x a)
	 (vect-x o))
      (- (vect-y b)
	 (vect-y o)))
     (* (- (vect-y a)
	   (vect-y o))
	(- (vect-x b)
	   (vect-x o)))))

; Returns the convex hull of a group of vertices in clockwise order.
(define (convex-hull vects)
  (let* ((sorted (%sort-vects vects))
	 (lower (list))
	 (upper (list)))
    (if (<= (length vects) 1) vects
	(begin
	 (map (lambda (v)
		(let loop ()
		  (when (and (>= (length lower) 2)
			     (<= (%cross (cadr lower) (car lower) v) 0))
		    (set! lower (cdr lower))
		    (loop)))
		(set! lower (cons v lower)))
	      vects)
	 (map (lambda (v)
	 	(let loop ()
	 	  (when (and (>= (length upper) 2)
	 		     (<= (%cross (cadr upper) (car upper) v) 0))
		    (set! upper (cdr upper))
	 	    (loop)))
	 	(set! upper (cons v upper)))
	      (reverse vects))
	 (reverse (append (cdr lower) (cdr upper)))))))

; Converts any polygon to a convex polygon.
(define (polygon-convex-hull vects)
  (convex-hull (polygon->vects vects)))

;-------------------------------------------------------
; Colour
;-------------------------------------------------------

; Creates a new RGB colour
(define (create-rgb r g b #!optional (a 1.0))
  (f64vector r g b a))

(define (rgb-r rgb)
  (f64vector-ref rgb 0))

(define (rgb-g rgb)
  (f64vector-ref rgb 1))

(define (rgb-b rgb)
  (f64vector-ref rgb 2))

(define (rgb-a rgb)
  (f64vector-ref rgb 3))

(define (rgb->hsv rgb)
  (let* ((r (rgb-r rgb))
	 (g (rgb-g rgb))
	 (b (rgb-b rgb))
	 (a (rgb-a rgb))
	 (mmin (min r g b))
	 (mmax (max r g b))
	 (c (- mmax mmin)))
    (if (not (= c 0.0))
	(let* ((v mmax)
	       (s (/ c v)))
	  (cond ((= mmax r)
		 (create-hsv (%wrap-degree
			    (* (fmod (/ (- g b) c) 6.0) 60.0)) s v a))
		((= mmax g)
		 (create-hsv (%wrap-degree
			    (* (+ (/ (- b r) c) 2.0) 60.0)) s v a))
		(else
		 (create-hsv (%wrap-degree
			    (* (+ (/ (- r g) c) 4.0) 60.0)) s v a))))
	(create-hsv 0 0 0 a))))

; Creates a new HSV colour
(define (create-hsv h s v #!optional (a 1.0))
  (f64vector h s v a))

(define (hsv-h hsv)
  (f64vector-ref hsv 0))

(define (hsv-s hsv)
  (f64vector-ref hsv 1))

(define (hsv-v hsv)
  (f64vector-ref hsv 2))

(define (hsv-a hsv)
  (f64vector-ref hsv 3))

(define (hsv->rgb hsv)
  (let* ((h (hsv-h hsv))
	 (s (hsv-s hsv))
	 (v (hsv-v hsv))
	 (a (hsv-a hsv))
	 (c (* v s))
	 (m (- v c))
	 (x (* c (- 1.0 (abs (- (fmod (/ h 60.0) 2) 1)))))
	 (m (- v c)))
    (cond
     ((and (>= h 0.0)
	   (<  h 60.0))  (create-rgb (+ c m) (+ x m) m a))
     ((and (>= h 60.0)
	   (<  h 120.0)) (create-rgb (+ x m) (+ c m) m a))
     ((and (>= h 120.0)
	   (<  h 180.0)) (create-rgb m (+ c m) (+ x m) a))
     ((and (>= h 180.0)
	   (<  h 240.0)) (create-rgb m (+ x m) (+ c m) a))
     ((and (>= h 240.0)
	   (<  h 300.0)) (create-rgb (+ x m) m (+ c m) a))
     ((and (>= h 300.0)
	   (<  h 360.0)) (create-rgb (+ c m) m (+ x m) a))
     (else (create-rgb m m m a)))))
)
