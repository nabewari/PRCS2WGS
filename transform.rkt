#lang racket
(require math)
(require math/base)
(require srfi/1)
;input x y
(define x 11573.375)
(define y 22694.980)
;直交座標系の原点
(define lo0 0)
(define la0 0)
;他のパラメータ
(define F (/ 1 298.257222101))
(define alfa 6378137)
(define m0 0.9999)
(define n (/ 1 (- (* 2 F) 1)))
;秒単位
(define p (* 3600 (/ 180 pi)))

;utils
(define(pow x n)
  (apply * (make-list n x)))

(define (sigma i j func)
  (apply + (map func (iota j i 1))))
;delta A Bを関数にしておく
(define (d j)
  (cond ((= j 1)
         (-
         (+
         (+
         (-
         (- (* 2 n) (* (/ 2 3) (pow n 3)))
         (* 2 (pow n 3)))
         (* (/ 116 45) (pow n 4)))
         (* (/ 26 45) (pow n 5)))
         (* (/ 2854 675) (pow n 6))))
        ((= j 2)
         (+
         (+
         (-
         (- (* (/ 7 3) (pow n 2)) (* (/ 8 5) (pow n 3)))
         (* (/ 227 45) (pow n 4)))
         (* (/ 2704 315) (pow n 5)))
         (* (/ 2323 945) (pow n 6))))
        ((= j 3)
         (+
         (-
         (- (* (/ 56 15) (pow n 3)) (* (/ 136 35) (pow n 4)))
         (* (/ 1262 105) (pow n 5)))
         (* (/ 73814 2835) (pow n 6))))
        ((= j 4)
         (-
         (- (* (/ 4279 630) (pow n 4)) (* (/ 332 35) (pow n 5)))
         (* (/ 399572 14175) (pow n 6))))
        ((= j 5)
         (- (* (/ 4174 315) (pow n 5)) (* (/ 144838 6237) (pow n 6))))
        ((= j 6)
         (* (/ 601676 22275) (pow n 6)))))

(define (A j)
  (cond ((= j 0)
         (+ 1
            (/ (pow n 2) 4)
            (/ (pow n 4) 64)))
        ((= j 1)
         (* -1 (/ 3 2) (- n (+
                             (/ (pow n 3) 8)
                             (/ (pow n 5) 64)))))
        ((= j 2)
         (* (/ 15 16) (- (pow n 2) (/ (pow n 4) 4))))
        ((= j 3)
         (* -1 (/ 35 48) (- (pow n 3) (* (/ 5 16) (pow n 5)))))
        ((= j 4)
         (* (/ 315 512) (pow n 4)))
        ((= j 5)
         (* -1 (/ 693 1280) (pow n 5)))))
(define (b j)
  (cond ((= j 1)
         (-
         (-
         (+
         (- (* (/ 1 2) n) (* (/ 2 3) (pow n 2)))
         (* (/ 37 96) (pow n 3)))
         (* (/ 1 360) (pow n 4)))
         (* (/ 81 512) (pow n 5))))
        ((= j 2)
         (+
         (-
         (+ (* (/ 1 48) (pow n 2)) (* (/ 1 15) (pow n 3)))
         (* (/ 437 1440) (pow n 4))
         (* (/ 46 105) (pow n 5)))))
        ((= j 3)
         (-
         (- (* (/ 17 480) (pow n 3)) (* (/ 37 840) (pow n 4)))
         (* (/ 209 4480) (pow n 5))))
        ((= j 4)
         (- (* (/ 4397 161280) (pow n 4)) (* (/ 11 504) (pow n 5))))
        ((= j 5)
         (* (/ 4583 161280) (pow n 5)))))

;有象無象の計算
(define Sphi (* (/ (* m0 alfa) (+ 1 n)) (+ (* (A 0) (/ la0 p))
                                           (sigma 1 5
                                                  (lambda (j)
                                                    (* (A j) (sin (* 2 j la0))))))))
(define A_ (* (/ (* m0 alfa) (+ 1 n)) (A 0)))

(define psi-0 (/ (+ x Sphi) A_))
(define nyu-0 (/ y A_))

(define psi (- psi-0
               (sigma 1 5
                      (lambda (j) (* (b j)
                                     (sin (* 2 j psi-0))
                                     (cosh (* 2 j nyu-0)))))))
(define nyu (- nyu-0
               (sigma 1 5
                      (lambda (j) (* (b j)
                                     (cos (* 2 j psi-0))
                                     (sinh (* 2 j nyu-0)))))))
(define theta (- 1
               (sigma 1 5
                      (lambda (j) (* 2
                                     j
                                     (b j)
                                     (cos (* 2 j psi-0))
                                     (cosh (* 2 j nyu-0)))))))         
(define tau (- 1
               (sigma 1 5
                      (lambda (j) (* 2
                                     j
                                     (b j)
                                     (sin (* 2 j psi-0))
                                     (sinh (* 2 j nyu-0)))))))
                                                                

(define kai (asin (/ (sin psi) (cosh nyu))))

;経度緯度
(define la (+ kai (* p (sigma 1 6 (lambda (j) (d j) (sin (* 2 j kai)))))))
(define lo (+ lo0 (atan (/ (sinh nyu) (cos psi)))))
