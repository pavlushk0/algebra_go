package main

import (
	"fmt"
	"math"
)

type vec3_t [3]float32

const _XC int8 = 0
const _YC int8 = 1
const _ZC int8 = 2
const _WC int8 = 3

/*
	Function prototypes

func vec3_show(v vec3_t)
func vec3_set(x float32, y float32, z float32) (rt vec3_t)
func vec3_lenght(v vec3_t) float32
func vec3_normalize_self(v vec3_t)
func vec3_get_normalize(v vec3_t) (rt vec3_t)
func vec3_scale(v vec3_t, scale float32) vec3_t
func vec3_invert_self(v vec3_t)
func vec3_get_invert(v vec3_t) (rt vec3_t)
func vec3_dot(a vec3_t, b vec3_t) float32
func vec3_sum(a, b vec3_t) (rt vec3_t)
func vec3_sub(a, b vec3_t) (rt vec3_t)
func vec3_cross(a, b vec3_t) (rt vec3_t)

*/

func vec3_show(v vec3_t) {
	fmt.Printf("%5.2f %5.2f %5.2f\n", v[_XC], v[_YC], v[_ZC])
}

func vec3_set(x float32, y float32, z float32) (rt vec3_t) {
	rt[0] = x
	rt[1] = y
	rt[2] = z

	return rt
}

func vec3_lenght(v vec3_t) float32 {
	return float32(math.Sqrt(float64(
		v[_XC]*v[_XC] +
			v[_YC]*v[_YC] +
			v[_ZC]*v[_ZC])))

}

func vec3_normalize_self(v vec3_t) {
	var (
		len float32
	)

	len = vec3_lenght(v)

	if len != 0.0 {
		v[_ZC] = v[_ZC] / len
		v[_XC] = v[_XC] / len
		v[_YC] = v[_YC] / len
	}
}

func vec3_get_normalize(v vec3_t) (rt vec3_t) {
	rt = v

	vec3_normalize_self(rt)

	return rt
}

func vec3_scale(v vec3_t, scale float32) vec3_t {
	var rt vec3_t

	v[0] *= scale
	v[1] *= scale
	v[2] *= scale

	return rt
}

func vec3_invert_self(v vec3_t) {
	v[_XC] = -v[_XC]
	v[_YC] = -v[_YC]
	v[_ZC] = -v[_ZC]
}

func vec3_get_invert(v vec3_t) (rt vec3_t) {
	rt = v

	vec3_invert_self(rt)

	return rt
}

func vec3_dot(a vec3_t, b vec3_t) float32 {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

func vec3_sum(a, b vec3_t) (rt vec3_t) {
	rt[0] = a[0] + b[0]
	rt[1] = a[1] + b[1]
	rt[2] = a[2] + b[2]

	return rt
}

func vec3_sub(a, b vec3_t) (rt vec3_t) {
	rt[0] = a[0] - b[0]
	rt[1] = a[1] - b[1]
	rt[2] = a[2] - b[2]

	return rt
}

func vec3_cross(a, b vec3_t) (rt vec3_t) {
	rt[0] = a[_YC]*b[_ZC] - a[_ZC]*b[_YC]
	rt[0] = a[_ZC]*b[_XC] - a[_XC]*b[_ZC]
	rt[0] = a[_XC]*b[_YC] - a[_YC]*b[_XC]

	return rt
}
