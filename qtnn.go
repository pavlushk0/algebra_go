package main

import (
	"fmt"
	"math"
)

type qtnn_t [4]float32

/*
	Function prototypes

func qtnn_show(q qtnn_t)
func qtnn_lenght(q qtnn_t) float32
func qtnn_normalize(q qtnn_t) (rt qtnn_t)
func qtnn_inverse(q qtnn_t) (rt qtnn_t)
func qtnn_scale(q qtnn_t, scale float32) (rt qtnn_t)
func qtnn_sum(a, b qtnn_t) (rt qtnn_t)
func qtnn_sub(a, b qtnn_t) (rt qtnn_t)
func qtnn_dot(a, b qtnn_t) float32
func qtnn_mult(a, b qtnn_t) (rt qtnn_t)
func qtnn_mult_vec3(a qtnn_t, b vec3_t) (rt qtnn_t)
func qtnn_from_vec3(v vec3_t) (rt qtnn_t)
func qtnn_from_axisangl(a vec3_t, phi float64) (rt qtnn_t)
func qtnn_from_euler(yaw, pitch, roll float64) (rt qtnn_t)
func qtnn_to_vec3(q qtnn_t) (rt vec3_t)
func qtnn_transform_vec3(a qtnn_t, b vec3_t) (rt vec3_t)

*/

func qtnn_show(q qtnn_t) {
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", q[_XC], q[_YC], q[_ZC], q[_WC])
}

func qtnn_lenght(q qtnn_t) float32 {
	return float32(math.Sqrt(float64(
		q[_XC]*q[_XC] +
			q[_YC]*q[_YC] +
			q[_ZC]*q[_ZC] +
			q[_WC]*q[_WC])))
}

func qtnn_normalize(q qtnn_t) (rt qtnn_t) {
	var (
		len float32
	)

	len = qtnn_lenght(q)

	if len != 0.0 {
		rt[_WC] = q[_WC] / len
		rt[_XC] = q[_XC] / len
		rt[_YC] = q[_YC] / len
		rt[_ZC] = q[_ZC] / len
	}

	return rt
}

func qtnn_inverse(q qtnn_t) (rt qtnn_t) {
	rt[_WC] = q[_WC]
	rt[_XC] = -q[_XC]
	rt[_YC] = -q[_YC]
	rt[_ZC] = -q[_ZC]

	return rt
}

func qtnn_scale(q qtnn_t, scale float32) (rt qtnn_t) {
	rt[_WC] = q[_WC] * scale
	rt[_XC] = q[_XC] * scale
	rt[_YC] = q[_YC] * scale
	rt[_ZC] = q[_ZC] * scale
	return rt
}

func qtnn_sum(a, b qtnn_t) (rt qtnn_t) {
	rt[0] = a[0] + b[0]
	rt[1] = a[1] + b[1]
	rt[2] = a[2] + b[2]
	rt[3] = a[3] + b[3]

	return rt
}

func qtnn_sub(a, b qtnn_t) (rt qtnn_t) {
	rt[0] = a[0] - b[0]
	rt[1] = a[1] - b[1]
	rt[2] = a[2] - b[2]
	rt[3] = a[3] - b[3]

	return rt
}

func qtnn_dot(a, b qtnn_t) float32 {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
}

func qtnn_mult(a, b qtnn_t) (rt qtnn_t) {
	rt[_WC] = a[_WC]*b[_WC] - a[_XC]*b[_XC] - a[_YC]*b[_YC] - a[_ZC]*b[_ZC]
	rt[_XC] = a[_WC]*b[_XC] + a[_XC]*b[_WC] + a[_YC]*b[_ZC] - a[_ZC]*b[_YC]
	rt[_YC] = a[_WC]*b[_YC] - a[_XC]*b[_ZC] + a[_YC]*b[_WC] + a[_ZC]*b[_XC]
	rt[_ZC] = a[_WC]*b[_ZC] + a[_XC]*b[_YC] - a[_YC]*b[_XC] + a[_ZC]*b[_WC]

	return rt
}

/* function is broken */
func qtnn_mult_vec3(a qtnn_t, b vec3_t) (rt qtnn_t) {
	rt[_WC] = -a[_WC]*b[_XC] - a[_YC]*b[_YC] - a[_ZC]*b[_ZC]
	rt[_XC] = a[_WC]*b[_XC] + a[_YC]*b[_ZC] - a[_ZC]*b[_YC]
	rt[_YC] = a[_WC]*b[_YC] - a[_XC]*b[_ZC] + a[_ZC]*b[_XC]
	rt[_ZC] = a[_WC]*b[_ZC] + a[_XC]*b[_YC] - a[_YC]*b[_XC]

	return rt
}

func qtnn_from_vec3(v vec3_t) (rt qtnn_t) {
	rt[_XC] = v[_XC]
	rt[_YC] = v[_YC]
	rt[_ZC] = v[_ZC]
	rt[_WC] = 0.0

	return rt
}

func qtnn_from_axisangl(a vec3_t, phi float64) (rt qtnn_t) {
	var (
		sinhalfphi float32
	)

	sinhalfphi = float32(math.Sin(phi * 0.5))

	rt[_WC] = float32(math.Cos(phi * 0.5))
	rt[_XC] = a[_XC] * sinhalfphi
	rt[_YC] = a[_YC] * sinhalfphi
	rt[_ZC] = a[_ZC] * sinhalfphi

	return rt
}

func qtnn_from_euler(yaw, pitch, roll float64) (rt qtnn_t) {
	var (
		qyaw, qpitch, qroll qtnn_t
	)

	qyaw = qtnn_from_axisangl(vec3_set(1.0, 0.0, 0.0), yaw)
	qpitch = qtnn_from_axisangl(vec3_set(0.0, 1.0, 0.0), pitch)
	qroll = qtnn_from_axisangl(vec3_set(0.0, 0.0, 1.0), roll)

	rt = qtnn_mult(qyaw, qpitch)

	rt = qtnn_mult(rt, qroll)

	return rt
}

func qtnn_to_vec3(q qtnn_t) (rt vec3_t) {
	return vec3_set(q[_XC], q[_YC], q[_ZC])
}

func qtnn_transform_vec3(a qtnn_t, b vec3_t) (rt vec3_t) {
	var (
		vq, tmp qtnn_t
	)

	vq = qtnn_from_vec3(b)

	tmp = qtnn_mult(a, vq)
	tmp = qtnn_mult(tmp, qtnn_inverse(a))

	return qtnn_to_vec3(tmp)
}
