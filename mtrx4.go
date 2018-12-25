package main

import "fmt"

type mtrx4_t [16]float32

/* func proto

func mtrx4_idtt() (rt mtrx4_t)
func mtrx4_set(m [16]float32) (rt mtrx4_t)
func mtrx4_set_float(a00, a01, a02,
					 a10, a11, a12,
					 a20, a21, a22 float32) (rt mtrx4_t)
func mtrx4_set_euler(yaw, pitch, roll float32) (rt mtrx4_t)
func mtrx4_set_axisangl(ax vec3_t, phi float32) (rt mtrx4_t)
func mtrx4_show(m mtrx4_t)

*/

func mtrx4_idtt() (rt mtrx4_t) {
	var (
		i, j int32
		n    int32 = 4
	)

	for i = 0; i < n; i++ {
		for j = 0; j < n; j++ {
			if i == j {
				rt[id_rw(i, j, n)] = 1.0
			} else {
				rt[id_rw(i, j, n)] = 0.0
			}
		}
	}

	return rt
}

func mtrx4_set(m [16]float32) (rt mtrx3_t) {
	var (
		i, j int32
		n    int32 = 4
	)

	for i = 0; i < n; i++ {
		for j = 0; j < n; j++ {
			rt[id_rw(i, j, n)] = m[id_rw(i, j, n)]
		}
	}
	return rt
}

func mtrx4_set_float(a00, a01, a02, a03,
	a10, a11, a12, a13,
	a20, a21, a22, a23,
	a30, a31, a32, a33 float32) (rt mtrx4_t) {

	rt[0] = a00
	rt[1] = a01
	rt[2] = a02
	rt[3] = a03

	rt[4] = a10
	rt[5] = a11
	rt[6] = a12
	rt[7] = a13

	rt[8] = a20
	rt[9] = a21
	rt[10] = a22
	rt[11] = a23

	rt[12] = a30
	rt[13] = a31
	rt[14] = a32
	rt[15] = a33

	return rt
}

func mtrx4_set_euler(yaw, pitch, roll float32) (rt mtrx4_t) {
	var (
		cosy, siny, cosp, sinp, cosr, sinr float32
	)

	cosy = cosf(yaw)
	siny = sinf(yaw)
	cosp = cosf(pitch)
	sinp = sinf(pitch)
	cosr = cosf(roll)
	sinr = sinf(roll)

	rt[0] = cosy*cosr - siny*cosp*sinr
	rt[1] = -cosy*sinr - siny*cosp*cosr
	rt[2] = siny * sinp
	rt[3] = 0.0

	rt[4] = siny*cosr + cosy*cosp*sinr
	rt[5] = -siny*sinr + cosy*cosp*cosr
	rt[6] = -cosy * sinp
	rt[7] = 0.0

	rt[8] = sinp * sinr
	rt[9] = sinp * cosr
	rt[10] = cosp
	rt[11] = 0.0

	rt[12] = 0.0
	rt[13] = 0.0
	rt[14] = 0.0
	rt[15] = 1.0

	return rt
}

func mtrx4_set_axisangl(ax vec3_t, phi float32) (rt mtrx4_t) {
	var (
		cosphi, sinphi, vxvy, vxvz, vyvz, vx, vy, vz float32
	)

	cosphi = cosf(phi)
	sinphi = sinf(phi)
	vxvy = ax[_XC] * ax[_YC]
	vxvz = ax[_XC] * ax[_ZC]
	vyvz = ax[_YC] * ax[_ZC]
	vx = ax[_XC]
	vy = ax[_YC]
	vz = ax[_ZC]

	rt[0] = cosphi + (1.0-cosphi)*vx*vx
	rt[1] = (1.0-cosphi)*vxvy - sinphi*vz
	rt[2] = (1.0-cosphi)*vxvz + sinphi*vy
	rt[3] = 0.0

	rt[4] = (1.0-cosphi)*vxvy + sinphi*vz
	rt[5] = cosphi + (1.0-cosphi)*vy*vy
	rt[6] = (1.0-cosphi)*vyvz - sinphi*vz
	rt[7] = 0.0

	rt[8] = (1.0-cosphi)*vxvz - sinphi*vy
	rt[9] = (1.0-cosphi)*vyvz + sinphi*vx
	rt[10] = cosphi + (1.0-cosphi)*vz*vz
	rt[11] = 0.0

	rt[12] = 0.0
	rt[13] = 0.0
	rt[14] = 0.0
	rt[15] = 1.0

	return rt
}

func mtrx4_show(m mtrx4_t) {
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[0], m[1], m[2], m[3])
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[4], m[5], m[6], m[7])
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[8], m[9], m[10], m[11])
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[12], m[13], m[14], m[15])
}
