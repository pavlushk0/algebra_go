package main

import (
	"fmt"
	"math"
)

type mtrx3_t [9]float32

/*	multidimensional array mapping, array[i][j]
	row-wise (C, C++):
	(0	1)
	(2	3)

	column-wise (Fortran, Matlab):
	(0	2)
	(1	3)
*/

/*
	Function prototypes

func id_rw(i int32, j int32, n int32) int32
func id_cw(i int32, j int32, n int32) int32

func mtrx3_set(m [9]float32) (rt mtrx3_t)
func mtrx3_set_euler(yaw, pitch, roll float64) (rt mtrx3_t)
func mtrx3_set_axisangl(ax vec3_t, phi float64) (rt mtrx3_t)
func mtrx3_show(m mtrx3_t)
func mtrx3_get_idtt() (rt mtrx3_t)
func mtrx3_det(m mtrx3_t) float32
func mtrx3_mult(a, b mtrx3_t) (rt mtrx3_t)
func mtrx3_mult_vec3(m mtrx3_t, v vec3_t) (rt vec3_t)
func mtrx3_lu(m mtrx3_t) (lm, um mtrx3_t)
func mtrx3_ldlt(m mtrx3_t) (lm mtrx3_t, dv vec3_t)
func mtrx3_get_transpose(m mtrx3_t) (rt mtrx3_t)
func mtrx3_tranpose_self(m mtrx3_t)

*/

func mtrx3_set(m [9]float32) (rt mtrx3_t) {
	var (
		i, j int32
		n    int32 = 3
	)

	for i = 0; i < n; i++ {
		for j = 0; j < n; j++ {
			rt[id_rw(i, j, n)] = m[id_rw(i, j, n)]
		}
	}

	return rt
}

func mtrx3_set_euler(yaw, pitch, roll float64) (rt mtrx3_t) {
	var (
		cosy, siny, cosp, sinp, cosr, sinr float32
	)
	cosy = float32(math.Cos(yaw))
	siny = float32(math.Sin(yaw))
	cosp = float32(math.Cos(pitch))
	sinp = float32(math.Sin(pitch))
	cosr = float32(math.Cos(roll))
	sinr = float32(math.Sin(roll))

	rt[0] = cosy*cosr - siny*cosp*sinr
	rt[1] = -cosy*sinr - siny*cosp*cosr
	rt[2] = siny * sinp

	rt[3] = siny*cosr + cosy*cosp*sinr
	rt[4] = -siny*sinr + cosy*cosp*cosr
	rt[5] = -cosy * sinp

	rt[6] = sinp * sinr
	rt[7] = sinp * cosr
	rt[8] = cosp

	return rt
}

func mtrx3_set_axisangl(ax vec3_t, phi float64) (rt mtrx3_t) {
	var (
		cosphi, sinphi, vxvy, vxvz, vyvz, vx, vy, vz float32
	)

	cosphi = float32(math.Cos(phi))
	sinphi = float32(math.Sin(phi))
	vxvy = ax[_XC] * ax[_YC]
	vxvz = ax[_XC] * ax[_ZC]
	vyvz = ax[_YC] * ax[_ZC]
	vx = ax[_XC]
	vy = ax[_YC]
	vz = ax[_ZC]

	rt[0] = cosphi + (1.0-cosphi)*vx*vx
	rt[1] = (1.0-cosphi)*vxvy - sinphi*vz
	rt[2] = (1.0-cosphi)*vxvz + sinphi*vy

	rt[3] = (1.0-cosphi)*vxvy + sinphi*vz
	rt[4] = cosphi + (1.0-cosphi)*vy*vy
	rt[5] = (1.0-cosphi)*vyvz - sinphi*vz

	rt[6] = (1.0-cosphi)*vxvz - sinphi*vy
	rt[7] = (1.0-cosphi)*vyvz + sinphi*vx
	rt[8] = cosphi + (1.0-cosphi)*vz*vz

	return rt
}

func mtrx3_show(m mtrx3_t) {
	fmt.Printf("%5.2f %5.2f %5.2f\n", m[0], m[1], m[2])
	fmt.Printf("%5.2f %5.2f %5.2f\n", m[3], m[4], m[5])
	fmt.Printf("%5.2f %5.2f %5.2f\n", m[6], m[7], m[8])
}

func mtrx3_get_idtt() (rt mtrx3_t) {
	var (
		i, j int32
		n    int32 = 3
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

func mtrx3_det(m mtrx3_t) float32 {
	return m[0]*m[4]*m[8] +
		m[6]*m[1]*m[5] +
		m[2]*m[3]*m[7] -
		m[0]*m[7]*m[5] -
		m[8]*m[3]*m[1]
}

func mtrx3_mult(a, b mtrx3_t) (rt mtrx3_t) {
	var (
		i, j int32
		n    int32 = 3
	)

	for i = 0; i < n; i++ {
		for j = 0; j < n; j++ {
			rt[id_rw(i, j, n)] =
				a[id_rw(0, j, n)]*b[id_rw(i, 0, n)] +
					a[id_rw(1, j, n)]*b[id_rw(i, 1, n)] +
					a[id_rw(2, j, n)]*b[id_rw(i, 2, n)]
		}
	}

	return rt
}

func mtrx3_mult_vec3(m mtrx3_t, v vec3_t) (rt vec3_t) {
	rt[_XC] = m[0]*v[_XC] + m[1]*v[_YC] + m[2]*v[_ZC]
	rt[_YC] = m[3]*v[_XC] + m[4]*v[_YC] + m[5]*v[_ZC]
	rt[_ZC] = m[6]*v[_XC] + m[7]*v[_YC] + m[8]*v[_ZC]

	return rt
}

/*
	Где-то здесь ошибка, долго искал
	ничего не вышло и взял код из сети
*/
/*
func mtrx3_lu(m mtrx3_t) (l, u mtrx3_t) {
	var (
		i, j, k int32
		lm, um  mtrx3_t
		sum     float32
	)

	for j = 0; j < 3; j++ {
		um[id_rw(0, j, 3)] = m[id_rw(0, j, 3)]
	}

	for j = 0; j < 3; j++ {
		lm[id_rw(j, 0, 3)] = m[id_rw(j, 0, 3)] / um[id_rw(0, 0, 3)]
	}

	for i = 1; i < 3; i++ {
		for j = i; j < 3; j++ {
			sum = 0.0
			for k = 0; k < i; k++ {
				sum = sum + lm[id_rw(i, k, 3)]*um[id_rw(k, j, 3)]
			}
			um[id_rw(i, j, 3)] = m[id_rw(i, j, 3)] - sum
		}
	}

	for i = 1; i < 3; i++ {
		for j = i; j < 3; j++ {
			if i > j {
				lm[id_rw(j, i, 3)] = 0.0
			} else {
				sum = 0.0
				for k = 0; k < i; k++ {
					sum = sum + lm[id_rw(j, k, 3)]*um[id_rw(k, i, 3)]
				}
				lm[id_rw(j, i, 3)] = (1.0 / um[id_rw(i, i, 3)]) * (m[id_rw(j, i, 3)] - sum)
			}
		}
	}

	return lm, um
}
*/

/*
	Нижнетреугольная (L, lm) матрица имеет единицы по диагонали
*/
func mtrx3_lu(m mtrx3_t) (lm, um mtrx3_t) {
	var (
		i, j, k int32
		n       int32 = 3
		sum     float32
	)

	for i = 0; i < n; i++ {
		for k = i; k < n; k++ {
			sum = 0
			for j = 0; j < i; j++ {
				sum += (lm[id_rw(i, j, n)] * um[id_rw(j, k, n)])
			}
			um[id_rw(i, k, n)] = m[id_rw(i, k, n)] - sum
		}

		for k = i; k < n; k++ {
			if i == k {
				lm[id_rw(i, i, n)] = 1.0
			} else {
				sum = 0
				for j = 0; j < i; j++ {
					sum += lm[id_rw(k, j, n)] * um[id_rw(j, i, n)]
				}
				lm[id_rw(k, i, n)] = (m[id_rw(k, i, n)] - sum) / um[id_rw(i, i, n)]
			}
		}
	}

	return lm, um
}

func mtrx3_ldlt(m mtrx3_t) (lm mtrx3_t, dv vec3_t) {
	var (
		i, j, k int32
		n       int32 = 3
		sum     float32
	)
	for i = 0; i < n; i++ {
		for j = i; j < n; j++ {
			sum = m[id_rw(j, i, n)]
			for k = 0; k < i; k++ {
				sum = sum - lm[id_rw(i, k, n)]*dv[k]*lm[id_rw(j, k, n)]
				if i == j {
					if sum <= 0 {
						fmt.Println("A is not positive deﬁnite")
						return mtrx3_get_idtt(), vec3_set(0.0, 0.0, 0.0)
					}
					dv[i] = sum
					lm[id_rw(i, i, n)] = 1.0
				} else {
					lm[id_rw(j, i, n)] = sum / dv[i]
				}
			}
		}
	}

	return lm, dv
}

func mtrx3_get_transpose(m mtrx3_t) (rt mtrx3_t) {
	var (
		i, j int32
		n    int32 = 3
		tmp  float32
	)

	rt = m

	for i = 0; i < n; i++ {
		for j = 0; j < i; j++ {
			tmp = rt[id_rw(i, i, n)]
			rt[id_rw(i, j, n)] = rt[id_rw(j, i, n)]
			rt[id_rw(j, i, n)] = tmp
		}
	}

	return rt
}

func mtrx3_tranpose_self(m mtrx3_t) {
	var (
		i, j int32
		n    int32 = 3
		tmp  float32
	)

	for i = 0; i < n; i++ {
		for j = 0; j < i; j++ {
			tmp = m[id_rw(i, i, n)]
			m[id_rw(i, j, n)] = m[id_rw(j, i, n)]
			m[id_rw(j, i, n)] = tmp
		}
	}

}
