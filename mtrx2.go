package main

import (
	"fmt"
	"math"
)

type mtrx2_t [4]float32

/*
	Function prototypes

func mtrx2_idtt() (rt mtrx2_t)
func mtrx2_set(m [9]float32) (rt mtrx2_t)
func mtrx2_set_float(a00, a01, a10, a11 float32) (rt mtrx2_t)
func mtrx2_set_rot(phi float32)
func mtrx2_show(m mtrx2_t)
func mtrx2_det(m mtrx2_t) (rt float32)
func mtrx2_det_lu(m mtrx2_t) (rt float32)
func mtrx2_mult(a, b mtrx2_t) (rt mtrx2_t)
func mtrx2_mult_vec3(m mtrx2_t, v vec3_t) (rt vec3_t)
func mtrx2_lu(m mtrx2_t) (lm, um mtrx2_t)
func mtrx2_ldlt(m mtrx2_t) (lm mtrx2_t, dv vec3_t)
func mtrx2_transpose(m mtrx2_t) (rt mtrx2_t)
func mtrx2_invert(m mtrx2_t) (rt mtrx2_t)
func mtrx2_solve_gauss(m mtrx2_t, v vec3_t) (rt vec3_t)
func mtrx2_insert_cmn(m mtrx2_t, v vec3_t, cmn int32) (rt mtrx2_t)
func mtrx2_solve_kramer(m mtrx2_t, v vec3_t) (rt vec3_t)

*/

func mtrx2_idtt() (rt mtrx2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j int32
	)

	for i = 0; i < mrange; i++ {
		for j = 0; j < mrange; j++ {
			if i == j {
				rt[id_rw(i, j, mrange)] = 1.0
			} else {
				rt[id_rw(i, j, mrange)] = 0.0
			}
		}
	}

	return rt
}

func mtrx2_set(m [4]float32) (rt mtrx2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j int32
	)

	for i = 0; i < mrange; i++ {
		for j = 0; j < mrange; j++ {
			rt[id_rw(i, j, mrange)] = m[id_rw(i, j, mrange)]
		}
	}

	return rt
}

func mtrx2_set_float(a00, a01, a10, a11 float32) (rt mtrx2_t) {
	rt[0] = a00
	rt[1] = a01
	rt[2] = a10
	rt[3] = a11

	return rt
}

func mtrx2_rot(phi float32) (rt mtrx2_t) {
	var (
		cosphi, sinphi float32
	)

	sinphi = sinf(phi)
	cosphi = cosf(phi)

	rt[0] = cosphi
	rt[1] = -sinphi
	rt[2] = -sinphi
	rt[3] = cosphi

	return rt
}

func mtrx2_show(m mtrx2_t) {
	fmt.Printf("%5.2f %5.2f\n", m[0], m[1])
	fmt.Printf("%5.2f %5.2f\n", m[2], m[3])
}

func mtrx2_det(m mtrx2_t) float32 {
	return m[0]*m[3] - m[1]*m[2]
}

func mtrx2_det_lu(m mtrx2_t) (rt float32) {
	const (
		mrange int32 = 2
	)

	var (
		i            int32
		l, u         mtrx2_t
		l_det, u_det float32
	)

	l, u = mtrx2_lu(m)

	l_det = l[0]
	u_det = u[0]

	for i = 1; i < mrange; i++ {
		l_det *= l[id_rw(i, i, mrange)]
		u_det *= u[id_rw(i, i, mrange)]
	}

	return l_det * u_det
}

func mtrx2_mult(a, b mtrx2_t) (rt mtrx2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j, k int32
		tmp     float32
	)

	for i = 0; i < mrange; i++ {
		for j = 0; j < mrange; j++ {
			tmp = 0.0
			for k = 0; k < mrange; k++ {
				tmp = tmp + a[id_rw(k, j, mrange)]*b[id_rw(i, k, mrange)]
			}
			rt[id_rw(i, j, mrange)] = tmp
		}
	}

	return rt
}

func mtrx2_mult_vec(m mtrx2_t, v vec2_t) (rt vec2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j int32
		tmp  float32
	)

	for i = 0; i < mrange; i++ {
		tmp = 0
		for j = 0; j < mrange; j++ {
			tmp = tmp + m[id_rw(i, j, mrange)]*v[j]
		}
		rt[i] = tmp
	}

	return rt
}

/*
	Нижнетреугольная (L, lm) матрица имеет единицы по диагонали
*/

func mtrx2_lu(m mtrx2_t) (lm, um mtrx2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j, k int32
		sum     float32
	)

	for i = 0; i < mrange; i++ {
		for k = i; k < mrange; k++ {
			sum = 0
			for j = 0; j < i; j++ {
				sum += (lm[id_rw(i, j, mrange)] * um[id_rw(j, k, mrange)])
			}
			um[id_rw(i, k, mrange)] = m[id_rw(i, k, mrange)] - sum
		}

		for k = i; k < mrange; k++ {
			if i == k {
				lm[id_rw(i, i, mrange)] = 1.0
			} else {
				sum = 0
				for j = 0; j < i; j++ {
					sum += lm[id_rw(k, j, mrange)] * um[id_rw(j, i, mrange)]
				}
				lm[id_rw(k, i, mrange)] = (m[id_rw(k, i, mrange)] - sum) / um[id_rw(i, i, mrange)]
			}
		}
	}

	return lm, um
}

func mtrx2_ldlt(m mtrx2_t) (lm mtrx2_t, dv vec2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j, k int32
		sum     float32
	)

	for i = 0; i < mrange; i++ {
		for j = i; j < mrange; j++ {
			sum = m[id_rw(j, i, mrange)]
			for k = 0; k < i; k++ {
				sum = sum - lm[id_rw(i, k, mrange)]*dv[k]*lm[id_rw(j, k, mrange)]
				if i == j {
					if sum <= 0 {
						fmt.Println("A is not positive deﬁnite")
						return mtrx2_idtt(), vec2_set(0.0, 0.0)
					}
					dv[i] = sum
					lm[id_rw(i, i, mrange)] = 1.0
				} else {
					lm[id_rw(j, i, mrange)] = sum / dv[i]
				}
			}
		}
	}

	return lm, dv
}

func mtrx2_transpose(m mtrx2_t) (rt mtrx2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j int32
		tmp  float32
	)

	rt = m

	for i = 0; i < mrange; i++ {
		for j = 0; j < i; j++ {
			tmp = rt[id_rw(i, i, mrange)]
			rt[id_rw(i, j, mrange)] = rt[id_rw(j, i, mrange)]
			rt[id_rw(j, i, mrange)] = tmp
		}
	}

	return rt
}

func mtrx2_invert(m mtrx2_t) mtrx2_t {
	var (
		det float32
	)

	det = mtrx2_det(m)

	if fabs(det) < f_eps {
		fmt.Println("mtrx_invert(): determinant is a zero!\n")
		return mtrx2_idtt()
	}

	return mtrx2_set_float(m[3], -m[1]/det, -m[2]/det, m[0]/det)
}

func mtrx2_solve_gauss(m mtrx2_t, v vec3_t) (rt vec3_t) {
	const (
		mrange int32 = 2
	)

	var (
		i, j, k int32
		t       float32
		a       [mrange][mrange + 1]float32
	)

	for i = 0; i < mrange; i++ { //было ++i
		for j = 0; j < mrange; j++ { //было ++j
			a[i][j] = m[id_rw(i, j, mrange)]
			a[i][mrange] = v[i]
		}
	}

	/* Pivotisation */
	for i = 0; i < mrange; i++ {
		for k = i + 1; k < mrange; k++ {
			if math.Abs(float64(a[i][i])) < math.Abs(float64(a[k][i])) {
				for j = 0; j <= mrange; j++ {
					t = a[i][j]
					a[i][j] = a[k][j]
					a[k][j] = t
				}
			}
		}
	}

	/* прямой ход */
	for k = 1; k < mrange; k++ {
		for j = k; j < mrange; j++ {
			t = a[j][k-1] / a[k-1][k-1]
			for i = 0; i < mrange+1; i++ {
				a[j][i] = a[j][i] - t*a[k-1][i]
			}
		}
	}

	/* обратный ход */
	for i = mrange - 1; i >= 0; i-- {
		rt[i] = a[i][mrange] / a[i][i]
		for j = mrange - 1; j > i; j-- {
			rt[i] = rt[i] - a[i][j]*rt[j]/a[i][i]
		}
	}

	return rt
}

func mtrx2_insert_cmn(m mtrx2_t, v vec3_t, cmn int32) (rt mtrx2_t) {
	const (
		mrange int32 = 2
	)

	var (
		i int32
		j int32 = 0
	)

	rt = m

	for i = cmn; i < mrange*mrange; i += mrange {
		rt[i] = v[j]
		j++
	}

	return rt
}

func mtrx2_solve_kramer(m mtrx2_t, v vec3_t) (rt vec3_t) {
	const (
		mrange int32 = 2
	)

	var (
		i       int32
		det     float32
		kr_mtrx mtrx2_t
	)

	det = mtrx2_det(m)

	if fabs(det) < f_eps {
		fmt.Println("mtrx_solve_kramer(): system has no solve\n")
		return vec3_set(0.0, 0.0, 0.0)
	}

	for i = 0; i < mrange; i++ {
		kr_mtrx = mtrx2_insert_cmn(m, v, i)
		rt[i] = mtrx2_det(kr_mtrx) / det
	}

	return rt
}
