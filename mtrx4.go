package main

import "fmt"

type mtrx4_t [16]float32

/*
	Function prototypes

func mtrx4_show(m mtrx4_t)
func mtrx4_set(m [16]float32) (rt mtrx3_t)
func mtrx4_get_idtt() (rt mtrx4_t)

*/

func mtrx4_show(m mtrx4_t) {
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[0], m[1], m[2], m[3])
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[4], m[5], m[6], m[7])
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[8], m[9], m[10], m[11])
	fmt.Printf("%5.2f %5.2f %5.2f %5.2f\n", m[12], m[13], m[14], m[15])
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

func mtrx4_get_idtt() (rt mtrx4_t) {
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
