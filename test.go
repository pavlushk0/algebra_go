package main

import (
	"fmt"
)

func main() {
	var (
		l, u     mtrx3_t
		mlu, tmp mtrx3_t
	)

	mlu = mtrx3_set([9]float32{
		10.0, -7.0, 0.0,
		-3.0, 6.0, 2.0,
		5.0, -1.0, 5.0})

	l, u = mtrx3_lu(mlu)

	fmt.Println("L = ")
	mtrx3_show(l)

	fmt.Println("U = ")
	mtrx3_show(u)

	fmt.Println("UxL = ")
	tmp = mtrx3_mult(u, l)
	mtrx3_show(tmp)

	v1 := vec2_set(2.4, 4.1)
	v2 := vec2_set(4.4, 1.1)

	fmt.Println(vec2_dot(v1, v2))
	fmt.Println(vec_dot(v1, v2))

}
