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
}
