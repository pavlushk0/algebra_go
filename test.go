package main

import (
	"fmt"
)

func main() {
	var (
		l, u     mtrx3_t
		mlu, tmp mtrx3_t
	)

	mlu = mtrx3Set([9]float32{
		10.0, -7.0, 0.0,
		-3.0, 6.0, 2.0,
		5.0, -1.0, 5.0})

	l, u = mtrx3LU(mlu)

	fmt.Println("L = ")
	mtrx3Show(l)

	fmt.Println("U = ")
	mtrx3Show(u)

	fmt.Println("UxL = ")
	tmp = mtrx3Mult(u, l)
	mtrx3Show(tmp)

	v1 := vec2Set(2.4, 4.1)
	v2 := vec2Set(4.4, 1.1)

	fmt.Println(vec2Dot(v1, v2))
	fmt.Println(vec2Dot(v1, v2))

}
