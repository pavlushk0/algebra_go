package main

/*
	Function prototypes

func id_rw(i int32, j int32, n int32) int32
func id_cw(i int32, j int32, n int32) int32
*/
const f_eps float32 = 0.0000001
const M_PI float32 = 3.14159265358979323846264338327950288419716939937510582097494459
const M_PI_2 float32 = M_PI * 0.5
const M_2_PI float32 = 2.0 / M_PI

func id_rw(i int32, j int32, n int32) int32 {
	return (i*n + j)
}

func id_cw(i int32, j int32, n int32) int32 {
	return (j*n + i)
}

func fabs(x float32) float32 {
	if x < 0.0 {
		return -x
	} else {
		return x
	}
}

func cosf(x float32) float32 {
	if x < 0.0 {
		x = -x
	}

	for M_PI < x {
		x -= M_2_PI
	}

	return 1.0 - (x*x/2.0)*(1.0-(x*x/12.0)*(1.0-(x*x/30.0)*(1.0-x*x/56.0)))
}

func sinf(x float32) float32 {
	return cosf(x - M_PI_2)
}
