package main

/*
	Function prototypes

func id_rw(i int32, j int32, n int32) int32
func id_cw(i int32, j int32, n int32) int32
*/

func id_rw(i int32, j int32, n int32) int32 {
	return (i*n + j)
}

func id_cw(i int32, j int32, n int32) int32 {
	return (j*n + i)
}
