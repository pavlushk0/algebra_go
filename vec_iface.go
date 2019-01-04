package main

type vec_t interface {
	vec_at(i int32) (bool, float32)
}

func (v vec2_t) vec_at(i int32) (bool, float32) {
	if (i < 0) || (i > 1) {
		return false, 0.0
	} else {
		return true, v[i]
	}
}

func (v vec3_t) vec_at(i int32) (bool, float32) {
	if (i < 0) || (i > 2) {
		return false, 0.0
	} else {
		return true, v[i]
	}
}

func vec_dot(a, b vec_t) float32 {
	var (
		rt     float32 = 0.0
		tmp    bool    = true
		i      int32   = 0
		af, bf float32
	)

	for tmp {
		tmp, af = a.vec_at(i)
		_, bf = b.vec_at(i)
		rt = rt + af*bf
		i++
	}

	return rt
}
