def coeffs(zeros, start = 0):
	degree = len(zeros) - start
	if degree == 0:
		return [0]
	if degree == 1:
		return [zeros[start], 1]

	total = [0 for i in range(degree + 1)]
	first = zeros[start]

	shifted = coeffs(zeros, start + 1)

	for i in range(degree):
		total[i] += first * shifted[i]
		total[i + 1] += shifted[i]

	return total