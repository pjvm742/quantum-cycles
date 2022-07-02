package main

import (
	"fmt"
	"os"
	"io"
	"strings"
	"encoding/csv"
	"math"
	"sort"
)

type qubo [][]float64
/* also used to represent a flow */
type capmat [][]int
type wgtmat [][]int
/* capacitated weighted directed graph */
type cwdgraph struct {
	arcc capmat
	arcw wgtmat
}
type translentry struct {
	start int
	end int
	bitval int
}
type transltable []translentry


func constructqubo(g cwdgraph, penmult float64, vdis bool, adis bool, rings []int, ringfactor float64) (qubo, transltable) {
	a := g.arcc
	n := len(a)

	varr := make([]translentry, n)
	voffsets := make([]int, n)
	vertsize := 0
	if vdis {
		for i := 0; i < n; i++ {
			varr[i] = translentry{i, i, 1}
			voffsets[i] = i
			vertsize++
		}
	} else {
		for i := 0; i < n; i++ {
			voffsets[i] = vertsize
			sum1 := 0
			sum2 := 0
			for j := 0; j < n; j++ {
				sum1 += a[i][j]
			}
			for j := 0; j < n; j++ {
				sum2 += a[j][i]
			}
			if sum1 < sum2 {
				varr[i] = translentry{i, i, sum1}
				vertsize += log2(sum1)
			} else {
				varr[i] = translentry{i, i, sum2}
				vertsize += log2(sum2)
			}
		}
	}

	underao := make([]int, n*n)
	aoffsets := make([][]int, n)
	for i := 0; i < n; i++ {
		aoffsets[i] = underao[i*n : (i+1)*n]
	}
	arcsize := 0
	m := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a[i][j] > 0 {
				if adis {
					aoffsets[i][j] = vertsize + m
					a[i][j] = 1
					arcsize++
					m++
				} else {
					aoffsets[i][j] = vertsize + arcsize
					arcsize += log2(a[i][j])
					m++
				}
			}
		}
	}
	aarr := make([]translentry, m)
	for i, k := 0, 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a[i][j] > 0 {
				aarr[k] = translentry{i, j, a[i][j]}
				k++
			}
		}
	}
	
	nvar := vertsize + arcsize

	vartable := make([]translentry, nvar)
	idx := 0
	for i := 0; i < n; i++ {
		val := varr[i].bitval
		nbits := log2(val)
		if nbits == 0 {
			continue
		}
		for p := 0; p < nbits -1; p++ {
			vartable[idx] = translentry{i, i, ( 1 << p )}
			idx++
		}
		vartable[idx] = translentry{i, i, val - ( 1 << (nbits -1)) +1 }
		idx++
	}
	for i := 0; i < m; i++ {
		orig := aarr[i].start
		dest := aarr[i].end
		val := aarr[i].bitval
		nbits := log2(val)
		if nbits == 0 {
			continue
		}
		for p := 0; p < nbits -1; p++ {
			vartable[idx] = translentry{orig, dest, ( 1 << p )}
			idx++
		}
		vartable[idx] = translentry{orig, dest, val - ( 1 << (nbits -1)) +1 }
		idx++
	}
	if idx < nvar {
		panic("Not enough variables in the translation table")
	}

	underQUBO := make([]float64, nvar*nvar)
 	/* first vertsize indices are for the vertices,
	 * the rest for the arcs */
	QUBO := make([][]float64, nvar)
	for i := 0; i < nvar; i++ {
		QUBO[i] = underQUBO[i*nvar : (i+1)*nvar]
	}

	/* first the base objective */
	aw := g.arcw
	for k := vertsize; k < nvar; k++ {
		arc := vartable[k]
		i := arc.start
		j := arc.end
		bv := arc.bitval
		QUBO[k][k] -= float64(aw[i][j] * bv)
	}

	/* now the penalties */
	for i := 0; i < n; i++ {
		localmult := penmult * math.Pow(ringfactor, float64(rings[i]))
		var indices []int
		var values []int
		for k := voffsets[i]; k < len(vartable); k++ {
			vbit := vartable[k]
			if vbit.start != i || vbit.end != i {
				break
			}
			indices = append(indices, k)
			values = append(values, vbit.bitval)
		}
		for j := 0; j < n; j++ {
			if a[i][j] == 0 {
				continue
			}
			for k := aoffsets[i][j]; k < len(vartable); k++ {
				abit := vartable[k]
				if abit.start != i || abit.end != j {
					break
				}
				indices = append(indices, k)
				values = append(values, - abit.bitval)
			}
		}
		addquadpen(QUBO, indices, values, localmult)
	}
	for j := 0; j < n; j++ {
		localmult := penmult * math.Pow(ringfactor, float64(rings[j]))
		var indices []int
		var values []int
		for k := voffsets[j]; k < len(vartable); k++ {
			vbit := vartable[k]
			if vbit.start != j || vbit.end != j {
				break
			}
			indices = append(indices, k)
			values = append(values, vbit.bitval)
		}
		for i := 0; i < n; i++ {
			if a[i][j] == 0 {
				continue
			}
			for k := aoffsets[i][j]; k < len(vartable); k++ {
				abit := vartable[k]
				if abit.start != i || abit.end != j {
					break
				}
				indices = append(indices, k)
				values = append(values, - abit.bitval)
			}
		}
		addquadpen(QUBO, indices, values, localmult)
	}

	return QUBO, vartable
}

func log2(k int) int {
	exp := 0
	for k >= ( 1 << exp ) {
		exp++
	}
	return exp
}

/* modifies qubomatrix */
func addquadpen(qubomatrix qubo, indices []int, values []int, mult float64) {
	k := len(indices)
	for p := 0; p < k; p++ {
		for q := 0; q < k; q++ {
			i := indices[p]
			j := indices[q]
			val1 := values[p]
			val2 := values[q]
			qubomatrix[i][j] += mult * float64(val1 * val2)
		}
	}
}

func mkrings(arcs capmat) ([]int, int) {
	n := len(arcs)
	rings := make([]int, n)
	rings[0] = 1

	r := 1
	full := false
	for !full {
		done := true
		for i:= 0; i < n; i++ {
			if rings[i] == r {
				for j := 0; j < n; j++ {
					if rings[j] == 0 && ( arcs[i][j] > 0 || arcs[j][i] > 0 ) {
						rings[j] = r + 1
						done = false
					}
				}
			}
		}

		if done {
			for i := 0; i < n; i++ {
				if rings[i] == 0 {
					rings[i] = 1
					break
				}
			}
			r = 1
		} else {
			r++
		}

		full = true
		for i := 0; i < n; i++ {
			if rings[i] == 0 {
				full = false
			}
		}
	}

	max := 0
	for i := 0; i < n; i++ {
		if rings[i] > max {
			max = rings[i]
		}
	}

	return rings, max
}

func simplify(g cwdgraph) (cwdgraph, []int) {
	a := g.arcc
	w := g.arcw
	n := len(a)

	vertcaps := make([]int, n)
	
	for i := 0; i < n; i++ {
		isum := 0
		osum := 0
		for j := 0; j < n; j++ {
			osum += a[i][j]
			isum += a[j][i]
		}
		if osum > isum {
			vertcaps[i] = isum
		} else {
			vertcaps[i] = osum
		}
	}

	change := true
	for change {
		change = false
		for i := 0; i < n; i++ {
			c := vertcaps[i]
			for j := 0; j < n; j++ {
				if a[i][j] > c {
					a[i][j] = c
				}
				if a[j][i] > c {
					a[j][i] = c
				}
			}
		}
		for i := 0; i < n; i++ {
			isum := 0
			osum := 0
			for j := 0; j < n; j++ {
				osum += a[i][j]
				isum += a[j][i]
			}
			if osum < vertcaps[i] {
				vertcaps[i] = osum
				change = true
			}
			if isum < vertcaps[i] {
				vertcaps[i] = isum
				change = true
			}
		}
	}

	nempty := 0
	for i := 0; i < n; i++ {
		if vertcaps[i] == 0 {
			nempty++
		}
	}
	if nempty == 0 {
		indices := make([]int, n)
		for i := 0; i < n; i++ {
			indices[i] = i
		}
		return g, indices
	}

	m := n - nempty
	old_indices := make([]int, m)
	for i, p := 0, 0; i < n; i++ {
		if vertcaps[i] > 0 {
			old_indices[p] = i
			p++
		}
	}

	undernewa := make([]int, m*m)
	underneww := make([]int, m*m)
	newa := make([][]int, m)
	neww := make([][]int, m)
	for i := 0; i < m; i++ {
		newa[i] = undernewa[i*m : (i+1)*m]
		neww[i] = underneww[i*m : (i+1)*m]
	}
	for p := 0; p < m; p++ {
		for q := 0; q < m; q++ {
			i, j := old_indices[p], old_indices[q]
			newa[p][q] = a[i][j]
			neww[p][q] = w[i][j]
		}
	}
	return cwdgraph{newa, neww}, old_indices
}

func adjusted_avg(vals []int) float64 {
	k := len(vals)
	sort.Ints(vals)
	v := make([]float64, k)
	for i := 0; i < k; i++ {
		v[i] = float64(vals[i])
	}
	var kicked = true
	var avg float64
	for kicked {
		kicked = false
		sum := float64(0)
		for i := 0; i < k; i++ {
			sum += v[i]
		}
		avg = sum / float64(k)
		varsum := float64(0)
		for i := 0; i < k; i++ {
			d := v[i] - avg
			varsum += d*d
		}
		variance := varsum / float64(k)
		stddev := math.Sqrt(variance)

		if v[k-1] > avg + 3*stddev {
			v = v[:k-1]
			k--
			kicked = true
		} else if v[0] < avg - 3*stddev {
			v = v[1:]
			k--
			kicked = true
		}
	}
	return avg
}

func solval(flows capmat, weights wgtmat) int {
	n := len(flows)
	value := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if flows[i][j] > 0 {
				value += flows[i][j] * weights[i][j]
			}
		}
	}
	return value
}

func isfeasible(vertices []int, arcs capmat) bool {
	feas := true
	n := len(vertices)
	for i := 0; i < n; i++ {
		flow := vertices[i]
		sum1 := 0
		sum2 := 0
		for j := 0; j < n; j++ {
			sum1 += arcs[i][j]
			sum2 += arcs[j][i]
		}
		if sum1 != flow || sum2 != flow {
			feas = false
		}
	}
	return feas
}

/* modifies flows */
func printdecomp(flows capmat, old_indices []int) {
	n := len(flows)
	marked := make([]bool, n)
	var path []int
	
	for {
		empty := true
		u := 0
initsearch:
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				if flows[i][j] > 0 {
					empty = false
					u = i
					break initsearch
				}
			}
		}
		if empty {
			break
		}

		for !marked[u] {
			marked[u] = true
			path = append(path, u)
			for j := 0; j < n; j++ {
				if flows[u][j] > 0 {
					u = j
					break
				}
			}
		}

		start := 0
		for i := range path {
			if path[i] == u {
				start = i
			}
		}
		k := len(path)
		min := flows[ path[k-1] ][u]
		for p := start; p < k-1; p++ {
			f := flows[ path[p] ][ path[p+1] ]
			if f < min {
				min = f
			}
		}

		flows[ path[k-1] ][u] -= min
		for p := start; p < k-1; p++ {
			flows[ path[p] ][ path[p+1] ] -= min
		}

		fmt.Printf("Flow of %d along ", min)
		for p := start; p < k; p++ {
			fmt.Printf("%d -> ", old_indices[ path[p] ])
		}
		fmt.Printf("%d\n", old_indices[u])

		for i := 0; i < n; i++ {
			marked[i] = false
		}
		path = path[0:0]
	}
}

func backtranslate(sol []bool, tlt transltable, n int) ([]int, capmat) {
	vertices := make([]int, n)
	underflow := make([]int, n*n)
	flow := make([][]int, n)
	for i := 0; i < n; i++ {
		flow[i] = underflow[i*n : (i+1)*n]
	}

	nvar := len(tlt)
	for p := 0; p < nvar; p++ {
		if sol[p] {
			vrbl := tlt[p]
			i := vrbl.start
			j := vrbl.end
			if i == j {
				vertices[i] += vrbl.bitval
			} else {
				flow[i][j] += vrbl.bitval
			}
		}
	}
	return vertices, flow
}

func main() {
	args := os.Args
	nargs := len(args)
	var filename string
	var penmult float64 = 1
	var absmult bool
	var expectmult bool
	var ringf float64 = 1
	var expectringf bool
	var adis bool
	var vdis bool
	for a := 1; a < nargs; a++ {
		if expectmult {
			_, err := fmt.Sscanf(args[a], "%f", &penmult)
			if err != nil {
				fmt.Fprintln(os.Stderr, "Malformed arguments: argument after -m isn't a number")
				os.Exit(1)
			}
			expectmult = false
		} else if expectringf {
			_, err := fmt.Sscanf(args[a], "%f", &ringf)
			if err != nil {
				fmt.Fprintln(os.Stderr, "Malformed arguments: argument after -r isn't a number")
				os.Exit(1)
			}
			expectringf = false
		} else if args[a] == "-a" {
			adis = true
		} else if args[a] == "-v" {
			vdis = true
			adis = true
		} else if args[a] == "-m" {
			expectmult = true
		} else if args[a] == "-M" {
			absmult = true
			expectmult = true
		} else if args[a] == "-r" {
			expectringf = true
		} else {
			filename = args[a]
		}
	}
	graph := readgraph(filename)
	reducedg, oldinds := simplify(graph)
	fmt.Printf("After pre-processing, the number of vertices is %d\n", len(reducedg.arcc))
	if !absmult {
		var weights []int
		n := len(reducedg.arcc)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				wgt := reducedg.arcw[i][j]
				if wgt > 0 {
					weights = append(weights, wgt)
				}
			}
		}
		avg := adjusted_avg(weights)
		fmt.Fprintf(os.Stderr, "Adjusted average  of the weights (scale for the penalties): %.2g\n", avg)
		penmult *= avg
	}
	rings, max := mkrings(reducedg.arcc)
	ringfactor := math.Pow(ringf, 1/float64(max))
	penmult *= 1/ringfactor
	qubomatrix, tlt := constructqubo(reducedg, penmult, vdis, adis, rings, ringfactor)
	outputfile := strings.Replace(filename, ".graph.tsv", ".qubo.tsv", 1)
	writeQUBO(outputfile, qubomatrix)
	fmt.Fprintln(os.Stderr, "Please press enter when the solution file is there")
	fmt.Scanf("\n")
	solfile := strings.Replace(filename, ".graph.tsv", ".sol.tsv", 1)
	processsolutions(solfile, tlt, oldinds, reducedg.arcw)
}

func processsolutions(filename string, tlt transltable, oldinds []int, weights wgtmat) {
	n := len(weights)
	nvar := len(tlt)
	
	f, err := os.Open(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	r := csv.NewReader(f)
	r.Comma = '\t'

	s := 0
	var firstsol []bool
	row, err := r.Read()
	for first := true; err == nil; row, err = r.Read() {
		sol := make([]bool, nvar)
		for i := 0; i < nvar; i++ {
			if row[i] == "0" {
				sol[i] = false
			} else if row[i] == "1" {
				sol[i] = true
			} else {
				fmt.Fprintf(os.Stderr, "Unexpected value in solution file (not 1 or 0) at index %d\n", i)
				os.Exit(1)
			}
		}
		if first {
			firstsol = sol
			first = false
		}
		vertflows, arcflows := backtranslate(sol, tlt, n)
		if isfeasible(vertflows, arcflows) {
			fmt.Fprintf(os.Stderr, "The %d-th solution is feasible\n", s)
			fmt.Printf("Solution value: %d\n", solval(arcflows, weights))
			printdecomp(arcflows, oldinds)
			os.Exit(0)
		}
		s++
	}
	if err != io.EOF {
		fmt.Fprintln(os.Stderr, "Something went wrong trying to read the solution file")
		os.Exit(1)
	}
	fmt.Fprintln(os.Stderr, "None of the solutions are feasible")
	fmt.Fprintln(os.Stderr, "Breaks in the first solution:")
	vertflows, arcflows := backtranslate(firstsol, tlt, n)
	showbreaks(vertflows, arcflows, oldinds)
}

func readgraph(filename string) cwdgraph {
	f, err := os.Open(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	r := csv.NewReader(f)
	r.Comma = '\t'

	firstrow, err := r.Read()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	n := len(firstrow)

	undercap := make([]int, n*n)
	capacities := make([][]int, n)
	pos := 0

	for j := 0; j < n; j++ {
		fmt.Sscanf(firstrow[j], "%d", &undercap[pos])
		pos++
	}
	capacities[0] = undercap[0:n]

	for i := 1; i < n; i++ {
		row, err := r.Read()
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		for j := 0; j < n; j++ {
			fmt.Sscanf(row[j], "%d", &undercap[pos])
			pos++
		}
		capacities[i] = undercap[i*n : (i+1)*n]
	}

	for i := 0; i < n; i++ {
		capacities[i][i] = 0
	}

	underwgt := make([]int, n*n)
	weights := make([][]int, n)
	pos = 0
	
	for i := 0; i < n; i++ {
		row, err := r.Read()
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		for j := 0; j < n; j++ {
			fmt.Sscanf(row[j], "%d", &underwgt[pos])
			pos++
		}
		weights[i] = underwgt[i*n : (i+1)*n]
	}

	g := cwdgraph{capacities, weights}
	return g
}

func writeQUBO(filename string, problem qubo) {
	f, err := os.Create(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	w := csv.NewWriter(f)
	w.Comma = '\t'

	n := len(problem)
	rep := make([][]string, n)
	for i := 0; i < n; i++ {
		row := make([]string, n)
		for j := 0; j < n; j++ {
			row[j] = fmt.Sprintf("%.6g", problem[i][j])
		}
		rep[i] = row
	}
	w.WriteAll(rep)
}

func showmat(mat [][]int) {
	n := len(mat)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			fmt.Printf("%d  ", mat[i][j])
		}
		fmt.Printf("\n")
	}
}

func showbreaks(vertices []int, arcs capmat, oldinds []int) {
	n := len(vertices)
	for i := 0; i < n; i++ {
		flow := vertices[i]
		sum1 := 0
		sum2 := 0
		for j := 0; j < n; j++ {
			sum1 += arcs[i][j]
			sum2 += arcs[j][i]
		}
		if sum1 != flow || sum2 != flow {
			fmt.Fprintf(os.Stderr, "Break at vertex %d; %d in simplified graph\n", oldinds[i], i)
			fmt.Fprintf(os.Stderr, "In: %d\n", sum2)
			fmt.Fprintf(os.Stderr, "Trough: %d\n", flow)
			fmt.Fprintf(os.Stderr, "Out: %d\n", sum1)
		}
	}
}
