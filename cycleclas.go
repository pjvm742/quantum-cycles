package main

import (
	"fmt"
	"os"
	"strings"
	"encoding/csv"
	"bufio"
)

/* also used to represent a flow */
type capmat [][]int
type wgtmat [][]int
/* capacitated weighted directed graph */
type cwdgraph struct {
	arcc capmat
	arcw wgtmat
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

func isfeasible(flows capmat) bool {
	feas := true
	n := len(flows)
	for i := 0; i < n; i++ {
		sum1 := 0
		sum2 := 0
		for j := 0; j < n; j++ {
			sum1 += flows[i][j]
			sum2 += flows[j][i]
		}
		if sum1 != sum2 || flows[i][i] > 0 {
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
		k := len(path)
		for i := 0; i < k; i++ {
			if path[i] == u {
				start = i
			}
		}
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
		path = nil
	}
}

func vdis2adis(g cwdgraph) cwdgraph {
	a := g.arcc
	w := g.arcw
	n := len(a)
	
	undernewa := make([]int, 4*n*n)
	newa := make([][]int, 2*n)
	for i := 0; i < 2*n; i++ {
		newa[i] = undernewa[i*2*n : (i+1)*2*n]
	}
	
	underneww := make([]int, 4*n*n)
	neww := make([][]int, 2*n)
	for i := 0; i < 2*n; i++ {
		neww[i] = underneww[i*2*n : (i+1)*2*n]
	}

	for i := 0; i < n; i++ {
		newa[2*i][2*i+1] = 1
		for j := 0; j < n; j++ {
			if a[i][j] == 0 {
				continue
			}
			newa[2*i+1][2*j] = a[i][j]
			neww[2*i+1][2*j] = w[i][j]
		}
	}

	return cwdgraph{newa, neww}
}

func adis2vdis(adisg cwdgraph) cwdgraph {
	adisflow := adisg.arcc
	adisw := adisg.arcw
	n := len(adisflow) / 2

	underflow := make([]int, n*n)
	flow := make([][]int, n)
	for i := 0; i < n; i++ {
		flow[i] = underflow[i*n : (i+1)*n]
	}

	underw := make([]int, n*n)
	w := make([][]int, n)
	for i := 0; i < n; i++ {
		w[i] = underw[i*n : (i+1)*n]
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			flow[i][j] = adisflow[2*i+1][2*j]
			w[i][j] = adisw[2*i+1][2*j]
		}
	}

	return cwdgraph{flow, w}
}

func setcaps1(a capmat) {
	n := len(a)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a[i][j] > 0 {
				a[i][j] = 1
			}
		}
	}
}

func main() {
	args := os.Args
	nargs := len(args)
	var filename string
	var adis bool
	var vdis bool
	for a := 1; a < nargs; a++ {
		if args[a] == "-a" {
			adis = true
		} else if args[a] == "-v" {
			vdis = true
			adis = true
		} else {
			filename = args[a]
		}
	}
	graph := readgraph(filename)
	if adis {
		setcaps1(graph.arcc)
	}
	reducedg, oldinds := simplify(graph)
	if vdis {
		reducedg = vdis2adis(reducedg)
	}
	fmt.Printf("After pre-processing, the number of vertices is %d\n", len(reducedg.arcc))
	outputfile := strings.Replace(filename, ".graph.tsv", ".graph.dimacs", 1)
	writedimacs(outputfile, reducedg)
	fmt.Fprintln(os.Stderr, "Please press enter when the solution file is there")
	fmt.Scanf("\n")
	solfile := strings.Replace(filename, ".graph.tsv", ".sol.dimacs", 1)
	readdimacs(solfile, reducedg.arcc)
	if vdis {
		reducedg = adis2vdis(reducedg)
	}
	arcflows := reducedg.arcc
	weights := reducedg.arcw
	if !isfeasible(arcflows) {
		fmt.Fprintln(os.Stderr, "The solution is invalid")
		os.Exit(1)
	}
	fmt.Printf("Solution value: %d\n", solval(arcflows, weights))
	printdecomp(arcflows, oldinds)
}

func writedimacs(filename string, g cwdgraph) {
	a := g.arcc
	w := g.arcw
	n := len(a)
	
	f, err := os.Create(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	m := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a[i][j] > 0 {
				m++
			}
		}
	}
	fmt.Fprintf(f, "p min %d %d\n", n, m)

	for i := 0; i < n; i++ {
		isum := 0
		osum := 0
		for j := 0; j < n; j++ {
			osum += a[i][j]
			isum += a[j][i]
		}
		s := osum - isum
		if s != 0 {
			fmt.Fprintf(f, "n %d %d\n", i+1, s)
		}
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a[i][j] > 0 {
				fmt.Fprintf(f, "a %d %d 0 %d %d\n", i+1, j+1, a[i][j], w[i][j])
			}
		}
	}
}

func readdimacs(filename string, a capmat) {
	f, err := os.Open(filename)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	s := bufio.NewScanner(f)

	var i int
	var j int
	var flow int
	var line string
	for s.Scan() {
		line = s.Text()
		if line[0] != 'f' {
			continue
		}
		_, err := fmt.Sscanf(line, "f %d %d %d", &i, &j, &flow)
		i--
		j--
		if err != nil {
			fmt.Fprintln(os.Stderr, "Malformed flow entry in input")
			os.Exit(1)
		}
		a[i][j] -= flow
	}
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
	//fmt.Print(undercap[1])
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
