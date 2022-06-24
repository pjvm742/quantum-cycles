package main

import (
	"math/rand"
	"math"
	"fmt"
	"os"
	"encoding/csv"
)

/* arc and vert are for the graph representation
 * that fits with Cui's data;
 * each arc is in one from-list and one to-list*/
type arc struct {
	start int
	end int
	cpty int
}
type vert struct {
	from []arc
	to []arc
}

type capmat [][]int
type wgtmat [][]int
/* capacitated weighted directed graph */
type cwdgraph struct {
	arcc capmat
	arcw wgtmat
}

func mkweights(a capmat) wgtmat {
	n := len(a)

	underw := make([]int, n*n)
	w := make([][]int, n)
	for i := 0; i < n; i++ {
		w[i] = underw[i*n : (i+1)*n]
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a[i][j] > 0 {
				w[i][j] = randwgt()
			}
		}
	}

	return w
}

func randwgt() int {
	options := [6]int{2,2,3,4,6,11}
	i := rand.Intn(6)
	return options[i]
}

func mkarcs(n int) capmat {
	undera := make([]int, n*n)
	a := make([][]int, n)
	for i := 0; i < n; i++ {
		a[i] = undera[i*n : (i+1)*n]
	}

	a[1][0] = randcap()
	a[2][0] = randcap()

	m := int(math.Sqrt(float64(n)))
	for i := 0; i < 3; i++ {
		for k := 0; k < m; k++ {
			j := rand.Intn(n)
			if i != j {
				a[i][j] = randcap()
			}
		}
	}

	for i := 3; i < n; i++ {
		for k := 0; k < 3; k++ {
			j := rand.Intn(n)
			if i != j {
				a[i][j] = randcap()
			}
		}
		j := rand.Intn(3)
		a[i][j] = randcap()
	}

	return a
}

func randcap() int {
	options := [6]int{1,1,2,3,4,6}
	i := rand.Intn(6)
	return options[i]
}

func simplifycui(v []vert) capmat {
	change := true
	n := len(v)
	var vtx vert

	for change {
		change = false
		for i := 0; i < n; i++ {
			vtx = v[i]
			af := vtx.from
			at := vtx.to
			if len(af) == 0 && len(at) == 0 {
				continue
			} else if len(af) == 0 {
				for _, a := range at {
					orig := a.start
					for p, b := range v[orig].from {
						if b.end == i {
							v[orig].from = append(v[orig].from[:p], v[orig].from[p+1:]...)
						}
					}
				}
				v[i].to = nil
				change = true
			} else if len(at) == 0 {
				for _, a := range af {
					dest := a.end
					for p, b := range v[dest].to {
						if b.start == i {
							v[dest].to = append(v[dest].to[:p], v[dest].to[p+1:]...)
						}
					}
				}
				v[i].from = nil
				change = true
			}
		}
	}

	new_indices := make([]int, n)
	var remaining []vert
	m := 0
	for i := 0; i < n; i++ {
		if len(v[i].to) == 0 || len(v[i].from) == 0 {
			continue
		}
		remaining = append(remaining, v[i])
		new_indices[i] = m
		m++
	}

	underarcs := make([]int, m*m)
	arcs := make([][]int, m)
	for i := 0; i < m; i++ {
		arcs[i] = underarcs[i*m : (i+1)*m]
	}
	
	for i := 0; i < m; i++ {
		for _, a := range remaining[i].from {
			j := new_indices[a.end]
			arcs[i][j] = a.cpty
		}
	}

	return arcs
}

func main() {
	args := os.Args
	nargs := len(args)
	var size int = 30
	var cui bool
	for a := 1; a < nargs; a++ {
		if args[a] == "-c" {
			cui = true
		} else {
			_, err := fmt.Sscanf(args[a], "%d", &size)
			if err != nil {
				fmt.Fprintln(os.Stderr, "Malformed arguments: remaining argument isn't an integer")
				os.Exit(1)
			}
			if size < 4 {
				fmt.Fprintln(os.Stderr, "Requested instance size is too small, need at least 4")
				os.Exit(1)
			}
		}
	}
	var arcs capmat
	if cui {
		arcs = readcui()
	} else {
		arcs = mkarcs(size)
	}
	weights := mkweights(arcs)
	printgraph(cwdgraph{arcs, weights})
}

func readcui() capmat {
	r := csv.NewReader(os.Stdin)
	r.Comma = '\t'

	var v []vert
	var ids []int
	
	var idi int
	var idj int
	var val float64
	var c int
	row, err := r.Read() // discard column names
	row, err = r.Read()
	for ; err == nil; row, err = r.Read() {
		fmt.Sscanf(row[0], "%d", &idi)
		fmt.Sscanf(row[1], "%d", &idj)
		fmt.Sscanf(row[2], "%g", &val)
		c = int(val)
		
		i := 0
		for ; i < len(v); i++ {
			if ids[i] == idi {
				break
			}
		}
		if i == len(v) {
			v = append(v, vert{})
			ids = append(ids, idi)
		}

		j := 0
		for ; j < len(v); j++ {
			if ids[j] == idj {
				break
			}
		}
		if j == len(v) {
			v = append(v, vert{})
			ids = append(ids, idj)
		}
		
		a := arc{i, j, c}
		v[i].from = append(v[i].from, a)
		v[j].to = append(v[j].to, a)
	}

	return simplifycui(v)
}

func printgraph(g cwdgraph) {
	w := csv.NewWriter(os.Stdout)
	w.Comma = '\t'

	arcs := g.arcc
	weights := g.arcw
	n := len(arcs)

	rep := make([][]string, 2*n)
	i := 0
	for ; i < n; i++ {
		row := make([]string, n)
		for j := 0; j < n; j++ {
			row[j] = fmt.Sprintf("%d", arcs[i][j])
		}
		rep[i] = row
	}
	for ; i < 2*n; i++ {
		row := make([]string, n)
		for j := 0; j < n; j++ {
			row[j] = fmt.Sprintf("%d", weights[i-n][j])
		}
		rep[i] = row
	}
	
	w.WriteAll(rep)
}
