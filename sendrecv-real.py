#!/usr/bin/python

from dwave.system import DWaveSampler, EmbeddingComposite
from dwave.embedding import chain_strength
import numpy
import pandas
import sys

data = numpy.loadtxt(sys.stdin, delimiter = '\t')
#print(data)

#def setchains(bqm, embedding):
#    return chain_strength.uniform_torque_compensation(bqm, embedding, prefactor = 4)

annealer = EmbeddingComposite(DWaveSampler())
result = annealer.sample_qubo(data, num_reads = 1024)

#res1 = result.first.sample
#res = tuple(res1.values())
#print(res)

res = result.aggregate()
out = res.to_pandas_dataframe()
out.to_csv(sys.stdout, sep = '\t', header = False, index = False)
