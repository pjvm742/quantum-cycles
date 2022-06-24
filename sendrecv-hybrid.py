#!/usr/bin/python

from dwave.system import LeapHybridSampler
import pandas
import sys

data = pandas.read_csv(sys.stdin, sep = '\t', header = None)
data = data.apply(pandas.to_numeric)
#print(data)

hybrid = LeapHybridSampler()
result = hybrid.sample_qubo(data)

#res1 = result.first.sample
#res = tuple(res1.values())
#print(res)

res = result.aggregate()
out = res.to_pandas_dataframe()
out.to_csv(sys.stdout, sep = '\t', header = False, index = False)
