# Use this space to brainstorm CLIPPy based graph interface

## Find the furthest 3 vertices from a seed in an unweighted graph
```python

  from clippy import Clippy
  c = Clippy(cmd_prefix="", loglevel=0)

  # Place to store the maptrix
  mymaptrix = "/dev/shm/foo"

  c.maptrix.create(mymaptrix, from_csv="/p/lustre1/pearce7/big_edge_list.csv")

  # Place to store vertex data to hold level information
  vldata = "/dev/shm/foo/foo_levels"
  
  c.algorithms.bfs(mymaptrix, source=0, out_level=vldata)

  c.vertex_data.topk(vldata, k=3, ascending=False)
 
```