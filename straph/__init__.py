from straph.generators import erdos_renyi, barabasi_albert
from straph.parser import parser, sort_csv
from straph.paths import Path
from straph.stream import (StreamGraph,
                           read_stream_graph,
                           stream_graph_from_events_list,
                           read_stream_graph_from_json,
                           DFS_iterative)
from straph.utils import nx_degree, hist_plot

name = "straph"

# Let users know if they're missing any of our hard dependencies
# From panda's github
hard_dependencies = ("numpy", "matplotlib", "dateutil", "networkx", "pandas","dpkt","sortedcontainers")
missing_dependencies = []

for dependency in hard_dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(f"{dependency}: {e}")

if missing_dependencies:
    raise ImportError(
        "Unable to import required dependencies:\n" + "\n".join(missing_dependencies)
    )
    
del hard_dependencies, dependency, missing_dependencies
