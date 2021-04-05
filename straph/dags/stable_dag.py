# Copyright (C) 2017-2020 Léo Rannou - Sorbonne Université/LIP6 - Thales
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from collections import defaultdict
from joblib import Parallel, delayed

from straph import components as comp
from straph.dags.dag import Dag


class StableDag(Dag):
    def __init__(self,
                 id=None,
                 times=None,
                 c_nodes=None,
                 c_links=None,
                 id_comp_to_comp=None,
                 node_to_id_comp=None,
                 segmented_node_to_id_comp=None,
                 adj_list=None,
                 ):
        """
        A basic constructor for the condensation DAG

        :param c_nodes : A list of St CC nodes (each component node represent a SCC : a set of nodes,
         a begin time, an end time)
        :param c_links : A list of directed link (each link represent two connected SCC)
        """
        super().__init__(id, times, c_nodes, c_links, id_comp_to_comp, node_to_id_comp,
                         segmented_node_to_id_comp, adj_list)

    def core_number(self, n_jobs=-1):
        L = defaultdict(list)

        def para_cores(cmp):
            return cmp.core_number()

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cores)(cmp) for cmp in self.c_nodes if cmp.size() > 1)
        for l in r:
            for k, v in l.items():
                L[k] += v
        return L

    def k_core(self, k, n_jobs=-1):
        L = []

        def para_cores(cmp):
            return cmp.k_core(k)

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cores)(cmp) for cmp in self.c_nodes if cmp.size() > 1)
        for l in r:
            L += l
        return L

    def all_cliques(self, n_jobs=-1):
        L = defaultdict(list)

        def para_cliques(cmp):
            return cmp.all_cliques()

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cliques)(cmp) for cmp in self.c_nodes if cmp.size() > 1)
        for l in r:
            for k, v in l.items():
                L[k] += v
        return L

    def k_clique(self, k, n_jobs=-1):
        L = []

        def para_cliques(cmp):
            return cmp.k_clique(k)

        r = Parallel(n_jobs=n_jobs, mmap_mode='r+')(
            delayed(para_cliques)(cmp) for cmp in self.c_nodes if cmp.size() > 1)
        for l in r:
            L += l
        return L

    ################################
    #       FORMAT                 #
    ################################

    def cluster_to_object(self):
        new_cnodes = []
        for id_cc, cc in self.id_comp_to_comp.items():
            assert type(cc) == list
            new_cnodes = comp.StableConnectedComponent(id=id_cc, times=(cc[0][0], cc[0][1]),
                                                       nodes=set([c[2] for c in cc]))
        self.c_nodes = new_cnodes
        self.id_comp_to_comp = {cc.id: cc for cc in new_cnodes}
