import matplotlib.pyplot as plt

from straph import stream as sg

plt.rcParams["figure.figsize"] = (15, 10)


class path:
    def __init__(self,
                 times=None,
                 links=None,
                 ):
        '''
        A basic constructor for a path object
        :param times : A list of times corresponding to the links (first time = beginning ; last time = ending)
        :param links : A list of links composing the path. (first node = source ; last node = destination)
                '''
        self.times = times
        self.links = links

    def add_link(self, l, t):
        self.times.append(t)
        self.links.append(l)

    def length(self):
        return len(self.links)

    def duration(self):
        return self.times[-1] - self.times[0]

    def plot(self, S, color="#18036f",
             markersize=10, dag=False):
        '''
        Draw a path on the Stream Graph S
        :param S: A Stream Graph
        :return:
        '''

        if dag:

            dag = S.condensation_dag()
            dag.plot(node_to_label=S.node_to_label)
        else:
            S.plot()
        # Plot Source
        id_source = S.nodes.index(self.links[0][0])
        plt.plot([self.times[0]], [id_source], color=color,
                 marker='o', alpha=0.8, markersize=markersize)
        # Plot Destination
        id_destination = S.nodes.index(self.links[-1][1])
        plt.plot([self.times[-1]], [id_destination], color=color,
                 marker='o', alpha=0.8, markersize=markersize)
        # Plot Path
        for i in range(self.length()):
            l = self.links[i]
            t = self.times[i]
            id1 = S.nodes.index(l[0])
            id2 = S.nodes.index(l[1])
            idmax = max(id1, id2)
            idmin = min(id1, id2)
            plt.vlines(t, ymin=idmin, ymax=idmax, linewidth=4, alpha=0.8, color=color)
            if i != self.length() - 1:
                plt.hlines(id2, xmin=t, xmax=self.times[i + 1],
                           linewidth=4, alpha=0.8, color=color)
                # Plot marker
                if t != self.times[i + 1]:
                    plt.plot([t], [id2], color=color,
                             marker='>', alpha=0.8, markersize=markersize)
            if i != 0 and (t,id1) != (self.times[0],id_source) != (self.times[-1],id_destination):
                # Plot marker
                if id1 == idmin:
                    plt.plot([t], [id1], color=color,
                             marker='^', alpha=0.8, markersize=markersize)
                else:
                    plt.plot([t], [id1], color=color,
                             marker='v', alpha=0.8, markersize=markersize)
        plt.tight_layout()

    def check_coherence(self, S):
        for i in range(self.length()):
            l = self.links[i]
            l_ = (self.links[i][1], self.links[i][0])  # Inverse the order of the link
            if l not in S.links and l_ not in S.links:
                raise ValueError("Link : " + str(l) + " does not exists in the Stream Graph !")
            else:
                t = self.times[i]
                if l in S.links:
                    id_link = S.links.index(l)
                else:
                    id_link = S.links.index(l_)
                is_present = False
                for lt0, lt1 in zip(S.link_presence[id_link][::2], S.link_presence[id_link][1::2]):
                    if lt0 <= t <= lt1:
                        is_present = True
                if not is_present:
                    raise ValueError("Link : " + str(l) + " does not exists at time " + str(t) + " !")
        print("Check Path Coherence ok !")
        return
