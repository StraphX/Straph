import matplotlib.cm as cmx
import matplotlib.colors as colors
import pstats,os,sys
import cProfile
import subprocess
import pathlib,re


import networkx as nx
import matplotlib.pyplot as plt

def get_cmap(N,cmap = 'hsv'):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color.'''
    color_norm = colors.Normalize(vmin=0, vmax=N - 1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=cmap)

    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)

    return map_index_to_rgb_color

def profile_shit(sf,snakeviz=False):
    # print(" VMPROF PROFILING")
    # dat_file = open("/home/leo/Documents/stats_profile.dat",'r+')
    # vmprof.enable(dat_file.fileno(),period=0.00099)
    # f()
    # vmprof.disable()
    # vmprof.read_profile(dat_file)

    #
    # print(" # YAPPI PROFILING #")
    # yappi.start()
    # f()
    # yappi.get_func_stats().print_all()
    # yappi.get_thread_stats().print_all()
    cProfile.run(str(sf),
                 'profiling_stats')
    print("# CProfile PROFILING #")
    p = pstats.Stats("profiling_stats")
    p.sort_stats('cumulative').print_stats(60)
    if snakeviz:
        subprocess.run([os.path.join(os.path.dirname(sys.executable),"snakeviz"),"profiling_stats"])


def describe_to_latex_landscape_old(df, latex_output_path, title=None):
    '''
    Take a pandas DataFrame and output a latex table.
    Requirements : latex, pdflatex
    :param df: A Pandas DataFrame (of a reasonable size, nb of columns).
    :param latex_output_path: Path to write the corresponding LaTeX table.
    :param title: Title of the output LaTeX table.
    :return:
    '''
    pathlib.Path(latex_output_path).mkdir(parents=True, exist_ok=True)
    if title:
        path_tex = latex_output_path + '/_describe_' + title + '.tex'
    else:
        path_tex = latex_output_path + '/_describe.tex'
    n_columns = len(df.columns)
    df_describe = df.describe(percentiles=[0.1, 0.50, 0.6, 0.7, 0.8, 0.90, 0.99, 0.999, 0.9999])
    with open(path_tex, 'w') as output:
        output.writelines(['\documentclass{article} \n',
                           '\\usepackage[landscape]{geometry} \n',
                           '\\begin{document} \n',
                           '\\newgeometry{margin=1.5cm} \n',
                           ])
        output.write('\\begin{table} \n')
        output.write('\centering \n')
        output.write('\t \\begin{tabular}[t]{' + '|c' * (n_columns + 1) + '|} \n')
        output.write('\t \t \hline \n')
        output.write('\t \t')
        output.writelines([' & ' + re.sub('_', '\_', c) for c in df.columns])
        output.write(' \\\\ \n')
        output.write('\t \t \hline \n')
        for i in df_describe.index:
            output.write('\t \t' + re.sub('%', '\%', i) + ' & ')
            output.writelines([str(round(df_describe.loc[i, c], 6)) + ' & '
                               if c != df_describe.columns[-1]
                               else str(df_describe.loc[i, c])
                               for c in df_describe.columns])
            output.write(' \\\\ \n')
            output.write('\t \t \hline \n')
        output.write('\t \end{tabular} \n')

        if title:
            output.write('\caption{Description Of ' + title + ' }')
        else:
            output.write('\caption{Description}')
        output.write('\end{table} \n')
        output.write('\n')
        output.write('\end{document}')
    path_pdf_directory = latex_output_path
    subprocess.run(["pdflatex", "-output-directory", path_pdf_directory, path_tex])


def dataframe_to_latex_landscape(df, latex_output_path, title=None,
                                 highlight_best_result = False):
    '''
    Take a pandas DataFrame and output a latex table.
    Requirements : latex, pdflatex
    :param df: A Pandas DataFrame (of a reasonable size, nb of columns).
    :param latex_output_path: Path to write the corresponding LaTeX table.
    :param title: Title of the output LaTeX table.
    :return:
    '''
    pathlib.Path(latex_output_path).mkdir(parents=True, exist_ok=True)
    if title:
        path_tex = latex_output_path + '/_describe_' + title + '.tex'
    else:
        path_tex = latex_output_path + '/_describe.tex'
    n_columns = len(df.columns)

    with open(path_tex, 'w') as output:
        output.writelines(['\documentclass{article} \n',
                           '\\usepackage[landscape]{geometry} \n',
                           '\\usepackage{tabularx} \n'
                           '\\begin{document} \n',
                           '\\newgeometry{margin=1.5cm} \n',
                           ])
        output.write('\\begin{figure*} \n')
        output.write('\\centering \n')
        output.write('\t \\begin{tabularx}{\\textwidth}{' + '|c'+'|X' * (n_columns) + '| } \n')
        output.write('\t \t \hline \n')
        output.write('\t \t')
        output.writelines([' & ' + re.sub('_', '\_', c) for c in df.columns])
        output.write(' \\\\ \n')
        output.write('\t \t \hline \n')
        for i in df.index:
            write_i = re.sub('%', '\%', str(i))
            write_i = re.sub('_', '\_', str(write_i))
            output.write('\t \t' + write_i + ' & ')
            if highlight_best_result:
                ##print("df :",df)
                for c in df.columns:
                    if df.loc[i, c] != df.loc[i, :].min():
                        if c != df.columns[-1]: # last columns of the table we don't need the ' & '
                            output.write(str(round(df.loc[i, c], 6)) + ' & ')
                        else:
                            output.write(str(round(df.loc[i, c], 6)))
                    else:
                        if c != df.columns[-1]:  # last columns of the table we don't need the ' & '
                            output.write('\\textbf{' + str(round(df.loc[i, c], 6)) + '}' + ' & ')
                        else:
                            output.write('\\textbf{' + str(round(df.loc[i, c], 6)) + '}')
            else:
                for c in df.columns:
                    if c != df.columns[-1]:
                        output.writelines([str(round(df.loc[i, c], 6)) + ' & '
                                   ])
                    else:
                        output.writelines([str(round(df.loc[i, c], 6))])
            output.write(' \\\\ \n')
            output.write('\t \t \hline \n')
        output.write('\t \end{tabularx} \n')
        if title:
            title_latex = re.sub('_','\_',title)
            output.write('\caption{' + title_latex + '}')
        else:
            output.write('\caption{Description}')
        output.write('\end{figure*} \n')
        output.write('\n')
        output.write('\end{document}')
    path_pdf_directory = latex_output_path
    subprocess.run(["pdflatex", "-output-directory", path_pdf_directory, path_tex])



def plot_adjacency_list(S, a_l,label =True):
    '''
    Plot the current adjacency list *a_l*.

    :param S: A stream graph (we get its labels)
    :param a_l: an adjacency list
    :return: Plot of adjacency list
    '''
    fig = plt.figure()
    ax = plt.axes()

    G = nx.Graph()
    set_edge = set()

    if type(list(a_l.keys())[0]) == int:
        for k, v in a_l.items():
            for i in v:
                if (k, i) not in set_edge and (i, k) not in set_edge:
                    if label:
                        G.add_edge(S.node_to_label[k], S.node_to_label[i])
                    else:
                        G.add_edge(k, i)
                    set_edge.add((k, i))
    else:
        for k, v in a_l.items():
            for i in v:
                if (k, i[1]) not in set_edge and (i[1], k) not in set_edge:
                    if label:
                        G.add_edge(S.node_to_label[k[2]], S.node_to_label[i[1][2]], t1=i[0])
                    else:
                        G.add_edge(k[2], i[1][2], t1=i[0])
                    set_edge.add((k, i[1]))
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=700,
                           node_color='#5a5ff5', alpha=0.5, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color='#2d5986',
                           alpha=0.3, width=5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=20, ax=ax)

    if type(list(a_l.keys())[0]) == tuple:
        nx.draw_networkx_edge_labels(G, pos, font_size=20, ax=ax)
    plt.xlabel("t", fontname='Ubuntu', fontsize=20, color='#476b6b')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top=False, bottom=False, right=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()
    # plt.show()


def nx_degree(G):
    return {n: G.degree[n] for n in G.nodes}
