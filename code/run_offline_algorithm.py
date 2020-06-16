import argparse as args
import os
import sys
import traceback
import constants
import time
import util
from interaction_graph import *
import algorithms

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))

basepath = os.path.dirname(os.path.abspath(__file__)) + "/../"
basePath_data = basepath + "datasets/"

def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false','f','n','0'):
            return False
        else:
            raise args.ArgumentTypeError('Boolean value expected.')

def create_parser():
    parser = args.ArgumentParser()
    parser.add_argument('-d', '--dataset', help="Dataset to use", type=str, required=True)
    parser.add_argument("-a", "--alg", help="Algorithm to use", choices=["MIL", "D-MIL"], required=True)
    parser.add_argument('-s', '--seed', help="Seed for sampling", type=int, default=100)
    parser.add_argument('-r', '--runs', help="Number of runs of the selected algorithm", type=int, default=50)
    parser.add_argument("--hill_climbing",type=str2bool, nargs='?',
                        const=True, default='f', help='It enables the hill climbing postprocessing step')
    parser.add_argument('-I', '--hill_climbing_iters', help="Maximum number of iterations for hill climbing algorithm", type=int, default=8)
    parser.add_argument("--save_only_avg",type=str2bool, nargs='?',
                        const=True, default='t', help='if True, only average results across several runs are saved instead of results for each run')
    return parser

def main(parsed):
    dataset = parsed.dataset
    dataset_path = basePath_data + dataset + "/processed/" + dataset + ".txt"
    if "syn" in dataset:
        if "W" in dataset:
            # -d synWatts-P-N, P is a float and N is an integer
            splitted = dataset.split("-")
            p = float(splitted[1])
            nei = int(splitted[2])
            dataset_path = basePath_data + "/syn_watts/p_" + str(p) + "/nei_" + str(nei) + "/" + "/syn.txt"
            dataset = "/syn_watts/p_" + str(p) + "/nei_" + str(nei)
        else:
            # -d syn-X, X is an integer
            splitted = dataset.split("-")
            id_net = int(splitted[1])
            dataset_path = basePath_data + "/syn/m_" + str(id_net) + "/syn.txt"
            dataset = "/syn/m_" + str(id_net)
    oracle = parsed.alg
    oracle_str = oracle
    relocation = parsed.hill_climbing
    save_only_avg = parsed.save_only_avg
    if relocation:
        oracle_str += "_R"
    # set the random seed for reproducibility, in algorithms.py the various algorithms use numpy library to sample integers
    seed = parsed.seed
    constants.seed = seed
    np.random.seed(seed)
    I = parsed.hill_climbing_iters
    constants.maximum_number_iterations_relocation = I
    n_runs = parsed.runs
    constants.num_iteration_offline = n_runs
    try:
        ig = InteractionGraph.get_instance_from_file(dataset_path)
        p_plus = ig.get_p_plus()
        p_minus = ig.get_p_minus()
        sum_loss = 0.0
        sum_nc = 0.0
        sum_time = 0.0
        for i in range(n_runs):
            file_res = basepath + "/output/"
            if not os.path.exists(file_res):
                os.makedirs(file_res)
            file_res += dataset + "/" + oracle_str
            log_path = file_res
            if not os.path.exists(file_res):
                os.makedirs(file_res)
            if not save_only_avg:
                file_res += "/i_" + str(i) + "/"
                if not os.path.exists(file_res):
                    os.makedirs(file_res)
                    
            start = time.time()
            if oracle == "D-MIL":
                if save_only_avg:
                    clustering = algorithms.min_cc_degree_ailon_heap_lazy(ig, p_plus, p_minus, relocation, file_path_times=None)
                else:
                    clustering = algorithms.min_cc_degree_ailon_heap_lazy(ig, p_plus, p_minus, relocation, file_path_times=file_res)
            elif oracle == "MIL":
                if save_only_avg:
                    clustering = algorithms.min_cc_ailon(ig, p_plus, p_minus, relocation, file_path_times=None)
                else:
                    clustering = algorithms.min_cc_ailon(ig, p_plus, p_minus, relocation, file_path_times=file_res)
            run_time = time.time() - start
            loss = ig.analytically_expected_loss(clustering, p_plus, p_minus)
            n_clusters = util.count_clusters(clustering)
            if not save_only_avg:
                # save the loss of the computed clustering
                file_loss = file_res + "loss.txt"
                with open(file_loss, "w+") as f:
                    f.write(str(loss)) 
                # save the number of clusters
                file_n_cl = file_res + "n_clusters.txt"
                with open(file_n_cl, "w+") as f:
                    f.write(str(n_clusters)) 
                # save the clustering 
                file_res += "class.txt"
                with open(file_res, "w+") as f:
                    for c in clustering:
                        f.write(str(c) + "\n") 
            print("Run #{} completed: \t loss = {:.2f}, \t #clusters = {}".format(i, loss, n_clusters))
            sum_time += run_time
            sum_loss += loss
            sum_nc += n_clusters
        avg_time = sum_time / n_runs
        avg_loss = sum_loss / n_runs
        avg_nc = float(sum_nc) / n_runs
        print("Avg. loss = {:.2f}, avg. #clusters = {}".format(avg_loss, avg_nc))
        file_avg_loss = basepath + "/output/" + dataset + "/" + oracle_str + "/avg_loss.txt"
        file_avg_nc = basepath + "/output/" + dataset + "/" + oracle_str + "/avg_n_clusters.txt"
        file_avg_time = basepath + "/output/" + dataset + "/" + oracle_str + "/avg_time.txt"
        with open(file_avg_loss, "w+") as f:
            f.write(str(avg_loss))
        with open(file_avg_nc, "w+") as f:
            f.write(str(avg_nc))
        with open(file_avg_time, "w+") as f:
            f.write(str(avg_time))
    except FileNotFoundError as e:
        print(str(e))
    except Exception as e:
        print(str(e))
        trace_str = traceback.format_exc()
        with open(log_path + "/error_log.txt", "w+") as f:
            f.write(trace_str)

if __name__ == '__main__':
    parsed = create_parser().parse_args()
    main(parsed)