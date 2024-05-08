import multiprocessing
import os
import sys
import time
import numpy as np
import argparse
import random
import pandas as pd
import os
import gc
TF_ENABLE_DEPRECATION_WARNINGS=1

def sweep(x):
    div = x.sum(axis=1, dtype='float')
    res = x/div
    where_are_NaNs = np.isnan(res)
    res[where_are_NaNs] = 0
    return res

def scale(y, c=True, sc=True):
    x = y.copy()

    if c:
        x -= x.mean()
    if sc and c:
        x /= x.std()
    elif sc:
        x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
    return x

def shape(tensor):
    s = tensor.get_shape()
    return tuple([s[i].value for i in range(0, len(s))])

def param_args():
    """
    Get the arguments for inference to run the script and store them in the args object.
    
    Args:
        None
    Return:
        args object with arguments
    """
    date_str = time.strftime('%m%d')
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--n', type=int, default=2, help="number of multi-processing")
    parser.add_argument('--input', type=str, default='sample_data/',help="path for input(enrichment, adjacency matrices)")
    parser.add_argument('--output', type=str, default='sample_data/out/',help="path for output")
    parser.add_argument('--order', type=int, default=1,help="number of order")
    args = parser.parse_args()
    return args

def ssc_cpu(num_processes,num,sample_name):
    df = pd.read_csv(args.input + 'enrich/' + sample_name + '.csv')
#         adj = np.asmatrix(pd.read_csv('csv/adj_' + l + '/' + st_gsea_all_go_names[a] + '.csv'))
    adj_csv = pd.read_csv(args.input + 'adj_' + str(l) + '/' + sample_name + '.csv')
    adj = np.asmatrix((adj_csv).drop(columns="Unnamed: 0").values)
#             print("Sample name is: " + st_gsea_all_go.names[a])
#             print("Adjacency matrix: " + order_01_trimmed_adj_matrices_list.names[b])

    if len(df.columns.values.tolist()) != len(adj_csv.columns.values.tolist()):
        raise NotImplementedError("Error: Length not same!")
    for i in list(range(len(df.columns.values.tolist()))):
        if df.columns.values.tolist()[i] != adj_csv.columns.values.tolist()[i]:
            raise NotImplementedError("Error: Barcode not same!")
            break

    pathways_name = list(np.concatenate(df[["Unnamed: 0"]].values.tolist()))
    pathways = np.asmatrix(df.drop(columns="Unnamed: 0").values)
    pl = len(pathways)
    row = list(range(pl))
    p_val = []
#         n=1
#         para = [row[i * n:(i + 1) * n] for i in range((len(row) + n - 1) // n )]

    num_groups = num_processes

    # Remove brackets and split elements
    elements = row

    # Calculate length of each group
    group_length = len(elements) // num_groups
    remainder = len(elements) % num_groups  # Calculate remainder

    # Initialize the starting index
    start_index = 0

    # Split into groups
    result = []
    for i in range(num_groups):
        # Calculate the end index for the current group
        end_index = start_index + group_length + (1 if i < remainder else 0)
        # Add the sublist to the result
        result.append(elements[start_index:end_index])
        # Update the starting index for the next group
        start_index = end_index

    for i in result[num]:
        try:
            df_tmp = pd.read_csv(args.output + str(l) + '/' + sample_name + '_' + str(l) + '/' + sample_name + '_' + str(l) + '_' +str(i) + '.csv')
            if len(df_tmp) == pow(len(df),1):
                pass
            else:
                for j in list(range(pl)):
                    pathway_i = pathways_name[i]
                    pathway_j = pathways_name[j]
                    combo_name = pathway_i + '_x_' + pathway_j
                    if combo_name in  df_tmp['combo_name'].tolist():
                        continue
                    else:
#                         config = tf.ConfigProto(graph_options=tf.GraphOptions(optimizer_options=tf.OptimizerOptions(opt_level=tf.OptimizerOptions.L0)))
#                         config.gpu_options.visible_device_list = "0"
#                         sess = tf.Session(config=config)
#                         sess.run(tf.global_variables_initializer())
                        obj1 = np.array(scale(pathways[i,]).getT(),dtype=np.float32)
                        obj2 = np.array(scale(pathways[j,]),dtype=np.float32)
                        prod = np.matmul(obj1, obj2)
                        adj2 = np.array(sweep(adj),dtype=np.float32)
                        local_scc = np.einsum('ij,ij->i', adj2,prod)
#                         local_scc_ = sess.run(local_scc)
                        local_scc_list = list(local_scc)
                        global_scc = np.float32(np.mean(local_scc_list))
                        compare = []
                        num_row = list(range(max(obj1.shape)))
                        del prod
                        del local_scc
                        gc.collect()
#                         sess.close()
            #                     print("================================================ Permutation Test ===================================================================")
                        #start = time.time()
                        for k in list(range(100)):
#                             sess = tf.compat.v1.Session(config=config)
                            rng = np.random.default_rng(seed = k)
                            new_order = rng.permutation(num_row,0)
                            x = obj1[new_order]
                            y = np.transpose(obj2)
                            y = y[new_order]
                            y = np.transpose(y)
                            prod_shuff = np.tensordot(x, y, axes=1)
                            final_shuff = np.einsum('ij,ij->i', adj2, prod_shuff)
                            scc_shuff = np.mean(final_shuff,0)
#                             scc_shuff_ = sess.run(scc_shuff)
                            scc_shuff_ = scc_shuff
                            compare.append(scc_shuff_)
#                             sess.close()
                            del rng
                            del new_order
                            del scc_shuff
                            del scc_shuff_
                            del final_shuff
                            del prod_shuff
                            del x
                            del y
                            gc.collect()
                        del obj1
                        del obj2
                        compare = list(compare)
                        pval1=0
                        pval2=0
                        if global_scc < 0:
                            for p in compare:
                                if p < global_scc:
                                    pval1 = pval1 + 1
                        else:
                            for p in compare:
                                if p >= global_scc:
                                    pval2 = pval2 + 1

                        p_val = (pval1 + pval2)/100
                        gc.collect()
#                         tf.keras.backend.clear_session()

                        path = args.output + str(l) + '/' + sample_name + '_' + str(l) 
                        isExist = os.path.exists(path)
                        if not isExist:
                            os.makedirs(path)
                        if j == 0:
                            df2 = 0
                            df2 = pd.DataFrame({"combo_name":[combo_name],'local_scc':[local_scc_list], 'global_scc':[global_scc], 'permutation':[compare],'p_val':[p_val]})
                            df2.to_csv(path + '/' + sample_name + '_' + str(l) + '_' +str(i) + '.csv',index = None)
                        else:
                            df_tmp.loc[len(df_tmp.index)] = [combo_name, local_scc_list, global_scc, compare, p_val]
                            df_tmp.to_csv(path + '/' + sample_name + "_" + str(l) + '_' +str(i) + '.csv',index = None)
                        #end = time.time()
                        #elapsed = (end - start)
                        #print(elapsed)
        except:

            for j in list(range(pl)):

                pathway_i = pathways_name[i]
                pathway_j = pathways_name[j]
                combo_name = pathway_i + '_x_' + pathway_j
    #                     print("====================================================================================================================================")
    #                     print("= Sample name: " + st_gsea_all_go_names[a])
    #                     print("= Adjacency matrix: " + order_01_trimmed_adj_matrices_list.names[b])
    #                     print("= Combo name: [ %s ]"%(combo_name))
    #                     print("====================================================================================================================================")
#                 config = tf.ConfigProto(graph_options=tf.GraphOptions(optimizer_options=tf.OptimizerOptions(opt_level=tf.OptimizerOptions.L0)))
#                 config.gpu_options.visible_device_list = "0"
#                 sess = tf.Session(config=config)
#                 sess.run(tf.global_variables_initializer())
                obj1 = np.array(scale(pathways[i,]).getT(),dtype=np.float32)
                obj2 = np.array(scale(pathways[j,]),dtype=np.float32)
                prod = np.matmul(obj1, obj2)
                adj2 = np.array(sweep(adj),dtype=np.float32)
                local_scc = np.einsum('ij,ij->i', adj2,prod)
#                 local_scc_ = sess.run(local_scc)
                local_scc_list = list(local_scc)
                global_scc = np.float32(np.mean(local_scc_list))
                compare = []
                num_row = list(range(max(obj1.shape)))
                del prod
                del local_scc
                gc.collect()
#                 sess.close()
    #                     print("================================================ Permutation Test ===================================================================")
                #start = time.time()
                for k in list(range(100)):
#                     sess = tf.compat.v1.Session(config=config)
                    rng = np.random.default_rng(seed = k)
                    new_order = rng.permutation(num_row,0)
                    x = obj1[new_order]
                    y = np.transpose(obj2)
                    y = y[new_order]
                    y = np.transpose(y)
                    prod_shuff = np.tensordot(x, y, axes=1)
                    final_shuff = np.einsum('ij,ij->i', adj2, prod_shuff)
                    scc_shuff = np.mean(final_shuff,0)
#                     scc_shuff_ = sess.run(scc_shuff)
                    scc_shuff_ = scc_shuff
                    compare.append(scc_shuff_)
#                     sess.close()
                    del rng
                    del new_order
                    del scc_shuff
                    del scc_shuff_
                    del final_shuff
                    del prod_shuff
                    del x
                    del y
                    gc.collect()

                del obj1
                del obj2
                compare = list(compare)
                pval1=0
                pval2=0
                if global_scc < 0:
                    for p in compare:
                        if p < global_scc:
                            pval1 = pval1 + 1
                else:
                    for p in compare:
                        if p >= global_scc:
                            pval2 = pval2 + 1

                p_val = (pval1 + pval2)/100

                gc.collect()
#                 tf.keras.backend.clear_session()

                path = args.output + str(l) + '/' + sample_name + '_' + str(l) 
                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                if j == 0:
                    df2 = 0
                    df2 = pd.DataFrame({"combo_name":[combo_name],'local_scc':[local_scc_list], 'global_scc':[global_scc], 'permutation':[compare],'p_val':[p_val]})
                    df2.to_csv(path + '/' + sample_name + "_" + str(l) + '_' +str(i) + '.csv',index = None)
                else:
                    df2.loc[len(df2.index)] = [combo_name, local_scc_list, global_scc, compare, p_val]
                    df2.to_csv(path + '/' + sample_name + "_" + str(l) + '_' +str(i) + '.csv',index = None)

                #end = time.time()
                #elapsed = (end - start)
                #print(elapsed)
                 
parser = argparse.ArgumentParser(description='')
args = param_args()

sample_name = os.listdir(args.input + 'enrich/')
sample_name = [sn.split('.csv', 1)[0] for sn in sample_name]
l = args.order
num_processes = args.n
df = pd.read_csv(args.input + 'enrich/' + sample_name[0] + '.csv')
if len(df) < num_processes:
    num_processes = len(df)
processes = []
print("-------- START --------")
print("SampleList: " + str(sample_name))
print("nCPU: " + str(num_processes))

for sn in sample_name:
    print("In progress: " + sn)
    for num in range(num_processes):
        process = multiprocessing.Process(target=ssc_cpu,args=(num_processes,num,sn))
        process.start()
        processes.append(process)

    # Wait for all processes to finish
    for process in processes:
        process.join()
    df_01 = pd.read_csv(args.output + str(l) + '/' + sn + "_" + str(l) + '/' + sn + "_" + str(l) + '_' +str(0) + '.csv')
    for i in range(1,len(df)):
        df_02 = pd.read_csv(args.output + str(l) + '/' + sn + "_" + str(l) + '/' + sn + "_" + str(l) + '_' +str(i) + '.csv')
        df_01 = pd.concat([df_01,df_02])
    isExist = os.path.exists(args.output + str(l) + '/' + sn + "_" + str(l) + '.csv')
    if isExist:
        print("The Output for [" + sn + ", order: " + str(l) + "] already saved as: [" + args.output + str(l) + '/' + sn + "_" + str(l) + '.csv]')
    if not isExist:    
        df_01.to_csv(args.output + str(l) + '/' + sn + "_" + str(l) + '.csv',index = None)
        print("The Output for [" + sn + ", order: " + str(l) + "] saved as: [" + args.output + str(l) + '/' + sn + "_" + str(l) + '.csv]')

