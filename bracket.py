from optparse import OptionParser
from sympy import simplify,expand,symbols
import networkx as nx
import re

A,v,u=symbols("A,v,u")



''' 

Build PDs from Gauss code

'''


def build_pd(gc):
    PD = []
    checked = []
    for x in range(len(gc[0])):
        inv = gc[0].index(-1 * (gc[0][x]))
        if abs(gc[0][x]) not in checked:
            if gc[0][x] > 0 and gc[1][abs(gc[0][x]) - 1] == 1:
                PD.append([inv, x + 1, inv + 1, x])
            if gc[0][x] > 0 and gc[1][abs(gc[0][x]) - 1] == -1:
                PD.append([inv, x, inv + 1, x + 1])
            if gc[0][x] < 0 and gc[1][abs(gc[0][x]) - 1] == 1:
                PD.append([x, inv + 1, x + 1, inv])
            if gc[0][x] < 0 and gc[1][abs(gc[0][x]) - 1] == -1:
                PD.append([x, inv, x + 1, inv + 1])
            checked.append(abs(gc[0][x]))
    return PD


'''

 Build graph of regions

'''

def regions(pd):
    checked = set()
    region = set()
    regionlist = []
    nb_endpoints = 0
    for i in range(len(pd)):
        for k in range(4):
            if (i, k) not in checked:
                region = []
                i1 = i
                k1 = k
                while True:
                    checked.add((i1, k1))
                    if pd[i1][k1] > pd[i1][(k1 + 2) % 4]:
                        region.append((pd[i1][k1], 1))
                    else:
                        region.append((pd[i1][k1], -1))

                    found = False
                    for j in range(len(pd)):
                        for l in range(4):
                            if pd[i1][k1] == pd[j][l] and not (i1 == j and k1 == l):
                                i1 = j
                                k1 = (l + 1) % 4
                                found = True
                                break
                        if found:
                            break
                    if not found:
                        nb_endpoints += 1
                        if pd[i1][k1] > pd[i1][(k1 + 2) % 4]:
                            region.append((pd[i1][k1], -1))
                        else:
                            region.append((pd[i1][k1], 1))
                        if nb_endpoints > 2:
                            print ("too many endpoints")
                            return []
                        k1 = (k1 + 1) % 4
                    if i1 == i and k1 == k:
                        break
                regionlist.append(region)
    return regionlist


'''

    Recursive Function for states computation
    
    
'''




def state_list(partial_state,remaining_crosses,cross_signs):
    state=[]
    if len(remaining_crosses)==0:
        return partial_state
    else:
        remaining_crosses.pop(0)
        for x in partial_state:
            brack1 = [A*x[0], x[1]+[0]]
            brack2 = [A**(-1)*x[0], x[1]+[1]]
            state.append(brack1)
            state.append(brack2)
        return state_list(state,remaining_crosses,cross_signs)
                
#            
        

'''

    arcs in region function
    

'''

def crossing_regions(pd,regions_list):
    cross_regs=[]
    for x in pd:
        for reg in regions_list:
            if any(x[0]== r[0] for r in reg) and any(x[1]==r[0] for r in reg):
                r0=regions_list.index(reg)
            if any(x[1]==r[0] for r in reg) and any(x[2]==r[0] for r in reg):
                r1=regions_list.index(reg)
            if any(x[2]==r[0] for r in reg) and any(x[3]==r[0] for r in reg):
                r2=regions_list.index(reg)
            if any(x[3]==r[0] for r in reg) and any(x[0]==r[0] for r in reg):
                r3=regions_list.index(reg)                
                  
        cross_regs.append([r0,r1,r2,r3])
    return cross_regs


'''

   Connectivity Graph

'''

def build_graph(stat_lst,crs_lst):
    graph_list=[]
    for stat in stat_lst:
        G =nx.MultiGraph()
        G.add_nodes_from(range(max([j for i in crs_lst  for j in i])))
        for i in range(len(stat[1])):
            if G.has_edge(crs_lst[i][0],crs_lst[i][1]):
                if G[crs_lst[i][0]][crs_lst[i][1]][0]['weight']>1:
                    G.remove_edge(crs_lst[i][0],crs_lst[i][1])
                    G.add_edge(crs_lst[i][0],crs_lst[i][1],weight=1)
            
            elif not G.has_edge(crs_lst[i][0],crs_lst[i][1]) and crs_lst[i][0]!=crs_lst[i][1]:
                G.add_edge(crs_lst[i][0],crs_lst[i][1],weight=1)
            
            if G.has_edge(crs_lst[i][1],crs_lst[i][2]):
                if G[crs_lst[i][1]][crs_lst[i][2]][0]['weight']>1:
                    G.remove_edge(crs_lst[i][1],crs_lst[i][2])
                    G.add_edge(crs_lst[i][1],crs_lst[i][2],weight=1)
            
            elif not G.has_edge(crs_lst[i][1],crs_lst[i][2]) and crs_lst[i][1]!=crs_lst[i][2]:
                G.add_edge(crs_lst[i][1],crs_lst[i][2],weight=1)
            
            if G.has_edge(crs_lst[i][2],crs_lst[i][3]):
                if G[crs_lst[i][2]][crs_lst[i][3]][0]['weight']>1:
                    G.remove_edge(crs_lst[i][2],crs_lst[i][3])
                    G.add_edge(crs_lst[i][2],crs_lst[i][3],weight=1)
            
            elif not G.has_edge(crs_lst[i][2],crs_lst[i][3]) and crs_lst[i][2]!=crs_lst[i][3]:
                G.add_edge(crs_lst[i][2],crs_lst[i][3],weight=1)
            
            if G.has_edge(crs_lst[i][3],crs_lst[i][0]):
                if G[crs_lst[i][3]][crs_lst[i][0]][0]['weight']>1:  
                    G.remove_edge(crs_lst[i][3],crs_lst[i][0])
                    G.add_edge(crs_lst[i][3],crs_lst[i][0],weight=1)
            
            elif not G.has_edge(crs_lst[i][3],crs_lst[i][0]) and crs_lst[i][3]!=crs_lst[i][0]:
                G.add_edge(crs_lst[i][3],crs_lst[i][0],weight=1)

            if stat[1][i]==0:
                if not G.has_edge(crs_lst[i][0],crs_lst[i][2]) and crs_lst[i][0]!=crs_lst[i][2]:
                    G.add_edge(crs_lst[i][0],crs_lst[i][2],weight=2)
                
                if G.has_edge(crs_lst[i][1],crs_lst[i][3]):
                    if G[crs_lst[i][1]][crs_lst[i][3]][0]['weight']>0:
                        G.remove_edge(crs_lst[i][1],crs_lst[i][3])
                        G.add_edge(crs_lst[i][1],crs_lst[i][3],weight=0)
                elif not G.has_edge(crs_lst[i][1],crs_lst[i][3]) and crs_lst[i][1]!=crs_lst[i][3]:        
                    G.add_edge(crs_lst[i][1],crs_lst[i][3],weight=0)

            elif stat[1][i]==1:

                if G.has_edge(crs_lst[i][0],crs_lst[i][2]):
                    if G[crs_lst[i][0]][crs_lst[i][2]][0]['weight']>0:
                        G.remove_edge(crs_lst[i][0],crs_lst[i][2])
                        G.add_edge(crs_lst[i][0],crs_lst[i][2],weight=0)
                elif not G.has_edge(crs_lst[i][0],crs_lst[i][2]) and crs_lst[i][0]!=crs_lst[i][2]:
                    G.add_edge(crs_lst[i][0],crs_lst[i][2],weight=0)
                    
                if not G.has_edge(crs_lst[i][1],crs_lst[i][3]) and crs_lst[i][1]!=crs_lst[i][3]:                      
                    G.add_edge(crs_lst[i][1],crs_lst[i][3],weight=2)
            
        graph_list=graph_list+[G]
    return graph_list


'''
   
   Component Counter
   
   
'''

def component_counter(graph_list,stat_lst):
    if len(graph_list)!=len(stat_lst):
        raise "ERROR : graph list not the same length as states list."
    else:
        components_list=[]
        for g in graph_list:
            H = nx.MultiGraph(((source, target, attr) for source, target, attr in g.edges(data=True) if attr['weight']==0))
            H.add_nodes_from(set(g.nodes())- set(H.nodes()))
            components_list.append(nx.number_connected_components(H))
    return components_list

'''

    Find Outside Region

'''

def find_outside(gc,regions_list):
    for i in range(len(regions_list)):
        if set(gc[2])==set([m[0] for m in regions_list[i]]):
            return i



'''

    Find endpoint Regions
    
'''
        

def find_endpoint_regions(gc,regions_list):
    head_reg,tail_reg=None,None
    max_cross= max(gc[0])
    for i in range(len(regions_list)):
        if 0 in set([m[0] for m in regions_list[i]]):
            head_reg=i
        if 2*max_cross in set([m[0] for m in regions_list[i]]):
            tail_reg=i
    return (head_reg,tail_reg)
            
        
    
'''

   Find Dijkstra path
   
   
'''


def find_path(comp_list, graph,outside_region,end_regs):
    nested_comps=[]
    sum_path1,sum_path2 = 0 , 0
    for i in range(len(graph)):
        if comp_list[i]==1:
            nested_comps.append(0)
        elif comp_list[i]>1:
            path1=nx.dijkstra_path(graph[i],outside_region,end_regs[0])
            path2=nx.dijkstra_path(graph[i],outside_region,end_regs[1])
            if len(path1)==1:
                sum_path1= min(path1)
            else:
                for k in range(len(path1)-1):
                    for key in graph[i].get_edge_data(path1[k],path1[k+1]):
                        sum_path1 = sum_path1 + graph[i].get_edge_data(path1[k],path1[k+1])[key]['weight']
            
            if len(path2)==1:
                sum_path2= min(path2)
            else:
                for k in range(len(path2)-1):
                    for key in graph[i].get_edge_data(path2[k],path2[k+1]):
                        sum_path2 = sum_path2 + graph[i].get_edge_data(path2[k],path2[k+1])[key]['weight']
            nested_comps.append(min(sum_path1,sum_path2))
    return nested_comps
          
'''

   Regions Graph
   
'''

def unsmoothed(crs_lst,pd):
        U =nx.MultiGraph()
        
        U.add_nodes_from(range(max([j for i in crs_lst  for j in i])))
        for i in range(len(crs_lst)):
            if not U.has_edge(crs_lst[i][0],crs_lst[i][1]) and crs_lst[i][0]!=crs_lst[i][1]:
                U.add_edge(crs_lst[i][0],crs_lst[i][1],weight=1,arc=pd[i][1])
            if not U.has_edge(crs_lst[i][1],crs_lst[i][2]) and crs_lst[i][1]!=crs_lst[i][2]:
                U.add_edge(crs_lst[i][1],crs_lst[i][2],weight=1,arc=pd[i][2])
            if not U.has_edge(crs_lst[i][2],crs_lst[i][3]) and crs_lst[i][2]!=crs_lst[i][3]:
                U.add_edge(crs_lst[i][2],crs_lst[i][3],weight=1,arc=pd[i][3])
            if not U.has_edge(crs_lst[i][3],crs_lst[i][0]) and crs_lst[i][3]!=crs_lst[i][0]:
                U.add_edge(crs_lst[i][3],crs_lst[i][0],weight=1,arc=pd[i][0]) 
        return U

'''

   Add closure arc for extented bracket
   
'''

def add_closure(graph, end_regs,reg_list):
    crossed_arcs=[]
    path=nx.dijkstra_path(graph,end_regs[0],end_regs[1])
    for i in range(len(path)-1):
        crossed_arcs.append(graph.get_edge_data(path[i],path[i+1])[0]['arc'])
    cr_len = len(crossed_arcs)
    cr_set=[None]*cr_len
    for i in range(cr_len):
        for j in range(len(reg_list[path[i]])):
            if crossed_arcs[i] ==reg_list[path[i]][j][0]:
                cr_set[i]=reg_list[path[i]][j]
    return cr_set

 
'''

   Trace states

'''

def trace_states(pd,reg_list,stat_lst,graph_list):
    checked = set()
    trace_lst = []
    start=None
    origin_start=None
    found=False
    max_arc=2*len(pd)
    for cr in range(len(pd)):
        for arc in range(4):
            if pd[cr][arc]==0:
                origin_start=(cr,arc)
                found=True
                break
        if found:
            break
    
    for state in stat_lst:
        start=origin_start
        checked=set()
        trace = [(0,1)]
        
        while True:
            if pd[start[0]][start[1]]==max_arc:
                break
                
            if start not in checked:
                checked.add(start)
                i1=start[0]
                k1=start[1]
                found=False

                for j in range(len(pd)):
                    for l in range(4):
                        if pd[i1][k1] == pd[j][l] and not (i1 == j and k1 == l):
                            if state[1][j]==0:
                                if l==0 :
                                    found=True
                                    if (j,1) not in checked:
                                        l=1
                                        
                                    else:
                                        continue
                                elif l==1 :
                                    found=True
                                    if (j,0) not in checked:
                                        l=0
                                    else:
                                        continue
                                elif l==2 :
                                    found=True
                                    if (j,3) not in checked:
                                        l=3
                                    else:
                                        continue
                                elif l==3 :
                                    found=True
                                    if (j,2 not in checked):
                                        l=2
                                    else:
                                        continue
                                else:
                                    raise Exception ('Error: No condition met.')
                                i1=j
                                k1=l
                            else:
                                if l==0 :
                                    found=True
                                    if (j,3) not in checked:
                                        l=3
                                        
                                    else:
                                        continue
                                elif l==1 :
                                    found=True
                                    if (j,2) not in checked:
                                        l=2
                                    else:
                                        continue
                                elif l==2 :
                                    found=True
                                    if (j,1) not in checked:
                                        l=1
                                    else:
                                        continue
                                elif l==3 :
                                    found=True
                                    if (j,0 not in checked):
                                        l=0
                                    else:
                                        continue
                                else:
                                    raise Exception ('Error: No condition met.')
                                i1=j
                                k1=l
                    if found:
                        break
                            
                if not found:
                    endpoint=True
                if endpoint:
                    if state[1][i1]==0:
                        if k1==0 :
                            k1=1
                        
                        elif k1==1 :
                            k1=0
                       
                        elif k1==2 :
                            k1=3
                        
                        elif k1==3 :
                            k1=2
                    else:
                        if k1==0 :
                            k1=3
                        
                        elif k1==1 :
                            k1=2
                        
                        elif k1==2 :
                            k1=1
                        
                        elif k1==3 :
                            k1=0
                    endpoint=False
       
                if pd[i1][k1]>pd[i1][(k1 + 2) % 4]:
                    trace.append((pd[i1][k1],1))
                elif pd[i1][k1]<pd[i1][(k1 + 2) % 4]:
                    trace.append((pd[i1][k1],-1))
                else:
                    raise Exception ('Error: pd[i1][k1] ==(pd[i1][k1]+2)%4')
                start=(i1,k1)

        trace_lst.append(trace)
    return trace_lst



     
'''

   Compute Writhe
   
'''


def writhe(pd):
    writhe=0
    for i in range(len(pd)):
        writhe = writhe + (pd[i][1]-pd[i][3])
    return writhe


'''

   Compute Crossings' signs
   
'''


def cross_signs(pd):
    cross_signs=[]
    for i in range(len(pd)):
        if pd[i][1]-pd[i][3] >0:
            cross_signs.append(1)
        elif pd[i][1]-pd[i][3] <0:
            cross_signs.append(-1)
    return cross_signs

                                
'''
  
   Find sign changes

'''

def sign_changes(pd,comp_list,trace_lst,cr_arcs):
    max_arc = 2*len(pd)
    all_arcs=set(range(max_arc+1))
    
    if len(comp_list)!=len(trace_lst):
        raise Exception ('Error: comp_list != trace_lst')
    for i in range(len(comp_list)):
        if comp_list[i]>1:
            arcs_in_state=set([x[0] for x in trace_lst[i]])
            remaining_arcs = all_arcs - arcs_in_state
            for r_arc in remaining_arcs:
                trace_lst[i].append((r_arc,0))
    changes=[]
    for tr in trace_lst:
        tmp=[]
        for x in tr:
            for y in cr_arcs:
                if x[0] == y[0]:
                    if x[1]==1:
                        tmp.append(-y[1])
                        break
                    else:
                        tmp.append(y[1])
        changes.append(tmp)
    return [sum(x) for x in changes]
            
            
            
    

    
       
'''

   Build regular polynomial
   
'''


def build_poly(stat_lst, comp_list, nest_comps):
    poly_list=[]

    for k in range(len(stat_lst)):
        if planar:
            delta_exp=comp_list[k]-1-nest_comps[k]
            v_exp=nest_comps[k]
            poly_list.append(stat_lst[k][0]* (-A**2 - A**(-2))**(delta_exp) * v**v_exp)
        else:
            delta_exp=comp_list[k]-1
            poly_list.append(stat_lst[k][0]* (-A**2 - A**(-2))**(delta_exp))
        
    return simplify(sum(poly_list))


'''

   Build extended polynomial
   
'''


def build_ext_poly(stat_lst, comp_list, nest_comps,changes,cr_arcs):
    poly_list=[]

    for k in range(len(stat_lst)):
        if planar:
            delta_exp=comp_list[k]-1-nest_comps[k]
            v_exp=nest_comps[k]
            poly_list.append(stat_lst[k][0]* (-A**2 - A**(-2))**(delta_exp) * v**v_exp * u**changes[k])
        else:
            delta_exp=comp_list[k]-1
            poly_list.append(stat_lst[k][0]* (-A**2 - A**(-2))**(delta_exp) * u**changes[k])
        
    return simplify(sum(poly_list))


'''

    Main Function

'''

def main(gc,extented):
    brack=[[1,[]]]
    pd=build_pd(gc)
    regions_list=regions(pd)

    writh = writhe(pd)
    crs_sgn=cross_signs(pd)
    
    if planar:
        outside_region=find_outside(gc,regions_list)
    
    end_regs= find_endpoint_regions(gc,regions_list)
    
    remaining_crosses= [ x for x in pd]

    stat_lst=state_list(brack,remaining_crosses,crs_sgn)
        
    crs_reg=crossing_regions(pd,regions_list)
    
    graph_list = build_graph(stat_lst,crs_reg)
    
    comp_list= component_counter(graph_list,stat_lst)
    
    if planar:
        print ('Evaluating planar knotoid.')
        nest_comps=find_path(comp_list,graph_list,outside_region,end_regs)
    else:
        print ('Evaluating knotoid in S^2.')
        nest_comps=0
    
    no_sm_graph=unsmoothed(crs_reg,pd)
    cr_arcs=add_closure(no_sm_graph, end_regs,regions_list)

    if not cr_arcs:
        extented=False
    
    if extented:
        print ('Using extended {type} bracket.'.format(type="Turaev" if planar else ""))
        trace_list=trace_states(pd,regions_list,stat_lst,graph_list)
        changes=sign_changes(pd,comp_list,trace_list,cr_arcs)
        polynomial = build_ext_poly(stat_lst, comp_list, nest_comps,changes,cr_arcs)
    
        return expand((-A**3)**(-writh)*(u**(sum([x[1] for x in cr_arcs])))*polynomial )


    else:
        print ('Using {type} bracket'.format(type="Turaev" if planar else ""))
        polynomial= build_poly(stat_lst, comp_list, nest_comps)
    
        return expand((-A**3)**(-writh)*polynomial )


def str_to_gcode(gcodestr):
    global planar
    gcode = gcodestr.split(' ')
    gcode[0] = [int(x) for x in gcode[0].split(',')]
    hlp = []
    planar=False
    for i in range(len(gcode[1])):
        if gcode[1][i] == '+':
            hlp.append(1)
        if gcode[1][i] == '-':
            hlp.append(-1)
    gcode[1] = hlp
    if len(gcode)==3:
        gcode[2] = [int(x) for x in gcode[2].split(',')]
        planar=True
    return gcode



'''

    Main part
    
'''

if __name__ == "__main__":

    parser = OptionParser()   
    parser.add_option("-f", "--file", help='Specify input filename.', action="store",type="string",dest="filename")
    parser.add_option("-g", "--gausscode", help='Input single oriented Gauss code.\
 					Example: "1,-2,-1,2 +++ 0,2,3"\
                    The Gauss code of a planar knotoid diagram consists of three parts.\
                    First part includes crossings. Positive for overcrossing, negative for undercrossing\
                    Second part includes signs of crossings. +1 for positive crossing, -1 for negative crossing.\
                    Third part includes labels of arcs that touch the outer region of the diagram.\
                    A knotoid diagram on the sphere has only two parts; the part that indicates the arcs that touch the outer region is ommited.\
                  	WARNING: The quote characters " are required for this to work.', type="string",dest="gausscode")
    parser.add_option("-e", "--extended", help="Extended version of (Turaev) bracket.",action="store_true", default=False, dest="extended")
    parser.add_option("-o", "--output", help="Name of the output file.", action="store", type="string", default="results.txt", dest="output")
    (options,args)=parser.parse_args()
    optionsdict = vars(options)
    inputfilename=options.filename
    gc_inpt=options.gausscode
    extended=options.extended
    outputfilename=options.output


    if inputfilename and gc_inpt:
        raise Exception ('Error: Options -f and -g cannot be used simultaneously.')
    elif not inputfilename and not gc_inpt:
        raise Exception ('Error: No input is given.')
    elif not inputfilename and gc_inpt:
        print ('Single Gauss code input.')
        print ('Knotoid diagram: {}'.format(gc_inpt))
        print (main(str_to_gcode(gc_inpt),extended))
    elif inputfilename and not gc_inpt:
        print ('Multiple Gauss codes input.')
        file = open(inputfilename, 'r')
        file1 = open (outputfilename, 'w')
        file1.write("Gauss Code"+"\t"+"\t"+"\t"+"Polynomial"+"\n")
        for l in file:
            line = l.strip('\n')
            gcode = str_to_gcode(line)
            print ('Knotoid diagram: {}'.format(line))
            file1.write(str(gcode)+"\t"+str(main(gcode,extended))+"\n")
            
            
            
        file.close()
    


