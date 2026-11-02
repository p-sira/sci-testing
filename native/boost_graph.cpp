// Public Boost.Graph-derived comparisons for saved sparse.csgraph contracts.
//
// Inputs are explicit COO edges (including explicit zero-weight edges).
// This file uses PUBLIC Boost.Graph algorithms from Boost 1.83-compatible
// headers; it does not call SciPy, Sci, or private Boost implementation APIs.
// The contract is intentionally restricted to saved distances (not paths),
// canonical component labels, and MST weight, not arbitrary SciPy options.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <deque>

namespace {
thread_local std::string last_error;
using DG = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
    boost::no_property, boost::property<boost::edge_weight_t, double>>;
using UG = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, double>>;
constexpr int MAX_VERTICES = 256;
constexpr double INF = std::numeric_limits<double>::max();

bool named(const char* got, const char* expect) {
    return got && std::strcmp(got, expect) == 0;
}
bool path_api(const char* name) {
    return named(name,"dijkstra") || named(name,"bellman_ford") ||
           named(name,"floyd_warshall") || named(name,"johnson") ||
           named(name,"shortest_path");
}
bool supported(const char* name) {
    return path_api(name) || named(name,"connected_components") ||
           named(name,"minimum_spanning_tree");
}
void validate_options(const char* name, int n, const double* p, int np,
                      int n_indices) {
    if (!supported(name)) throw std::invalid_argument("unsupported public Boost.Graph comparison");
    if (n < 1 || n > MAX_VERTICES)
        throw std::domain_error("graph size must be between 1 and 256");
    if (!p || (path_api(name) && np != 2 &&
               !(named(name,"dijkstra") && np==4)) ||
        (named(name,"connected_components") && np != 2) ||
        (named(name,"minimum_spanning_tree") && np != 2))
        throw std::invalid_argument("invalid graph option signature");
    for(int i=0;i<np;++i)
        if(!std::isfinite(p[i]))
            throw std::domain_error("nonfinite graph option");
    if (p[0]!=0. && p[0]!=1.)
        throw std::domain_error("directed flag must be boolean");
    if (np>1 && p[1]!=0. && p[1]!=1.)
        throw std::domain_error("unweighted/connection flag must be boolean");
    if (np==4 && (p[2]!=0. && p[2]!=1. || p[3]<0.))
        throw std::domain_error("invalid Dijkstra min_only or limit");
    if (n_indices<0 || n_indices>n)
        throw std::domain_error("invalid number of source indices");
}
int output_count(const char* name,int n,const double* p,int np,int n_indices) {
    validate_options(name,n,p,np,n_indices);
    if(named(name,"connected_components")) return n+1;
    if(named(name,"minimum_spanning_tree")) return 1;
    bool min_only=named(name,"dijkstra") && np==4 && p[2]!=0.;
    return 2+n*(min_only ? 1 : (n_indices ? n_indices : n));
}
using EdgeMap = std::map<std::pair<int,int>,double>;
EdgeMap gather_edges(int n,const int* row,const int*col,
                     const double* value,int count,bool undirected,bool unweighted) {
    if(count<0 || count>n*n+256 || (count>0 && (!row||!col||!value)))
        throw std::invalid_argument("missing or invalid COO graph buffers");
    EdgeMap result;
    for(int i=0;i<count;++i) {
        if(row[i]<0 || row[i]>=n || col[i]<0 || col[i]>=n ||
           !std::isfinite(value[i]))
            throw std::domain_error("invalid graph edge coordinate/weight");
        int r=row[i],c=col[i];
        if(undirected && r>c) std::swap(r,c);
        const auto key=std::make_pair(r,c);
        const double w=unweighted ? 1. : value[i];
        // For undirected SciPy csgraph, symmetric weights take the minimum.
        const auto it=result.find(key);
        if(it==result.end()) result.emplace(key,w);
        else if(undirected) it->second=std::min(it->second,w);
        else it->second += w; // duplicate directed COO coordinates sum
    }
    return result;
}
DG make_directed(int n,const EdgeMap& edges,bool undirected) {
    DG graph(n);
    for(const auto& edge:edges) {
        const auto r=edge.first.first,c=edge.first.second;
        boost::add_edge(r,c,edge.second,graph);
        if(undirected && r!=c)boost::add_edge(c,r,edge.second,graph);
    }
    return graph;
}
UG make_undirected(int n,const EdgeMap& edges) {
    UG graph(n);
    for(const auto& edge:edges)
        if(edge.first.first!=edge.first.second)
            boost::add_edge(edge.first.first,edge.first.second,edge.second,graph);
    return graph;
}
std::vector<double> one_source(const DG& graph,int n,int source,
                               bool bellman,double limit) {
    std::vector<double> distances(static_cast<size_t>(n), INF);
    distances[source]=0.;
    auto index=get(boost::vertex_index,graph);
    auto dist=boost::make_iterator_property_map(distances.begin(),index);
    if(bellman) {
        if(!boost::bellman_ford_shortest_paths(graph, n,
            boost::weight_map(get(boost::edge_weight,graph))
                  .distance_map(dist)
                  .distance_inf(INF).distance_zero(0.)))
            throw std::domain_error("negative graph cycle");
    } else {
        boost::dijkstra_shortest_paths(graph,source,
            boost::weight_map(get(boost::edge_weight,graph))
                  .distance_map(dist)
                  .distance_inf(INF).distance_zero(0.));
    }
    for(double& value:distances)
        if(value>=INF || value>limit)
            value=std::numeric_limits<double>::infinity();
    return distances;
}
std::vector<std::vector<double>> all_pairs(const DG& graph,int n,
                                           const char* name) {
    std::vector<std::vector<double>> d(n,std::vector<double>(n,INF));
    bool ok;
    if(named(name,"floyd_warshall"))
        ok=boost::floyd_warshall_all_pairs_shortest_paths(graph,d,
             boost::weight_map(get(boost::edge_weight,graph))
                   .distance_inf(INF).distance_zero(0.));
    else
        ok=boost::johnson_all_pairs_shortest_paths(graph,d,
             boost::weight_map(get(boost::edge_weight,graph))
                   .distance_inf(INF).distance_zero(0.));
    if(!ok)throw std::domain_error("negative graph cycle");
    for(auto& row:d)
        for(double& x:row)
            if(x>=INF)x=std::numeric_limits<double>::infinity();
    return d;
}

// The second entry point preserves the seven established ABI contracts.
// Its input is a bounded rectangular dense view, because saved maximum
// bipartite matching and structural rank include genuine 2x3/3x2 matrices.
// Output length is reported separately for BFS/DFS's reachable-only order.
bool extra_supported(const char* name) {
    return named(name,"breadth_first_order") ||
           named(name,"breadth_first_tree") ||
           named(name,"depth_first_order") ||
           named(name,"depth_first_tree") ||
           named(name,"maximum_flow") ||
           named(name,"maximum_bipartite_matching") ||
           named(name,"structural_rank") ||
           named(name,"laplacian");
}
bool bfs(const char* name) {
    return named(name,"breadth_first_order") || named(name,"breadth_first_tree");
}
bool traversal(const char* name) {
    return bfs(name) || named(name,"depth_first_order") ||
           named(name,"depth_first_tree");
}
bool tree(const char* name) {
    return named(name,"breadth_first_tree") || named(name,"depth_first_tree");
}
struct VisitState {
    std::vector<int>* order;
    std::vector<int>* pred;
};
struct BFSRecorder : boost::default_bfs_visitor {
    VisitState state;
    explicit BFSRecorder(VisitState s):state(s) {}
    template<class Vertex,class Graph>
    void discover_vertex(Vertex v,const Graph&) const {
        state.order->push_back(static_cast<int>(v));
    }
    template<class Edge,class Graph>
    void tree_edge(Edge e,const Graph& g) const {
        (*state.pred)[static_cast<size_t>(boost::target(e,g))] =
            static_cast<int>(boost::source(e,g));
    }
};
struct DFSRecorder : boost::default_dfs_visitor {
    VisitState state;
    explicit DFSRecorder(VisitState s):state(s) {}
    template<class Vertex,class Graph>
    void discover_vertex(Vertex v,const Graph&) const {
        state.order->push_back(static_cast<int>(v));
    }
    template<class Edge,class Graph>
    void tree_edge(Edge e,const Graph& g) const {
        (*state.pred)[static_cast<size_t>(boost::target(e,g))] =
            static_cast<int>(boost::source(e,g));
    }
};
using Capacity = long long;
using FlowGraph = boost::adjacency_list<boost::vecS,boost::vecS,boost::directedS,
    boost::no_property,
    boost::property<boost::edge_capacity_t,Capacity,
    boost::property<boost::edge_residual_capacity_t,Capacity,
    boost::property<boost::edge_reverse_t,
        boost::adjacency_list_traits<boost::vecS,boost::vecS,
                                      boost::directedS>::edge_descriptor>>>>;
void add_flow_arc(FlowGraph& g,int from,int to,Capacity value) {
    const auto f=boost::add_edge(from,to,g).first;
    const auto r=boost::add_edge(to,from,g).first;
    auto caps=boost::get(boost::edge_capacity,g);
    auto reverse=boost::get(boost::edge_reverse,g);
    boost::put(caps,f,value);
    boost::put(caps,r,0);
    boost::put(reverse,f,r);
    boost::put(reverse,r,f);
}
void extra_validate(const char* name,int nr,int nc,
                    const double* matrix,int matrix_count,
                    const double* p,int np) {
    if(!extra_supported(name))throw std::invalid_argument("unsupported Boost.Graph extended API");
    if(nr<1 || nc<1 || nr>MAX_VERTICES || nc>MAX_VERTICES ||
       static_cast<long long>(nr)*nc!=matrix_count || !matrix)
        throw std::invalid_argument("invalid bounded rectangular graph matrix");
    const int wanted=traversal(name) || named(name,"maximum_flow") ? 2 :
        named(name,"laplacian") ? 3 : 1;
    if(!p || np!=wanted)
        throw std::invalid_argument("invalid extended graph option arity");
    for(int i=0;i<np;++i)
        if(!std::isfinite(p[i]))throw std::domain_error("nonfinite extended graph option");
    for(int i=0;i<matrix_count;++i)
        if(!std::isfinite(matrix[i]))
            throw std::domain_error("nonfinite extended graph matrix edge");
    if(traversal(name) || named(name,"maximum_flow") || named(name,"laplacian")) {
        if(nr!=nc)throw std::invalid_argument("this graph algorithm requires a square matrix");
    }
    if(traversal(name)) {
        if(p[0]<0 || p[0]>=nr || std::trunc(p[0])!=p[0] ||
           (p[1]!=0. && p[1]!=1.))
            throw std::domain_error("invalid start vertex or directed flag");
    } else if(named(name,"maximum_flow")) {
        if(std::trunc(p[0])!=p[0] || std::trunc(p[1])!=p[1] ||
           p[0]<0 || p[0]>=nr || p[1]<0 || p[1]>=nr || p[0]==p[1])
            throw std::domain_error("invalid source/sink");
    } else if(named(name,"laplacian")) {
        for(int i=0;i<3;++i)
            if(p[i]!=0. && p[i]!=1.)
                throw std::domain_error("invalid laplacian boolean option");
    } else if(p[0]!=0. && p[0]!=1.) {
        throw std::domain_error("invalid matching permutation option");
    }
}
int extra_capacity(const char* name,int nr,int nc) {
    if(traversal(name))
        return tree(name) ? 2+nr*nr : 2*nr;
    if(named(name,"laplacian"))return 2+nr*nr+nr;
    return 1;
}

std::vector<double> compute_extra(const char* name,int nr,int nc,
                                  const double* a,const double* p) {
    if(traversal(name)) {
        const bool undirected=p[1]==0.;
        // SciPy CSR traversal scans the source vertex's outgoing row first,
        // *then* its unseen incoming neighbors when directed=False. Sorting
        // the symmetric union instead changes BFS/DFS tie order and parent
        // edges on even small graphs. Populate each BGL vertex in that order.
        DG g(nr);
        for(int i=0;i<nr;++i) {
            for(int j=0;j<nc;++j)
                if(a[i*nc+j]!=0.)boost::add_edge(i,j,a[i*nc+j],g);
            if(undirected)
                for(int j=0;j<nr;++j)
                    if(j!=i && a[i*nc+j]==0. && a[j*nc+i]!=0.)
                        boost::add_edge(i,j,a[j*nc+i],g);
        }
        std::vector<int> order,pred(static_cast<size_t>(nr),-9999);
        const auto start=static_cast<size_t>(p[0]);
        const VisitState state{&order,&pred};
        if(bfs(name)) {
            std::vector<boost::default_color_type> colors(nr);
            boost::breadth_first_search(g,start,
                boost::visitor(BFSRecorder(state)).color_map(
                    boost::make_iterator_property_map(colors.begin(),
                        get(boost::vertex_index,g))));
        } else {
            std::vector<boost::default_color_type> colors(nr,boost::white_color);
            boost::depth_first_visit(g,start,DFSRecorder(state),
                boost::make_iterator_property_map(colors.begin(),
                    get(boost::vertex_index,g)));
        }
        if(tree(name)) {
            std::vector<double> result(static_cast<size_t>(2+nr*nr),0.);
            result[0]=result[1]=static_cast<double>(nr);
            for(int j=0;j<nr;++j)if(pred[j]>=0) {
                const int i=pred[j];
                // For directed=False SciPy takes the lower weight when
                // antiparallel input edges both exist. The traversal still
                // follows outgoing-first CSR enumeration, independently of
                // which direction supplies the tree edge's stored weight.
                const double outgoing=a[i*nc+j];
                const double incoming=a[j*nc+i];
                const double w=undirected && outgoing!=0. && incoming!=0.
                    ? std::min(outgoing,incoming)
                    : outgoing!=0. ? outgoing : incoming;
                result[2+i*nr+j]=w;
            }
            return result;
        }
        std::vector<double> result;
        result.reserve(order.size()+nr);
        for(int v:order)result.push_back(static_cast<double>(v));
        for(int predecessor:pred)result.push_back(static_cast<double>(predecessor));
        return result;
    }
    if(named(name,"maximum_bipartite_matching") ||
       named(name,"structural_rank")) {
        UG g(nr+nc);
        for(int i=0;i<nr;++i)
            for(int j=0;j<nc;++j)
                if(a[i*nc+j]!=0.)boost::add_edge(i,nr+j,g);
        using V=boost::graph_traits<UG>::vertex_descriptor;
        std::vector<V> mate(static_cast<size_t>(nr+nc));
        auto map=boost::make_iterator_property_map(mate.begin(),
                                                   get(boost::vertex_index,g));
        boost::edmonds_maximum_cardinality_matching(g,map);
        return {static_cast<double>(boost::matching_size(g,map))};
    }
    if(named(name,"maximum_flow")) {
        FlowGraph graph(nr);
        for(int i=0;i<nr;++i)
            for(int j=0;j<nc;++j) {
                const double val=a[i*nc+j];
                if(val<0. || std::trunc(val)!=val ||
                   val>static_cast<double>(std::numeric_limits<Capacity>::max()/nr))
                    throw std::domain_error("max flow needs nonnegative integral capacities");
                if(val>0.)add_flow_arc(graph,i,j,static_cast<Capacity>(val));
            }
        const Capacity value=boost::push_relabel_max_flow(
            graph,static_cast<int>(p[0]),static_cast<int>(p[1]));
        return {static_cast<double>(value)};
    }
    // BGL-derived directed graph degree and edge traversal. For normalized
    // Laplacian, SciPy's return_diag uses sqrt(degree), with isolated nodes
    // represented by 1 for safe normalization; self-loops cancel in L.
    DG g(nr);
    for(int i=0;i<nr;++i)
        for(int j=0;j<nc;++j)
            if(a[i*nc+j]!=0.)boost::add_edge(i,j,a[i*nc+j],g);
    std::vector<double> matrix(static_cast<size_t>(nr*nr),0.);
    auto weights=get(boost::edge_weight,g);
    auto edges=boost::edges(g);
    for(auto e=edges.first;e!=edges.second;++e) {
        const int i=static_cast<int>(boost::source(*e,g));
        const int j=static_cast<int>(boost::target(*e,g));
        const double w=get(weights,*e);
        matrix[i*nr+j]-=w;
        if(p[2]!=0.)matrix[j*nr+i]-=w;
    }
    std::vector<double> degree(static_cast<size_t>(nr),0.);
    for(int i=0;i<nr;++i)
        for(int j=0;j<nr;++j) {
            const double v=matrix[i*nr+j];
            if(i!=j)degree[p[1]!=0. ? i:j]-=v;
        }
    // Self-loops cancel between D and A; only off-diagonal degree remains.
    for(int i=0;i<nr;++i)matrix[i*nr+i]=degree[i];
    std::vector<double> diagonal=degree;
    if(p[0]!=0.) {
        for(double& d:diagonal) {
            if(d<0.)throw std::domain_error("negative degree prevents normalized Laplacian");
            d=d==0. ? 1. : std::sqrt(d);
        }
        for(int i=0;i<nr;++i)
            for(int j=0;j<nr;++j)
                matrix[i*nr+j]/=diagonal[i]*diagonal[j];
    }
    std::vector<double> result;
    result.reserve(static_cast<size_t>(nr*nr+nr+2));
    result.push_back(nr);result.push_back(nr);
    result.insert(result.end(),matrix.begin(),matrix.end());
    result.insert(result.end(),diagonal.begin(),diagonal.end());
    return result;
}
} // namespace

extern "C" {
const char* boost_graph_error(){return last_error.c_str();}
int boost_graph_supported(const char* name){return supported(name)?1:0;}
int boost_graph_output_count(const char* name,int n,
                             const double* opts,int nopts,int nindices) {
    try {return output_count(name,n,opts,nopts,nindices);}
    catch(...) {return 0;}
}
int boost_graph_evaluate(const char* name,int n,
                          const int* rows,const int* cols,const double* values,int nedges,
                          const double* opts,int nopts,const int* indices,int nindices,
                          double* out,int capacity) {
    last_error.clear();
    try {
        const int need=output_count(name,n,opts,nopts,nindices);
        if(!out||capacity<need)throw std::invalid_argument("insufficient graph output capacity");
        if(nindices>0 && !indices)throw std::invalid_argument("missing source indices");
        for(int i=0;i<nindices;++i)
            if(indices[i]<0||indices[i]>=n)
                throw std::domain_error("invalid graph source vertex");
        const bool is_mst=named(name,"minimum_spanning_tree");
        const bool is_connected=named(name,"connected_components");
        const bool undirected=is_mst || (!is_connected && opts[0]==0.) ||
            (is_connected && (opts[0]==0. || opts[1]==0.));
        const bool unweighted=path_api(name) && opts[1]!=0.;
        auto edges=gather_edges(n,rows,cols,values,nedges,undirected,unweighted);
        std::vector<double> answer(static_cast<size_t>(need));
        if(is_mst) {
            auto graph=make_undirected(n,edges);
            using E=boost::graph_traits<UG>::edge_descriptor;
            std::vector<E> tree;
            boost::kruskal_minimum_spanning_tree(graph,std::back_inserter(tree),
                     boost::weight_map(get(boost::edge_weight,graph)));
            double total=0.;
            for(const E& e:tree)total+=get(boost::edge_weight,graph,e);
            answer[0]=total;
        } else if(is_connected) {
            std::vector<int> labels(n);
            int count=0;
            if(undirected) {
                auto graph=make_undirected(n,edges);
                count=boost::connected_components(graph,&labels[0]);
            } else {
                auto graph=make_directed(n,edges,false);
                count=boost::strong_components(graph,
                      boost::make_iterator_property_map(labels.begin(),get(boost::vertex_index,graph)));
            }
            std::map<int,int> normalized;
            int next=0;
            answer[0]=static_cast<double>(count);
            for(int i=0;i<n;++i) {
                auto inserted=normalized.emplace(labels[i],next);
                if(inserted.second)++next;
                answer[1+i]=static_cast<double>(inserted.first->second);
            }
        } else {
            auto graph=make_directed(n,edges,undirected);
            const bool min_only=named(name,"dijkstra") && nopts==4 && opts[2]!=0.;
            const double limit=nopts==4?opts[3]:std::numeric_limits<double>::infinity();
            const int nrows=min_only?1:(nindices?nindices:n);
            answer[0]=static_cast<double>(nrows);
            answer[1]=static_cast<double>(n);
            const bool negative=std::any_of(edges.begin(),edges.end(),
                [](const auto& entry){return entry.second<0.;});
            if(named(name,"dijkstra") && negative)
                throw std::domain_error("negative edge weights unsupported by Dijkstra");
            if(named(name,"floyd_warshall") || named(name,"johnson")) {
                const auto distances=all_pairs(graph,n,name);
                for(int i=0;i<nrows;++i) {
                    const int source=nindices?indices[i]:i;
                    for(int j=0;j<n;++j)answer[2+i*n+j]=distances[source][j];
                }
            } else {
                for(int i=0;i<(nindices?nindices:n);++i) {
                    const int source=nindices?indices[i]:i;
                    auto distances=one_source(graph,n,source,
                        named(name,"bellman_ford") || negative,limit);
                    for(int j=0;j<n;++j) {
                        double& destination=answer[2+(min_only?0:i)*n+j];
                        if(min_only) {
                            if(i==0 || distances[j]<destination)destination=distances[j];
                        } else destination=distances[j];
                    }
                }
            }
        }
        for(int i=0;i<need;++i)out[i]=answer[i];
        return 0;
    } catch(const std::exception& ex){last_error=ex.what();return -1;}
      catch(...){last_error="unknown public Boost.Graph evaluation failure";return -1;}
}

int boost_graph_extra_supported(const char* name) {
    return extra_supported(name) ? 1 : 0;
}
int boost_graph_extra_output_capacity(const char* name,int nrows,int ncols) {
    if(!extra_supported(name) || nrows<1 || ncols<1 ||
       nrows>MAX_VERTICES || ncols>MAX_VERTICES)
        return 0;
    if((traversal(name) || named(name,"laplacian") ||
        named(name,"maximum_flow")) && nrows!=ncols)
        return 0;
    return extra_capacity(name,nrows,ncols);
}
int boost_graph_extra_evaluate(const char* name,int nrows,int ncols,
                              const double* dense,int dense_count,
                              const double* opts,int nopts,
                              double* output,int capacity,int* written) {
    last_error.clear();
    try {
        extra_validate(name,nrows,ncols,dense,dense_count,opts,nopts);
        const int required=extra_capacity(name,nrows,ncols);
        if(!output || capacity<required || !written)
            throw std::invalid_argument("insufficient extended graph output capacity");
        const auto result=compute_extra(name,nrows,ncols,dense,opts);
        if(result.size()>static_cast<size_t>(capacity))
            throw std::runtime_error("extended graph output exceeds capacity");
        std::copy(result.begin(),result.end(),output);
        *written=static_cast<int>(result.size());
        return 0;
    } catch(const std::exception& ex) {
        last_error=ex.what();return -1;
    } catch(...) {
        last_error="unknown extended public Boost.Graph failure";return -1;
    }
}
} // extern "C"
