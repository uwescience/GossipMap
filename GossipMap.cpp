#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <string>

#include "malloc.h"
#include "google/malloc_extension.h"

#include <graphlab.hpp>
#include <graphlab/rpc/dht.hpp>

inline double pLogP(double pr) {
	return (pr > 0.0 ? pr * log(pr) : 0.0);
}


void findAssignedPart(int* start, int* end, int numAll, int numProc, int myID);

// Global variables...
double RESET_PROB = 0.15;	// random reset probability.
double TOLERANCE = 1.0E-15;
double THRESHOLD = 1.0E-3;
double TOTALWEIGHT = 1.0;	// Should be updated after reading graph.
long N_NODE = 1;			// Should be updated after reading graph. number of nodes.
double SUM_DANGLING_SIZE = 1.0;
double SUM_ALL_SIZE = 1.0;		// sum of all size before normalize, and will be used for normalize node-size.
int MAX_ITER = 10;			// maximum iteration for graph findCommunity vertex program.
int MAX_SP_ITER = 3;		// maximum iteration for sp-graph findCommunity vertex program.
int MAX_ITER_PGRANK = 200;	// maximum iteration for calculating ergodic state (pagerank) of vertices.
int INTERVAL = 3;			// time interval for checking whether movement message valid or not.

int nextEmptyModID = 100000000;		// temporary emptyMod target ID. 

//  This will be a container of the information of each module.
//	This ModuleInfo will be stored in ditributed hash table (dht),
//	so the index will be the key of dht.
struct ModuleInfo : public graphlab::IS_POD_TYPE {
	int index;		// Index number of this module.
						// Index will be a key of dht for accessing corresponding ModuleInfo object.
	double exitPr;		// exit probability of this module.
	double stayPr;		// stay probability of this module, which is sum of p_alpha + exit probability.
	double sumPr;		// sum of p_alpha.  alpha is in this module i.
	double sumTPWeight;	// sum of teleport weight.   SUM (tau_a).
	double sumDangling;	// sum of dangling nodes weight.
	int numMembers;		// number of members of this Module.

	double modOutFlow;	// outFlow from the module. == sum of outFlow to other modules from members.

	//std::vector<int> members;	// This vector will keep the id() of this module's members.

	ModuleInfo(): index(1), exitPr(0.0), stayPr(0.0), sumPr(0.0), sumTPWeight(0.0), sumDangling(0.0), numMembers(0), modOutFlow(0.0) {}
	
	explicit ModuleInfo(int idx, double exitProb, double sumProb, double sumTPW, double sumDang, int numMem):
		index(idx), exitPr(exitProb), sumPr(sumProb), stayPr(exitProb + sumProb),
		sumTPWeight(sumTPW), sumDangling(sumDang), numMembers(numMem), modOutFlow(0.0) {}
};


struct ModuleInfos {
	std::map<int, ModuleInfo> modInfoHash;

	void save(graphlab::oarchive& oarc) const {
		oarc << modInfoHash;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> modInfoHash;
	}

	ModuleInfos& operator+=(const ModuleInfos& other) {
		for (std::map<int, ModuleInfo>::const_iterator it = other.modInfoHash.begin(); it != other.modInfoHash.end(); it++) {
			//if (modInfoHash.count(it->first) > 0) {
			if (modInfoHash.find(it->first) != modInfoHash.end()) {
				modInfoHash[it->first].sumPr += it->second.sumPr;
				modInfoHash[it->first].sumTPWeight += it->second.sumTPWeight;
				modInfoHash[it->first].sumDangling += it->second.sumDangling;
				modInfoHash[it->first].numMembers += it->second.numMembers;
				modInfoHash[it->first].exitPr += it->second.exitPr;
				modInfoHash[it->first].stayPr += it->second.stayPr;
				modInfoHash[it->first].modOutFlow += it->second.modOutFlow;
			}
			else {
				ModuleInfo mInfo;
				
				mInfo.index = it->second.index;
				mInfo.sumPr = it->second.sumPr;
				mInfo.sumTPWeight = it->second.sumTPWeight;
				mInfo.sumDangling = it->second.sumDangling;
				mInfo.numMembers = it->second.numMembers;
				mInfo.exitPr += it->second.exitPr;
				mInfo.stayPr += it->second.stayPr;
				mInfo.modOutFlow += it->second.modOutFlow;

				modInfoHash[it->first] = mInfo;
			}
		}

		return *this;
	}

};


ModuleInfos* modInfoPtr;


class distributed_ModInfos {
  private:
	graphlab::dc_dist_object<distributed_ModInfos> rmi;	// Local RMI object.
	size_t nProcs;
	size_t myID;
	graphlab::mutex lock;
 
  public:
	ModuleInfos localModInfos;							// local storage for EdgeMap.

	distributed_ModInfos(graphlab::distributed_control& dc) : rmi(dc, this) {
		nProcs = rmi.dc().numprocs();
		myID = rmi.dc().procid();
		rmi.barrier();				// make sure all machines finish constructing this object.
	}

	void combineModInfos(size_t procTo, ModuleInfos other) {
		if (procTo == myID) {
			lock.lock();
			localModInfos += other;
			lock.unlock();
		}
		else {
			// this function will call remote_call to the procTo.
			rmi.remote_call(procTo, &distributed_ModInfos::combineModInfos, procTo, other);
		}
	}

	// Must be called by all machine simultaneously.
	void clearLocalModInfos() {
		rmi.barrier();
		//std::map<int, ModuleInfo>().swap(localModInfos.modInfoHash);
		localModInfos.modInfoHash.clear();
		rmi.barrier();
	}
};





/*
 * Module State struct, which is based on gossiping protocol.
 */
struct ModState : public graphlab::IS_POD_TYPE {
	int modID;
	int numMembers;			// the number of members in this module.
	double modPr;			// sum of pagerank of member vertices.
	double sumTPWeight;		// sum of teleport weight of member vertices.
	double sumDangSize;		// sum of pagerank of dangling vertices in the module.

	double modOutFlow;		// out-flow from this module to other modules.
	//double inModFlow;		// in-flow from other modules to this module.

	ModState(): modID(1), numMembers(1), modPr(0.0), sumTPWeight(0.0), sumDangSize(0.0), 
				modOutFlow(0.0) {} 

	//explicit ModState(int id, int nMem, double sumPr, double sumTPW, double sumDang, double outF, double inF):
	explicit ModState(int id, int nMem, double sumPr, double sumTPW, double sumDang, double outF):
			modID(id), numMembers(nMem), modPr(sumPr), sumTPWeight(sumTPW), sumDangSize(sumDang),
			modOutFlow(outF) {}

	explicit ModState(const ModuleInfo modInfo) {
		modID = modInfo.index;
		numMembers = modInfo.numMembers;
		modPr = modInfo.sumPr;
		sumTPWeight = modInfo.sumTPWeight;
		sumDangSize = modInfo.sumDangling;
		modOutFlow = modInfo.modOutFlow;
	}

	// This will be called during the gather phase of findCommunity vertex-program
	// to sum ModState values for averaging them.
	ModState& operator+=(const ModState& other) {
		//this->numMembers += other.numMembers;
		//this->modPr += other.modPr;
		//this->sumTPWeight += other.sumTPWeight;
		//this->sumDangSize += other.sumDangSize;
		//this->modOutFlow += other.modOutFlow;
		if (this->modPr < other.modPr) {
			this->numMembers = other.numMembers;
			this->modPr = other.modPr;
			this->sumTPWeight = other.sumTPWeight;
			this->sumDangSize = other.sumDangSize;
			this->modOutFlow = other.modOutFlow;
		}
		return *this;
	}
};



/*
 *	move message struct, for a vertex (v) from srcMod to dstMod.
 *	NOTE that the move-unit could be superNode later..
 */
struct MoveMsg : public graphlab::IS_POD_TYPE {
	int vid;		// vertex id.
	int atIter;		// the iteration number at the movement of the vertex. (can be used as a timestamp!)
	int srcMod;		// source mod ID.
	int dstMod;		// destination mod ID

	int nVertices;		// the number of vertices in this movement.
	double ndSize;		// the size (pagerank) of the moved vertex.
	double ndTPWeight;	// the teleport weight of the vertex.  (Usually, it is propotional to the number of vertices in this move-unit.
	double ndDangSize;	// the size of dangling vertex.

	double deltaSrcOutFlow;		// the change of outflow of the srcMod.
	double deltaDstOutFlow;		// the change of outflow of the dstMod.
	
	MoveMsg(): vid(1), srcMod(1), dstMod(1), atIter(1), nVertices(1),
				ndSize(1.0), ndTPWeight(1.0), ndDangSize(0.0), deltaSrcOutFlow(0.0), deltaDstOutFlow(0.0) {}
	
	explicit MoveMsg(int id, int src, int dst, int iter, int nv, double size, double tpw, double dSize, double dSrc, double dDst):
				vid(id), srcMod(src), dstMod(dst), atIter(iter), nVertices(nv),
				ndSize(size), ndTPWeight(tpw), ndDangSize(dSize), deltaSrcOutFlow(dSrc), deltaDstOutFlow(dDst) {}

};



typedef std::pair<int, int> msg_key;		// msg_key contains (vid, iter).
typedef std::map<msg_key, MoveMsg> message_map;	// this will be the storage of the messages received from neighbors in node_data.

/*
 *	This will be actual message which is sent by signaling,
 *	so this should implement +operator function, for adding all messages since last run.
 *	Simply, we will have a message_map as a container for combined MoveMsg objects.
 */
struct Messages {
	message_map msgs;

	void save(graphlab::oarchive& oarc) const {
		oarc << msgs;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> msgs;
	}

	// Addition of two Messages objects.
	Messages& operator+=(const Messages& other) {
		// add ONLY MoveMsg objects which are not in msgs container.
		for (message_map::const_iterator it = other.msgs.begin(); it != other.msgs.end(); it++) {
			if (msgs.find(it->first) == msgs.end())
				msgs.insert((*it));
		}

		return *this;
	}
};




// storage for the msg_keys which represent the messages already seen.
typedef std::set<msg_key> received_keys;	// first : vID, second: iter.

// Define vertex data structure.
struct node_data {
	double size;			// p_alpha.
	int nMembers;		// number of members in this node, which is for sp-node case. otherwise, 1.

	double exitPr;		// exitPr for superNode or node...
	double ndOutFlow;	// sum of out-flow of this node...

	double nodeWeight;	// node weight for non-uniform teleport weight..
	int modIdx;				// module index number.
	int prevModId;			// previous module id. It will be used for rolling back to the previous module assignment.

	bool isSuper;			// true: if it is superNode (it is a group of nodes in the current level)
							// false: otherwise

	double teleportWeight;	// teleportWeight = normalized nodeWeight (nodeWeight/totalNodeWeight)
	
	double danglingSize;		// danglingSize will be equivalent to sumDangliing of corresponding module.
	bool isDangling;		// true: if it is a dangling node, false: otherwise.

	int pgIter;			// iteration counter for pagerank calculation.
	int iter;			// iteration counter for finding community vertex program.

	ModState myModStat;		// this will be updated based on received messages.


	void save(graphlab::oarchive& oarc) const {
		oarc << size << nMembers << exitPr << ndOutFlow << nodeWeight << modIdx << prevModId << isSuper << teleportWeight << danglingSize << isDangling << pgIter << iter << myModStat;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> size >> nMembers >> exitPr >> ndOutFlow >> nodeWeight >> modIdx >> prevModId >> isSuper >> teleportWeight >> danglingSize >> isDangling >> pgIter >> iter >> myModStat;
	}
	

	node_data(): nodeWeight(1.0), nMembers(1), isSuper(false), danglingSize(0.0), isDangling(false), pgIter(0), iter(0) {}
	explicit node_data(double ndWeight) : nodeWeight(ndWeight), nMembers(1), isSuper(false), danglingSize(0.0), isDangling(false), pgIter(0), iter(0) {}

	explicit node_data(const ModuleInfo& modInfo) {
		size = modInfo.sumPr;
		nMembers = modInfo.numMembers;
		ndOutFlow = modInfo.modOutFlow;		
		exitPr = modInfo.exitPr;
		nodeWeight = modInfo.sumTPWeight;
		modIdx = modInfo.index;
		prevModId = modIdx;
		isSuper = true;
		teleportWeight = modInfo.sumTPWeight;
		danglingSize = modInfo.sumDangling;
		isDangling = size == danglingSize ? true : false;
		pgIter = iter = 0;

		myModStat = ModState(modInfo);
	}
};



typedef std::map<int, double> flowmap;		// <modID, flow>
typedef std::map<int, double> edgeFlow;		// <vID, edgeFlow>

// Define flow_type data structure as the gather_type of the core-algorithm...
struct flow_type {
	flowmap outFlowToMod;	// outFlowToMod[modID] = outflow
	flowmap inFlowFromMod;	// inFlowFromMod[modID] = inflow

	std::map<int, ModState> neighborModStates;	// <modID, ModState> 
	std::map<int, int> numModOccurs;			// <modID, #occurrence>

	flow_type() {}

	void save(graphlab::oarchive& oarc) const {
		oarc << outFlowToMod << inFlowFromMod << neighborModStates << numModOccurs;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> outFlowToMod >> inFlowFromMod >> neighborModStates >> numModOccurs;
	}

	flow_type& operator+=(const flow_type& other) {
		for (flowmap::const_iterator it = other.outFlowToMod.begin(); it != other.outFlowToMod.end(); it++) {
			if(outFlowToMod.find(it->first) != outFlowToMod.end()) {
				outFlowToMod[it->first] += it->second;
			}
			else {
				outFlowToMod[it->first] = it->second;
			}
		}

		for (flowmap::const_iterator it = other.inFlowFromMod.begin(); it != other.inFlowFromMod.end(); it++) {
			if(inFlowFromMod.find(it->first) != inFlowFromMod.end()) {
				inFlowFromMod[it->first] += it->second;
			}
			else {
				inFlowFromMod[it->first] = it->second;
			}
		}

		for (std::map<int, ModState>::const_iterator it = other.neighborModStates.begin(); it != other.neighborModStates.end(); it++) {
			if (neighborModStates.find(it->first) != neighborModStates.end()) {
				neighborModStates[it->first] += it->second;
				numModOccurs[it->first] += other.numModOccurs.find(it->first)->second;
			}
			else {
				neighborModStates[it->first] = it->second;
				numModOccurs[it->first] = other.numModOccurs.find(it->first)->second;
			}
		}

		return *this;
	}
};


typedef node_data vertex_data_type;

typedef double edge_data_type;	// edge-weight is an edge data (double).

// the graph type is determined by the vertex and edge data types.
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;


struct OutFlows {
	std::map<int, double> outFlow;

	void save(graphlab::oarchive& oarc) const {
		oarc << outFlow;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> outFlow;
	}

	void addEdge(const graph_type::edge_type& e) {
		int src = e.source().data().modIdx;
		int dst = e.target().data().modIdx;

		if (src != dst) {
			if (outFlow.find(src) != outFlow.end())
				outFlow[src] += e.data();
			else 
				outFlow[src] = e.data();

			if (outFlow.find(dst) == outFlow.end())
				outFlow[dst] = 0.0;
		}
		else {
			if (outFlow.find(src) == outFlow.end())
				outFlow[src] = 0.0;
		}
	}

	OutFlows& operator+=(const OutFlows& other) {
		for (std::map<int, double>::const_iterator it = other.outFlow.begin(); it != other.outFlow.end(); it++) {
			if (outFlow.find(it->first) != outFlow.end()) {
				outFlow[it->first] += it->second;
			}
			else {
				outFlow[it->first] = it->second;
			}
		}

		return *this;
	}
};



class distributed_OutFlows {
  private:
	graphlab::dc_dist_object<distributed_OutFlows> rmi;	// Local RMI object.
	size_t nProcs;
	size_t myID;
	graphlab::mutex lock;
 
  public:
	OutFlows localOutFlows;							// local storage for EdgeMap.

	distributed_OutFlows(graphlab::distributed_control& dc) : rmi(dc, this) {
		nProcs = rmi.dc().numprocs();
		myID = rmi.dc().procid();
		rmi.barrier();				// make sure all machines finish constructing this object.
	}

	void combineOutFlows(size_t procTo, OutFlows other) {
		if (procTo == myID) {
			lock.lock();
			localOutFlows += other;
			lock.unlock();
		}
		else {
			// this function will call remote_call to the procTo.
			rmi.remote_call(procTo, &distributed_OutFlows::combineOutFlows, procTo, other);
		}
	}

	// Must be called by all machine simultaneously.
	void clearLocalOutFlowsMap() {
		rmi.barrier();
		//std::map<int, double>().swap(localOutFlows.outFlow);
		localOutFlows.outFlow.clear();
		rmi.barrier();
	}
};




/*
 * A simple function used by graph.transform_vertices(init_vertex);
 * to initialize the vertex data.
 */
void init_vertex(graph_type::vertex_type& v) { 
	v.data().modIdx = (int) v.id(); 
	v.data().size = 1.0 / N_NODE;
	v.data().teleportWeight = v.data().nodeWeight / TOTALWEIGHT;
	if (v.num_out_edges() == 0) {
		v.data().isDangling = true;
		v.data().danglingSize = v.data().size;
	}
}

void reset_modInfo(graph_type::vertex_type& v) { 
	v.data().modIdx = (int) v.id(); 
	v.data().iter = 0;

	v.data().myModStat = ModState((int)v.id(), v.data().nMembers, v.data().size, v.data().teleportWeight, v.data().danglingSize, v.data().ndOutFlow);
}


bool is_dangling(const graph_type::vertex_type& v) {
	return v.data().isDangling;
}

// This map-reduce function called for only dangling-vertex set.
double dangling_vertex_size(const graph_type::vertex_type& vertex) {
	//double dSize = 0.0;
	//if (vertex.data().isDangling)
	//	dSize = vertex.data().size;
	
	//return dSize;
	return  vertex.data().size;		// This is the case when restricted the dangling-vertex set.
}

double sum_size(const graph_type::vertex_type& vertex) {
	return vertex.data().size;
}

void normalize_vertex_size(graph_type::vertex_type& vertex) {
	vertex.data().size /= SUM_ALL_SIZE;
}


// This should be updated for general cases, but for SNAP dataset, this will work correctly.
void normalize_edge_value(graph_type::edge_type& e) {
	e.data() = 1.0 / e.source().num_out_edges();	// e.weight / sumEdgeWeights ....
}


void update_edge_flow(graph_type::edge_type& edge) {
	edge.data() *= edge.source().data().size;		// p_a * e.weight  <-- for faster running...
}


/********************
 * map-functions for the related map-reduce functions...
 ********************/
double sum_plogp_vertex(const graph_type::vertex_type& vertex) {
	return pLogP(vertex.data().size);
}

double sumExitPr(const graph_type::vertex_type& v) {
	return v.data().exitPr;
}

double sum_exitLogExit(const graph_type::vertex_type& v) {
	return pLogP(v.data().exitPr);
}

double sum_stayLogStay(const graph_type::vertex_type& v) {
	return pLogP(v.data().exitPr + v.data().size);
}


void addVertexToModInfos(ModuleInfos& mInfos, const graph_type::vertex_type& v) {
	int modIdx = v.data().modIdx;

	if (mInfos.modInfoHash.find(modIdx) != mInfos.modInfoHash.end()) {
		ModuleInfo& mInfoRef = mInfos.modInfoHash[modIdx];
	
		mInfoRef.sumPr += v.data().size;
		mInfoRef.sumTPWeight += v.data().teleportWeight;
		mInfoRef.sumDangling += v.data().danglingSize;
		mInfoRef.numMembers += v.data().nMembers;
	}
	else {
		ModuleInfo mInfo;

		mInfo.sumPr = v.data().size;
		mInfo.index = v.data().modIdx;
		mInfo.sumTPWeight = v.data().teleportWeight;
		mInfo.sumDangling = v.data().danglingSize;
		mInfo.numMembers = v.data().nMembers;

		mInfos.modInfoHash[modIdx] = mInfo;
	}
}



// update mod-state after converged.
void updateModState(graph_type::vertex_type& v) {
	ModState& myModRef = v.data().myModStat;
	int mid = v.data().modIdx;

	v.data().prevModId = v.data().modIdx;

	myModRef.modID = mid;
	myModRef.numMembers = modInfoPtr->modInfoHash[mid].numMembers;
	myModRef.modPr = modInfoPtr->modInfoHash[mid].sumPr;
	myModRef.sumTPWeight = modInfoPtr->modInfoHash[mid].sumTPWeight;
	myModRef.sumDangSize = modInfoPtr->modInfoHash[mid].sumDangling;
	myModRef.modOutFlow = modInfoPtr->modInfoHash[mid].modOutFlow;

	v.data().iter = 0;	// iteration number should be reset for running next loop.
}

// roll-back module assignment when MDL is increased...
void rollBackModID(graph_type::vertex_type& v) {
	v.data().modIdx = v.data().prevModId;
}


struct ModIDs {
	std::map<int, int> modIDs;		// <vid, mod_id>

	void save(graphlab::oarchive& oarc) const {
		oarc << modIDs;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> modIDs;
	}

	ModIDs& operator+=(const ModIDs other) {
		for (std::map<int, int>::const_iterator it = other.modIDs.begin(); it != other.modIDs.end(); it++) {
			// This would concatenate two different maps of <vid, mod_id>s, and each vid will appear only once,
			// so don't need to check existence test.
			modIDs.insert((*it));
		}

		return *this;
	}
};

ModIDs* modAssignPtr;	// used for hold global pointer for ModIDs to use within transform_vertex() functor.

ModIDs aggregateModAssign(const graph_type::vertex_type& v) {
	ModIDs mID;

	mID.modIDs[(int)v.id()] = v.data().modIdx;

	return mID;
}


void updateModAssign(graph_type::vertex_type& v) {
	// The previous modIdx is a vertex id of corresponding sp_node.
	int sp_id = v.data().modIdx;
	v.data().modIdx = modAssignPtr->modIDs[sp_id];
}



typedef std::map<std::pair<int, int>, edge_data_type> edgeList;		//map<<src, dst>, edgeflow>

struct EdgeMap {
	edgeList eList;		//edge list.

	void save(graphlab::oarchive& oarc) const {
		oarc << eList;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> eList;
	}

	void add(const graph_type::local_edge_type& e) {
		int src = e.source().data().modIdx;
		int dst = e.target().data().modIdx;

		if (src != dst) {
			if (eList.find(std::make_pair(src, dst)) != eList.end()) {
				eList[std::make_pair(src, dst)] += e.data();
			}
			else {
				eList[std::make_pair(src, dst)] = e.data();
			}
		}
	}

	EdgeMap& operator+=(const EdgeMap& other) {
		for (edgeList::const_iterator it = other.eList.begin(); it != other.eList.end(); it++) {
			if (eList.find(it->first) != eList.end()) {
				eList[it->first] += it->second;
			}
			else {
				eList[it->first] = it->second;
			}
		}

		return *this;
	}
};


class global_EdgeMap {
  private:
	graphlab::dc_dist_object<global_EdgeMap> rmi;	// Local RMI object.
	size_t nProcs;
	size_t myID;
	graphlab::mutex lock;
 
  public:
	EdgeMap localEdgeMap;							// local storage for EdgeMap.

	global_EdgeMap(graphlab::distributed_control& dc) : rmi(dc, this) {
		nProcs = rmi.dc().numprocs();
		myID = rmi.dc().procid();
		rmi.barrier();				// make sure all machines finish constructing this object.
	}

	void combineEdgeMap(size_t procTo, EdgeMap other) {
		if (procTo == myID) {
			lock.lock();
			localEdgeMap += other;
			lock.unlock();
		}
		else {
			// this function will call remote_call to the procTo.
			rmi.remote_call(procTo, &global_EdgeMap::combineEdgeMap, procTo, other);
		}
	}

	// Must be called by all machine simultaneously.
	void clearLocalEdgeMap() {
		rmi.barrier();
		//edgeList().swap(localEdgeMap.eList);
		localEdgeMap.eList.clear();
		rmi.barrier();
	}
};


// This is the function of calculating ergodic node size (probability).
// The ergodic node probability is practically the same as the normalized PageRank value of each node.
class calculateErgodicNodeSize :
	public graphlab::ivertex_program<graph_type, double>,
	public graphlab::IS_POD_TYPE {

	//double sumDanglingSize = 0.0;
	double last_size;

  public:
	// gather: we are going to gather information from in-edges
	edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
		return graphlab::IN_EDGES;
	}

	// Gather the weighted rank of the in-edge adjacent page .
	double gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
		return edge.source().data().size * edge.data();
	}

	// Use the total rank of the adjacent pages from gather() phase to update this page.
	void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
		double newval = (1.0 - RESET_PROB) * total;
		newval += (RESET_PROB + (1.0 - RESET_PROB) * SUM_DANGLING_SIZE) * vertex.data().teleportWeight;
		last_size = vertex.data().size;
		vertex.data().size = newval;
		vertex.data().pgIter++;
	}

	// The scatter edges depends on the last_change value..
	edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
		//if (std::fabs(last_size - vertex.data().size) > TOLERANCE && context.iteration() < MAX_ITER_PGRANK)
		if (std::fabs(last_size - vertex.data().size) > TOLERANCE && vertex.data().pgIter < MAX_ITER_PGRANK)
			return graphlab::OUT_EDGES;
		else
			return graphlab::NO_EDGES;
	}

	// the scatter function just signal adjacent pages.
	void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
		context.signal(edge.target());
	}
};





// This vertex program calculates initial exit-probability of each vertex,
// after calculating the ergodic probability of each vertex.
class update_exitPr :
	public graphlab::ivertex_program<graph_type, double>,
	public graphlab::IS_POD_TYPE {

  public:
	edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
		return graphlab::OUT_EDGES;
	}

	double gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
		return edge.data();
	}
	
	void apply(icontext_type& context, vertex_type& v, const gather_type& sumExitFlow) {
		if (!v.data().isDangling) {
			v.data().exitPr = RESET_PROB * (1.0 - v.data().teleportWeight) * v.data().size
								+ (1.0 - RESET_PROB) * sumExitFlow;
			v.data().ndOutFlow = sumExitFlow;
		} 
		else {
			v.data().exitPr = (1.0 - v.data().teleportWeight) * v.data().size;
			v.data().danglingSize = v.data().size;
			v.data().ndOutFlow = 0.0;		// dangling vertex, no direct out-flow.
		}	

		// initialize myModStat in v.data()
		//ModState(modID, nMember,sumPr,sumTPW, sumDang, outF);
		v.data().myModStat = ModState((int)v.id(), 1, v.data().size, v.data().teleportWeight, v.data().danglingSize, v.data().ndOutFlow);
	}

	edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
		return graphlab::NO_EDGES;
	}

	void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
	}
};




/*
 *	This is the vertex-program for the core algorithm of GossipMap.
 */
class findCommunity :
	public graphlab::ivertex_program<graph_type, flow_type, Messages> {

	bool moved;		// if this vertex moved, it will set to true.  Otherwise, false.
	Messages myMessages;	// This will contain messages sent to neighbors.
	Messages oneMsg;		// This will contain only the message about current movement.

	int myOldMod;
	MoveMsg myMoveMsg;
	msg_key mykey;

  public:
	void save(graphlab::oarchive& oarc) const {
		oarc << moved << myMessages << myOldMod << myMoveMsg << mykey;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> moved >> myMessages >> myOldMod >> myMoveMsg >> mykey;
	}

	// init() : will receive messages from signaling neighbors, and concatenating them.
	void init(icontext_type& context, const vertex_type& v, const message_type& msg) {

		//message_map().swap(myMessages.msgs);	// initialize myMessages container.
		myMessages.msgs.clear();

		if (v.data().iter >= (v.data().isSuper ? MAX_SP_ITER : MAX_ITER))
			return;

		// check the message was already seen or not.
		for (message_map::const_iterator it = msg.msgs.begin(); it != msg.msgs.end(); it++) {
			msg_key key = it->first;
			MoveMsg mvMsg = it->second;

			// first look at the message is related to myMod and is not from vertex itself.
			if ((mvMsg.srcMod == v.data().modIdx || mvMsg.dstMod == v.data().modIdx) && key.first != v.id() && (mvMsg.atIter + INTERVAL > v.data().iter)) {
				myMessages.msgs[key] = mvMsg;
			}
		}
	}


	edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
		return (vertex.data().iter >= (vertex.data().isSuper? MAX_SP_ITER : MAX_ITER)) ? graphlab::NO_EDGES : graphlab::ALL_EDGES;
	}

	flow_type gather(icontext_type& context, const vertex_type& v, edge_type& edge) const {
		// Calculate flow and return it.  If edge == OUT_EDGE, then it will be out-flow. Otherwise, in-flow.
		flow_type flow;

		if( edge.source().id() == v.id() ) {
			// OUT_EDGE case
			int modID = edge.target().data().modIdx;
			flow.outFlowToMod[modID] = (1.0 - RESET_PROB) * edge.data();
			flow.inFlowFromMod[modID] = 0.0;	// since the edge related to out-flow.
			
			ModState& nghbrStateRef = edge.target().data().myModStat;
			ModState cpState;
			cpState.modID = modID;
			cpState.numMembers = nghbrStateRef.numMembers;
			cpState.modPr = nghbrStateRef.modPr;
			cpState.sumTPWeight = nghbrStateRef.sumTPWeight;
			cpState.sumDangSize = nghbrStateRef.sumDangSize;
			cpState.modOutFlow = nghbrStateRef.modOutFlow;

			flow.neighborModStates[modID] = cpState;
			flow.numModOccurs[modID] = 1;
		}
		else {
			// IN_EDGE case
			int modID = edge.source().data().modIdx;
			flow.inFlowFromMod[modID] = (1.0 - RESET_PROB) * edge.data();
			flow.outFlowToMod[modID] = 0.0;		// since the edge related to in-flow.

			ModState& nghbrStateRef = edge.source().data().myModStat;
			ModState cpState;
			cpState.modID = modID;
			cpState.numMembers = nghbrStateRef.numMembers;
			cpState.modPr = nghbrStateRef.modPr;
			cpState.sumTPWeight = nghbrStateRef.sumTPWeight;
			cpState.sumDangSize = nghbrStateRef.sumDangSize;
			cpState.modOutFlow = nghbrStateRef.modOutFlow;

			flow.neighborModStates[modID] = cpState;
			flow.numModOccurs[modID] = 1;
		}

		return flow;
	}

	// Based on in-/out-flow information, and other module information,
	// this method will find a new module which could improve the quality of communities.
	void apply(icontext_type& context, vertex_type& v, const gather_type& inOutFlows) {
		moved = false;		// Initially set as false, and it will be set 'true' if this vertex moves.

		if (v.data().iter >= (v.data().isSuper? MAX_SP_ITER : MAX_ITER))
			return;

		flowmap outFlowToMod = inOutFlows.outFlowToMod;		// copy from gather result.
		flowmap inFlowFromMod = inOutFlows.inFlowFromMod;	// copy from gather result.
		//flowmap moduleFlow = inOutFlows.modFlow;

		// 1st strategy: average neighbor ModStates
		std::map<int, ModState> nghbrModStates = inOutFlows.neighborModStates;
		std::map<int, int> numOccurs = inOutFlows.numModOccurs;

		const int emptyTarget = context.num_vertices() + 1;	// This will be an indicator of moving to emptyModule.

		////////////////////////////////////////////////
		// update messages from signaling neighbors.  //
		////////////////////////////////////////////////
		ModState& myState = v.data().myModStat;

		for (message_map::iterator it = myMessages.msgs.begin(); it != myMessages.msgs.end(); it++) {
			//case 1 - myMod = it->srcMod : removing the vertex from myModStat.
			if ((it->second).srcMod == myState.modID) {
				myState.numMembers -= (it->second).nVertices;
				myState.modPr -= (it->second).ndSize;
				myState.sumTPWeight -= (it->second).ndTPWeight;
				myState.sumDangSize -= (it->second).ndDangSize;
				myState.modOutFlow += (it->second).deltaSrcOutFlow;
			}
			else if ((it->second).dstMod == myState.modID) {
				myState.numMembers += (it->second).nVertices;
				myState.modPr += (it->second).ndSize;
				myState.sumTPWeight += (it->second).ndTPWeight;
				myState.sumDangSize += (it->second).ndDangSize;
				myState.modOutFlow += (it->second).deltaDstOutFlow;
			}
		}


		// copy node specific values
		double ndSize = v.data().size;
		double ndExitPr = v.data().exitPr;
		double ndTPWeight = v.data().teleportWeight;
		double ndDanglingSize = v.data().danglingSize;

		int srcMod = v.data().modIdx;
		double srcModFlow = (outFlowToMod.find(srcMod) != outFlowToMod.end()) ? outFlowToMod[srcMod] + inFlowFromMod[srcMod] : 0.0;

		// add information of the current vertex 'v' here.
		if (nghbrModStates.find(srcMod) != nghbrModStates.end()) {
			nghbrModStates[srcMod] += v.data().myModStat;
			numOccurs[srcMod]++;
		}
		else {
			nghbrModStates[srcMod] = v.data().myModStat;
			numOccurs[srcMod] = 1;
		}

		ModState& srcModState = nghbrModStates[srcMod];
		//if (numOccurs[srcMod] > 1) {
		//	int nOcc = numOccurs[srcMod];
		//	srcModState.numMembers /= nOcc;
		//	srcModState.modPr /= nOcc;
		//	srcModState.sumTPWeight /= nOcc;
		//	srcModState.sumDangSize /= nOcc;
		//	srcModState.modOutFlow /= nOcc;
		//}


		// update myModState w.r.t. neighbors info..
		myState.numMembers = srcModState.numMembers;
		myState.modPr = srcModState.modPr;
		myState.sumTPWeight = srcModState.sumTPWeight;
		myState.sumDangSize = srcModState.sumDangSize;
		myState.modOutFlow = srcModState.modOutFlow;


		double alpha = RESET_PROB;
		double beta = 1.0 - RESET_PROB;

		// add approximated teleportation flow for srcModule...
		// if srcModFlow == 0, it means that no neighbors are in the same module.
		// And we assume that is the case of only one vertex is in the srcMod --> teleportation flow == 0.0.
		if (srcModFlow > 0) {
			srcModFlow += (alpha * ndSize + beta * ndDanglingSize) * (srcModState.sumTPWeight - ndTPWeight);
			srcModFlow += (alpha * (srcModState.modPr - ndSize) + beta * (srcModState.sumDangSize - ndDanglingSize)) * ndTPWeight;
		}

		/////////// Calculate oldExitPr w.r.t. neighbors' information ////////////
		double oldExitPrSrcMod = ndExitPr;

		if (srcModState.numMembers > 1) {
			oldExitPrSrcMod = alpha * (1.0 - srcModState.sumTPWeight) * srcModState.modPr 
								+ beta * (srcModState.modOutFlow + (1.0 - srcModState.sumTPWeight) * srcModState.sumDangSize);
		}

		double newExitPrSrcMod = oldExitPrSrcMod - ndExitPr + srcModFlow;

		ModState emptyModState;

		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (srcModState.numMembers > 1) {	//&& emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;

			emptyModState.numMembers = 0;
			emptyModState.modID = emptyTarget;

			numOccurs[emptyTarget] = 0;
		}

		std::pair<int, double> bestModFlow;		// <bestModIdx, modFlow>
		bestModFlow.second = 0.0;

		int bestMod = srcMod;
		double bestDiffCodeLen = 0.0;

		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) {
			int dstMod = it->first;
			double dstModFlow = it->second + inFlowFromMod[it->first];
			
			if (dstMod != srcMod) {
				ModState& dstModState = nghbrModStates[dstMod];
				//if (numOccurs[dstMod] > 1) {
				//	int nOcc = numOccurs[dstMod];
				//	dstModState.numMembers /= nOcc;
				//	dstModState.modPr /= nOcc;
				//	dstModState.sumTPWeight /= nOcc;
				//	dstModState.sumDangSize /= nOcc;
				//	dstModState.modOutFlow /= nOcc;
				//}

				// add approximated teleportation flow for dstModule...
				dstModFlow += (alpha * ndSize + beta * ndDanglingSize) * dstModState.sumTPWeight;
				dstModFlow += (alpha * dstModState.modPr + beta * dstModState.sumDangSize) * ndTPWeight;

				double oldExitPrDstMod = alpha * (1.0 - dstModState.sumTPWeight) * dstModState.modPr 
									+ beta * (dstModState.modOutFlow + (1.0 - dstModState.sumTPWeight) * dstModState.sumDangSize);
			
				double newExitPrDstMod = oldExitPrDstMod + ndExitPr - dstModFlow;

				double delta_exit_log_exit = pLogP(newExitPrSrcMod) + pLogP(newExitPrDstMod) - pLogP(oldExitPrSrcMod) - pLogP(oldExitPrDstMod);
				double delta_stay_log_stay = pLogP(newExitPrSrcMod + srcModState.modPr - ndSize) + pLogP(newExitPrDstMod + dstModState.modPr + ndSize) \
											 - pLogP(oldExitPrSrcMod + srcModState.modPr) - pLogP(oldExitPrDstMod + dstModState.modPr);

				double diffCodeLen = delta_stay_log_stay - 2.0 * delta_exit_log_exit;

				if (diffCodeLen < bestDiffCodeLen) {
					bestMod = dstMod;
					bestDiffCodeLen = diffCodeLen;
				}
				
			}
		}

		v.data().iter++;

		myOldMod = srcMod;

		if (bestMod != srcMod) {
			v.data().modIdx = bestMod;

			if (bestMod == emptyTarget) {
				v.data().modIdx = nextEmptyModID++;
			}

			double dSrcOutFlow = inFlowFromMod[srcMod] + outFlowToMod[srcMod] - v.data().ndOutFlow;
			double dDstOutFlow = v.data().ndOutFlow - outFlowToMod[bestMod] - inFlowFromMod[bestMod];

			myMoveMsg = MoveMsg(v.id(), srcMod, v.data().modIdx, v.data().iter, v.data().nMembers, ndSize, 
								ndTPWeight, ndDanglingSize, dSrcOutFlow, dDstOutFlow);
			mykey = std::make_pair(v.id(), v.data().iter);	// msg_key for myMsgs.

			myMessages.msgs[mykey] = myMoveMsg;	// inserting my move message into myMsgs.
			oneMsg.msgs[mykey] = myMoveMsg;

			// update myModStat
			ModState& bestStat = nghbrModStates[bestMod];

			myState.modID = v.data().modIdx;
			myState.numMembers = bestStat.numMembers + v.data().nMembers;
			myState.modPr = bestStat.modPr + ndSize;
			myState.sumTPWeight = bestStat.sumTPWeight + ndTPWeight;
			myState.sumDangSize = bestStat.sumDangSize + ndDanglingSize;
			myState.modOutFlow = bestStat.modOutFlow + dDstOutFlow;

			moved = true;
		}

	}	// END OF apply()


	// scatter-edges
	edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
		return (moved) ? graphlab::ALL_EDGES : graphlab::NO_EDGES;
	}

	// scatter
	void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
		if( edge.source().id() == vertex.id() ) {
			// out-edges.
			if (edge.target().data().iter < (edge.target().data().isSuper ? MAX_SP_ITER : MAX_ITER)) {
				context.signal(edge.target(), oneMsg);
			}
		}
		else {
			// in-edges.
			if (edge.source().data().iter < (edge.source().data().isSuper ? MAX_SP_ITER : MAX_ITER)) {
				context.signal(edge.source(), oneMsg);
			}
		}
	}
};






/*
 * We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", relaxmap_writer()) to save the graph.
 *	COPY FROM demoapps/pagerank/simple_pagerank.cpp.
 */
struct relaxmap_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data().modIdx << "\t" << v.data().size << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of relaxmap writer




// function declaration...
void getDistModInfos(distributed_ModInfos& distModInfos, graph_type& graph, graphlab::distributed_control& dc);
void getDistOutFlows(distributed_OutFlows& distOutFlows, graph_type& graph, graphlab::distributed_control& dc);
double calculateExitProb(const std::map<int, double>& outFlow, std::map<int, ModuleInfo>& modInfo, double& new_sum_log_exit, double& new_sum_log_stay);



int main(int argc, char** argv) {
	// Initialize control plain using mpi...
	graphlab::mpi_tools::init(argc, argv);
	graphlab::distributed_control dc;
	global_logger().set_log_level(LOG_INFO);

	// Parse command line options ---
	graphlab::command_line_options clopts("GossipMap Algorithm.");

	std::string graph_file;				// graph input file name.
	std::string format = "snap";		// graph file format.
	std::string exec_type = "async";	// engine running type: synchronous or asynchronous.
	int trials = 1;						// the number of trials for this test.
	int mode = 1;						// the running mode flag. 1 - coreOnce, 2 - coreRepeat.  Default = coreOnce.
	bool repeatMode = false;				// repeat mode is unset in default.
	int outerMode = 2;						// the running outerloop mode flag. 1 - outerOnce, 2 - outerRepeat.  Default = outerRepeat.
	bool repeatOuterMode = true;				// repeat mode is set in default.


	clopts.attach_option("graph", graph_file, "The graph file. Required.");
	clopts.add_positional("graph");
	clopts.attach_option("format", format, "The graph file format.");
	//clopts.attach_option("sync", exec_type, "The engine running type: synchronous or asynchronous.");
	clopts.attach_option("thresh", THRESHOLD, "The threshold for convergence condition.");
	clopts.attach_option("tol", TOLERANCE, "The threshold for pagerank (ergodic state) convergence condition.");
	clopts.attach_option("maxiter", MAX_ITER, "The maximum of the iteration for finding community.");
	clopts.attach_option("maxspiter", MAX_SP_ITER, "The maximum of the iteration of sp-graph for finding community.");
	clopts.attach_option("trials", trials, "The number of trials for finding community repeatedly.");
	clopts.attach_option("interval", INTERVAL, "The time interval for checking whether the received message is valid or not.");
	clopts.attach_option("mode", mode, "The running mode of finding community: 1 - coreOnce, 2 - coreRepeat.");
	clopts.attach_option("outmode", outerMode, "The running outerloop mode of finding community: 1 - outerOnce, 2 - outerRepeat.");

	std::string prefix;
	clopts.attach_option("prefix", prefix, "If set, this app will save the community detection result.");

	if (!clopts.parse(argc, argv)) {
		dc.cout() << "Error in parsing command line arguments." << std::endl;
		clopts.print_description();
		return EXIT_FAILURE;
	}

	if (graph_file == "") {
		dc.cout() << "No graph file specified. Cannot continue." << std::endl;
		clopts.print_description();
		return EXIT_FAILURE;
	}
	
	// Build the graph
	graph_type graph(dc, clopts);
	dc.cout() << "Loading graph in format: " << format << std::endl;
	graph.load_format(graph_file, format);
	graph.finalize();	// We should call finalize() function before querying the graph.

	dc.cout() << "# vertices: " << graph.num_vertices()
			  << ", # edges: " << graph.num_edges() << std::endl;

	N_NODE = graph.num_vertices();
	
	if (format == "snap")
		TOTALWEIGHT = (double) graph.num_vertices();
	else
		TOTALWEIGHT = (double) graph.num_vertices();		// TODO: need to update.

	if (mode == 2) {
		dc.cout() << "mode == 2, so set repeatMode = true." << std::endl;
		repeatMode = true;
	}

	if (outerMode == 1) {
		dc.cout() << "outerMode == 1, so set repeatOuterMode = false." << std::endl;
		repeatOuterMode = false;
	}

		
	// Initialize the vertex data..
	graph.transform_vertices(init_vertex);
	graphlab::vertex_set danglingNodes = graph.select(is_dangling);
	int numDangling = graph.vertex_set_size(danglingNodes);

	// Calculate initial sum-dangling-size ....
	SUM_DANGLING_SIZE = graph.map_reduce_vertices<double>(dangling_vertex_size, danglingNodes);

	dc.cout() << "The number of dangling nodes = " << numDangling << std::endl;
	dc.cout() << "Initial sum of dangling size = " << SUM_DANGLING_SIZE << std::endl;


	graph.transform_edges(normalize_edge_value, graph.complete_set(), graphlab::OUT_EDGES);

	// Running the engine...
	std::string pg_exec_type = "async";
	graphlab::omni_engine<calculateErgodicNodeSize> engine(dc, graph, pg_exec_type, clopts);
	engine.signal_all();
	engine.start();

	SUM_ALL_SIZE = graph.map_reduce_vertices<double>(sum_size);
	graph.transform_vertices(normalize_vertex_size);

	const double runtime = engine.elapsed_seconds();
	dc.cout() << "Runnining time of Calculating ergodic state probability for GLab_infomap: " << runtime << " seconds." << std::endl;




	graph.transform_edges(update_edge_flow, graph.complete_set(), graphlab::OUT_EDGES);	// edge.data() will store p_a * w_ab for reducing repeated computation.

	graphlab::omni_engine<update_exitPr> engine2(dc, graph, exec_type, clopts);
	engine2.signal_all();
	engine2.start();
	const double initCodeLengthTime = engine2.elapsed_seconds();
	dc.cout() << "Running time of calculating initial exit probability of the input graph: " << initCodeLengthTime << " seconds." << std::endl;



	double sum_size_log_size = graph.map_reduce_vertices<double>(sum_plogp_vertex);
	double sumAllExit = graph.map_reduce_vertices<double>(sumExitPr);
	double sumExit_log_sumExit = pLogP(sumAllExit);
	double sum_exit_log_exit = graph.map_reduce_vertices<double>(sum_exitLogExit);
	double sum_stay_log_stay = graph.map_reduce_vertices<double>(sum_stayLogStay);
	double initCodeLength = (sumExit_log_sumExit - 2.0 * sum_exit_log_exit + sum_stay_log_stay - sum_size_log_size) / log(2.0);

	dc.cout() << "sumAllExitPr = " << sumAllExit << std::endl;
	dc.cout() << "sumExit_log_sumExit = " << sumExit_log_sumExit / log(2.0) << std::endl;
	dc.cout() << "- sum_size_log_size = " << - (sum_size_log_size / log(2.0)) << std::endl;
	dc.cout() << "- 2.0 * sum_exit_log_exit = " << - 2.0 * (sum_exit_log_exit / log(2.0)) << std::endl;
	dc.cout() << "sum_stay_log_stay = " << sum_stay_log_stay / log(2.0) << std::endl;

	dc.cout() << "Initial Code Length = " << initCodeLength << std::endl;

	
	// initialize nextEmptyModID
	nextEmptyModID += 1000000 * dc.procid();

	double outer_thresh = THRESHOLD * 3.0;


	int numProcs = dc.numprocs();

	dc.full_barrier();

	/////////////////////////////////////////
	// Now, ready to run a vertex program  //
	// for finding new community of each   //
	// node...							   //
	/////////////////////////////////////////


	for (int i = 1; i <= trials; i++) {
		dc.cout() << "Trial: " << i << std::endl;

		graph.transform_vertices(reset_modInfo);		// This will set v.data().modIdx = v.id() and reset myModStat, correspondingly.

		dc.full_barrier();
	
		graphlab::timer timer;
		timer.start();			// This starts the time for measuring elapsed time...

		double currentMDL = initCodeLength;
		double prevMDL = initCodeLength;
		double outerPrevMDL = initCodeLength;

		int numModule = 0;
	
		double new_sum_log_exit = 0.0;
		double new_sum_log_stay = 0.0;
		double sumAllExitPr = 0.0;		

		double sumFindCommTime = 0.0;

		int superStep = 1;

		graphlab::timer timer2;		// This timer will be used for measuring sub-routine times...

		double exitCalcTime = 0.0;
		double getModInfoTime = 0.0;
		double combineOutFlowsTime = 0.0;
		double aggrEdgeMapTime = 0.0;
		double spGraphGenTime = 0.0;
		double spGraphFinalizeTime = 0.0;

		ModuleInfos modInfos;
		ModuleInfos oldModInfos;

		distributed_ModInfos distModInfos(dc);
		distributed_OutFlows distOutFlows(dc);

		global_EdgeMap distEdgeMap(dc);

		bool swapModInfos = false;

		do {
			outerPrevMDL = currentMDL;

			do {
				prevMDL = currentMDL;

				numModule = 0;
	
				new_sum_log_exit = 0.0;
				new_sum_log_stay = 0.0;
				sumAllExitPr = 0.0;		

				if (swapModInfos) {
					oldModInfos.modInfoHash.swap(modInfos.modInfoHash);
					modInfos.modInfoHash.clear();
				}


				graphlab::omni_engine<findCommunity> engine4(dc, graph, exec_type, clopts);
				engine4.signal_all();
				engine4.start();


				const double findCommTime = engine4.elapsed_seconds();
				sumFindCommTime += findCommTime;
				dc.cout() << "Running time of finding communities of the input graph: " << findCommTime << " seconds." << std::endl;

				dc.full_barrier();

				timer2.start();

				distModInfos.clearLocalModInfos();

				getDistModInfos(distModInfos, graph, dc);

				getModInfoTime += timer2.current_time();

				dc.full_barrier();

				timer2.start();


				distOutFlows.clearLocalOutFlowsMap();

				getDistOutFlows(distOutFlows, graph, dc);

				combineOutFlowsTime += timer2.current_time();


				dc.full_barrier();

				timer2.start();

				
				sumAllExitPr = calculateExitProb(distOutFlows.localOutFlows.outFlow, distModInfos.localModInfos.modInfoHash, new_sum_log_exit, new_sum_log_stay);

				dc.full_barrier();

				exitCalcTime += timer2.current_time();

				// Here we have only partial modInfos and OutFlows...
				// So, we need to combine them all here!
				modInfos.modInfoHash.swap(distModInfos.localModInfos.modInfoHash);
				dc.all_reduce(modInfos);
				dc.all_reduce(sumAllExitPr);
				dc.all_reduce(new_sum_log_exit);
				dc.all_reduce(new_sum_log_stay);

				modInfoPtr = &modInfos;

				// Now we don't need to keep localOutFlows info, since those are copied to modInfos.
				// so, clear it to reduce memory usage.
				distOutFlows.localOutFlows.outFlow.clear();

				numModule = modInfos.modInfoHash.size();

				dc.cout() << "The number of active modules: " << numModule << std::endl;

				dc.cout() << "sumAllExitPr = " << sumAllExitPr << std::endl;
				dc.cout() << "sumExit_log_sumExit = " << pLogP(sumAllExitPr) / log(2.0) << std::endl;
				dc.cout() << "- sum_size_log_size = " << - (sum_size_log_size / log(2.0)) << std::endl;
				dc.cout() << "- 2.0 * sum_exit_log_exit = " << - 2.0 * (new_sum_log_exit / log(2.0)) << std::endl;
				dc.cout() << "sum_stay_log_stay = " << new_sum_log_stay / log(2.0) << std::endl;

				currentMDL = (pLogP(sumAllExitPr) - 2.0 * new_sum_log_exit + new_sum_log_stay - sum_size_log_size) / log(2.0);
				dc.cout() << "Calculated Code Length = " << currentMDL << std::endl;


				if (currentMDL > prevMDL) {
					dc.cout() << "CurrentMDL is LARGER than previous MDL, so we will rollback to the previous module structure.\n";
					graph.transform_vertices(rollBackModID);	// This will revert the movement made by this vertex program.
					modInfoPtr = &oldModInfos;
					currentMDL = prevMDL;
					swapModInfos = false;
				}
				else {
					swapModInfos = true;
				}

				// Update module state ...
				graph.transform_vertices(updateModState);
			
			} while (repeatMode && (prevMDL - currentMDL > THRESHOLD));


			//////////////////////////////
			// Super Node Movement!!!  ///
			//////////////////////////////
			do {
				prevMDL = currentMDL;


				timer2.start();
				/////////////////////////
				// generate superNode  //
				/////////////////////////
				distEdgeMap.clearLocalEdgeMap();

				EdgeMap modEdgeMap;

#ifdef _OPENMP
#pragma omp parallel
#endif
				{
					EdgeMap tmpEdgeMap = EdgeMap();

#ifdef _OPENMP
					#pragma omp for
#endif
					// below is from distributed_aggregator.hpp and modified it.
					for (int i = 0; i <(int)graph.num_local_vertices(); ++i) {
						for (graph_type::local_edge_list_type::iterator it = graph.l_vertex(i).in_edges().begin(); it != graph.l_vertex(i).in_edges().end(); it++) {
							tmpEdgeMap.add(*it);
						}
					}

#ifdef _OPENMP
					#pragma omp critical
#endif
					{
						modEdgeMap += tmpEdgeMap;
					}
				}

				std::vector<EdgeMap> partEdgeMapToProc(numProcs);

				for (edgeList::iterator it = modEdgeMap.eList.begin(); it != modEdgeMap.eList.end(); it++) {
					int toProc = (it->first.first + it->first.second) % numProcs;

					partEdgeMapToProc[toProc].eList[it->first] = it->second;	// don't need to check key.
				}
	
				// copy my-own-part of edgeMap from what I currently have.
				distEdgeMap.localEdgeMap.eList.swap(partEdgeMapToProc[dc.procid()].eList);
				
				size_t target = (dc.procid()+1) % numProcs;
				size_t myID = dc.procid();

				dc.full_barrier();	// make sure that each localEdgeMap has its own part, here.

				while (target != myID) {
					// calling this will combine EdgeMap part belong to other procs.
					distEdgeMap.combineEdgeMap(target, partEdgeMapToProc[target]);
					target = (target + 1) % numProcs;
				}

				dc.full_barrier();

				modEdgeMap.eList.clear();
				for(std::vector<EdgeMap>::iterator it = partEdgeMapToProc.begin(); it != partEdgeMapToProc.end(); it++) {
					it->eList.clear();
				}
				partEdgeMapToProc.clear();


				aggrEdgeMapTime += timer2.current_time();


				timer2.start();
				graph_type sp_graph(dc, clopts);

				// generate sp_graph in parallel !!
				{
					int start, end;
					findAssignedPart(&start, &end, modInfoPtr->modInfoHash.size(), numProcs, dc.procid());

					int k = 0;
					std::map<int, ModuleInfo>::const_iterator it = modInfoPtr->modInfoHash.begin();
					for (k = 0; k < start; k++) {
						it++;
					}

					//for (std::map<int, ModuleInfo>::const_iterator it = modInfos.modInfoHash.begin(); it != modInfos.modInfoHash.end(); it++) {
					for (k = start; k < end && it != modInfoPtr->modInfoHash.end(); it++, k++) {
						sp_graph.add_vertex(it->first, node_data(it->second));
					}

					dc.full_barrier();		// make sure all the vertices are added before starting to add edges.
	
					edgeList& localEdgeList = distEdgeMap.localEdgeMap.eList;

					for (edgeList::const_iterator eit = localEdgeList.begin(); eit != localEdgeList.end(); eit++)
						sp_graph.add_edge(eit->first.first, eit->first.second, eit->second);
				}
				// don't need to keep localEdgeList, so clear it!
				distEdgeMap.localEdgeMap.eList.clear();



				spGraphGenTime += timer2.current_time();

				timer2.start();
				sp_graph.finalize();
				spGraphFinalizeTime += timer2.current_time();


				dc.cout() << "The number of vertices = " << sp_graph.num_vertices() << std::endl;
				dc.cout() << "The number of edges in sp_graph = " << sp_graph.num_edges() << std::endl;


				// After generating sp-graph, we can swap oldModInfos and modInfos...
				if (swapModInfos) {
					oldModInfos.modInfoHash.swap(modInfos.modInfoHash);
					modInfos.modInfoHash.clear();
				}


				/////////////////////////////////////
				// run findCommunity with sp-graph //
				/////////////////////////////////////
				dc.full_barrier();

				graphlab::omni_engine<findCommunity> engine5(dc, sp_graph, exec_type, clopts);
				engine5.signal_all();
				engine5.start();

				const double findCommTime_spGraph = engine5.elapsed_seconds();
				sumFindCommTime += findCommTime_spGraph;
				dc.cout() << "Running time of finding communities of the sp-graph: " << findCommTime_spGraph << " seconds.\n";

				dc.full_barrier();


				// call map-reduce-verte for aggregating all module-assignment information.
				ModIDs sp_mods = sp_graph.map_reduce_vertices<ModIDs>(aggregateModAssign);
				modAssignPtr = &sp_mods;

				graph.transform_vertices(updateModAssign);


				numModule = 0;
	
				new_sum_log_exit = 0.0;
				new_sum_log_stay = 0.0;
				sumAllExitPr = 0.0;		

				timer2.start();

				distModInfos.clearLocalModInfos();

				getDistModInfos(distModInfos, graph, dc);

				getModInfoTime += timer2.current_time();

				dc.full_barrier();

				timer2.start();
				distOutFlows.clearLocalOutFlowsMap();

				getDistOutFlows(distOutFlows, graph, dc);	

				combineOutFlowsTime += timer2.current_time();

				dc.full_barrier();


				timer2.start();

	
				sumAllExitPr = calculateExitProb(distOutFlows.localOutFlows.outFlow, distModInfos.localModInfos.modInfoHash, new_sum_log_exit, new_sum_log_stay);

				dc.full_barrier();

				exitCalcTime += timer2.current_time();

				// Here we have only partial modInfos and OutFlows...
				// So, we need to combine them all here!
				//modInfos = ModuleInfos();
				modInfos.modInfoHash.swap(distModInfos.localModInfos.modInfoHash);
				dc.all_reduce(modInfos);
				dc.all_reduce(sumAllExitPr);
				dc.all_reduce(new_sum_log_exit);
				dc.all_reduce(new_sum_log_stay);
	
				modInfoPtr = &modInfos;

				// Now we don't need to keep localOutFlows info, since those are copied to modInfos.
				// so, clear it to reduce memory usage.
				distOutFlows.localOutFlows.outFlow.clear();
				

				numModule = modInfos.modInfoHash.size();


				dc.cout() << "The number of active modules: " << numModule << std::endl;
	
				dc.cout() << "sumAllExitPr = " << sumAllExitPr << std::endl;
				dc.cout() << "sumExit_log_sumExit = " << pLogP(sumAllExitPr) / log(2.0) << std::endl;
				dc.cout() << "- sum_size_log_size = " << - (sum_size_log_size / log(2.0)) << std::endl;
				dc.cout() << "- 2.0 * sum_exit_log_exit = " << - 2.0 * (new_sum_log_exit / log(2.0)) << std::endl;
				dc.cout() << "sum_stay_log_stay = " << new_sum_log_stay / log(2.0) << std::endl;
	
				currentMDL = (pLogP(sumAllExitPr) - 2.0 * new_sum_log_exit + new_sum_log_stay - sum_size_log_size) / log(2.0);
				dc.cout() << "Calculated Code Length = " << currentMDL << std::endl;


				if (currentMDL > prevMDL) {
					dc.cout() << "CurrentMDL is LARGER than previous MDL, so we will rollback to the previous module structure.\n";
					graph.transform_vertices(rollBackModID);	// This will revert the movement made by this vertex program.
					
					modInfoPtr = &oldModInfos;		// rollback to the old modInfos.
					currentMDL = prevMDL;

					swapModInfos = false;
				}
				else {
					swapModInfos = true;
				}


				// Update module state ...
				graph.transform_vertices(updateModState);

			} while ((prevMDL - currentMDL) > THRESHOLD);

			dc.cout() << "SuperStep[" << superStep << "] - Code Length = " << currentMDL << ",\t";
			dc.cout() << "current elapsed time = " << timer.current_time() << std::endl;

			superStep++;

		} while (repeatOuterMode && ((outerPrevMDL - currentMDL) > outer_thresh));
		
		dc.cout() << "Elapsed time for finding Communities = " << timer.current_time() << std::endl;
		dc.cout() << "Sum of findCommTime = " << sumFindCommTime << std::endl;
		dc.cout() << "Sum of exitCalcTime = " << exitCalcTime << std::endl;
		dc.cout() << "Sum of getModInfoTime = " << getModInfoTime << std::endl;
		dc.cout() << "Sum of combineOutFlowsTime = " << combineOutFlowsTime << std::endl;
		dc.cout() << "Sum of aggrEdgeMapTime = " << aggrEdgeMapTime << std::endl;
		dc.cout() << "Sum of spGraphGenTime = " << spGraphGenTime << std::endl;
		dc.cout() << "Sum of spGraphFinalizeTime = " << spGraphFinalizeTime << std::endl;
		dc.cout() << "Final Calculated Code Length = " << currentMDL << std::endl;

		// Save the final community detection result...
		if (prefix != "") {
			std::stringstream ss;
			ss << i;

			//sp_graph.save(prefix+"_"+ss.str()+"_sp", relaxmap_writer(),
			//			false,	// do not gzip
			//			true,	// save vertices
			//			false);	// do not save edges

			graph.save(prefix+"_"+ss.str(), relaxmap_writer(),
						false,	// do not gzip
						true,	// save vertices
						false);	// do not save edges
		}
	}

	// Finalize mpi layer and quit...
	graphlab::mpi_tools::finalize();

	return EXIT_SUCCESS;
} // END OF MAIN.



// Miscellaneous function..
void findAssignedPart(int* start, int* end, int numAll, int numProc, int myID) {
	if (numAll % numProc == 0) {
		*start = (numAll / numProc) * myID;
		*end = (numAll/ numProc) * (myID + 1);
	}
	else {
		int block = numAll / numProc;
		int modular = numAll % numProc;

		if (myID < modular) {
			*start = (block + 1) * myID;
			*end = (block + 1) * (myID + 1);
		}
		else {
			*start = block * myID + modular;
			*end = block * (myID + 1) + modular;
		}
	}
}

/*
 * This function will generate distributed_OutFlows from local edges.
 * Each machine will have a part of OutFlows information.
 * computation and communication cost is much lower than using map_reduce_edges().
 */
void getDistOutFlows(distributed_OutFlows& distOutFlows, graph_type& graph, graphlab::distributed_control& dc) {
	int numProcs = dc.numprocs();
	size_t myID = dc.procid();

	OutFlows localOutFlows = OutFlows();

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		OutFlows tmpOutFlows = OutFlows();

#ifdef _OPENMP
		#pragma omp for
#endif
		for (int i = 0; i <(int)graph.num_local_vertices(); ++i) {
			for (graph_type::local_edge_list_type::iterator it = graph.l_vertex(i).in_edges().begin(); it != graph.l_vertex(i).in_edges().end(); it++) {
				//localOutFlows.addEdge(*it);
				tmpOutFlows.addEdge(*it);
			}
		}

#ifdef _OPENMP
		#pragma omp critical
#endif
		{
			localOutFlows += tmpOutFlows;
		}
	}

	std::vector<OutFlows> partOutFlowToProc(numProcs);

	for (std::map<int, double>::const_iterator it = localOutFlows.outFlow.begin(); it != localOutFlows.outFlow.end(); it++) {
		int toProc = it->first % numProcs;

		partOutFlowToProc[toProc].outFlow[it->first] = it->second;	// don't need to check key duplication.
	}

	// copy my-own-part of outFlow from what I currently have.
	distOutFlows.localOutFlows.outFlow.swap(partOutFlowToProc[myID].outFlow);

	size_t targetProc = (myID+1) % numProcs;

	dc.full_barrier();	// make sure that each localEdgeMap has its own part, here.

	while (targetProc != myID) {
		// calling this will combine OutFlows part assigned to myID currently belongs to other procs.
		distOutFlows.combineOutFlows(targetProc, partOutFlowToProc[targetProc]);
		targetProc = (targetProc + 1) % numProcs;
	}

	dc.full_barrier();
}






void getDistModInfos(distributed_ModInfos& distModInfos, graph_type& graph, graphlab::distributed_control& dc) {
	int numProcs = dc.numprocs();
	size_t myID = dc.procid();

	ModuleInfos localModInfos = ModuleInfos();

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		ModuleInfos tmpModInfos = ModuleInfos();

#ifdef _OPENMP
		#pragma omp for
#endif
		for (int i = 0; i <(int)graph.num_local_vertices(); ++i) {
			if (graph.l_vertex(i).owner() == myID) {
				//tmpModInfos.addVertex(graph.l_vertex(i));
				addVertexToModInfos(tmpModInfos, graph.l_vertex(i));
			}
		}

#ifdef _OPENMP
		#pragma omp critical
#endif
		{
			localModInfos += tmpModInfos;
		}
	}

	std::vector<ModuleInfos> partModInfoToProc(numProcs);

	for (std::map<int, ModuleInfo>::const_iterator it = localModInfos.modInfoHash.begin(); it != localModInfos.modInfoHash.end(); it++) {
		int toProc = it->first % numProcs;

		partModInfoToProc[toProc].modInfoHash[it->first] = it->second;	// don't need to check key duplication.
	}

	// copy my-own-part of outFlow from what I currently have.
	distModInfos.localModInfos.modInfoHash.swap(partModInfoToProc[myID].modInfoHash);

	size_t targetProc = (myID+1) % numProcs;

	dc.full_barrier();	// make sure that each localEdgeMap has its own part, here.

	while (targetProc != myID) {
		// calling this will combine EdgeMap part belong to other procs.
		distModInfos.combineModInfos(targetProc, partModInfoToProc[targetProc]);
		targetProc = (targetProc + 1) % numProcs;
	}

	dc.full_barrier();
}




double calculateExitProb(const std::map<int, double>& outFlow, std::map<int, ModuleInfo>& modInfo, double& new_sum_log_exit, double& new_sum_log_stay) {

	double sumAllExitPr = 0.0;
	new_sum_log_exit = 0.0;
	new_sum_log_stay = 0.0;

	for (std::map<int, double>::const_iterator it = outFlow.begin(); it != outFlow.end(); it++) {
		int mId = it->first;

		double exitPr = RESET_PROB * (1.0 - modInfo[mId].sumTPWeight) * modInfo[mId].sumPr;
		exitPr += (1.0 - RESET_PROB) * (it->second + (1.0 - modInfo[mId].sumTPWeight) * modInfo[mId].sumDangling);

		sumAllExitPr += exitPr;
		new_sum_log_exit += pLogP(exitPr);
		new_sum_log_stay += pLogP(exitPr + modInfo[mId].sumPr);

		//numModule++;

		//////////////////////////
		// Now update exitPr and outFlow of the corresponding modInfo.
		modInfo[mId].exitPr = exitPr;
		modInfo[mId].stayPr = exitPr + modInfo[mId].sumPr;
		modInfo[mId].modOutFlow = it->second;
	}

	return sumAllExitPr;
}
