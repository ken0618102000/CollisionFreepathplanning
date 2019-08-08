/*
https://github.com/yan-qi/k-shortest-paths-cpp-version
*/


#include "stdafx.h"
#include <set>
#include <map>
#include <queue>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "GraphElements.h"
#include "Graph.h"
#include "DijkstraShortestPathAlg.h"
#include "YenTopKShortestPathsAlg.h"

using namespace std;
const double Graph::DISCONNECT = (numeric_limits<double>::max)();


//----------------------DijkstraShortestPathAlg.cpp-------------------------
BasePath* DijkstraShortestPathAlg::get_shortest_path(BaseVertex* source, BaseVertex* sink)
{
	determine_shortest_paths(source, sink, true);

	std::vector<BaseVertex*> vertex_list;
	std::map<BaseVertex*, double>::const_iterator pos =
		m_mpStartDistanceIndex.find(sink);
	double weight = pos != m_mpStartDistanceIndex.end() ? pos->second : Graph::DISCONNECT;

	if (weight < Graph::DISCONNECT)
	{
		BaseVertex* cur_vertex_pt = sink;
		do
		{
			vertex_list.insert(vertex_list.begin(), cur_vertex_pt);

			std::map<BaseVertex*, BaseVertex*>::const_iterator pre_pos =
				m_mpPredecessorVertex.find(cur_vertex_pt);

			if (pre_pos == m_mpPredecessorVertex.end()) break;

			cur_vertex_pt = pre_pos->second;

		} while (cur_vertex_pt != source);

		vertex_list.insert(vertex_list.begin(), source);
	}
	return new BasePath(vertex_list, weight);
}

void DijkstraShortestPathAlg::determine_shortest_paths(BaseVertex* source, BaseVertex* sink, bool is_source2sink)
{
	//1. clear the intermediate variables
	clear();

	//2. initiate the local variables
	BaseVertex* end_vertex = is_source2sink ? sink : source;
	BaseVertex* start_vertex = is_source2sink ? source : sink;
	m_mpStartDistanceIndex[start_vertex] = 0;
	start_vertex->Weight(0);
	m_quCandidateVertices.insert(start_vertex);

	//3. start searching for the shortest path
	while (!m_quCandidateVertices.empty())
	{
		multiset<BaseVertex*, WeightLess<BaseVertex> >::const_iterator pos = m_quCandidateVertices.begin();

		BaseVertex* cur_vertex_pt = *pos; //m_quCandidateVertices.top();
		m_quCandidateVertices.erase(pos);

		if (cur_vertex_pt == end_vertex) break;

		m_stDeterminedVertices.insert(cur_vertex_pt->getID());

		improve2vertex(cur_vertex_pt, is_source2sink);
	}
}

void DijkstraShortestPathAlg::improve2vertex(BaseVertex* cur_vertex_pt, bool is_source2sink)
{
	// 1. get the neighboring vertices 
	set<BaseVertex*>* neighbor_vertex_list_pt = new set<BaseVertex*>();

	if (is_source2sink)
	{
		m_pDirectGraph->get_adjacent_vertices(cur_vertex_pt, *neighbor_vertex_list_pt);
	}
	else
	{
		m_pDirectGraph->get_precedent_vertices(cur_vertex_pt, *neighbor_vertex_list_pt);
	}

	// 2. update the distance passing on the current vertex
	for (set<BaseVertex*>::iterator cur_neighbor_pos = neighbor_vertex_list_pt->begin();
		cur_neighbor_pos != neighbor_vertex_list_pt->end(); ++cur_neighbor_pos)
	{
		//2.1 skip if it has been visited before
		if (m_stDeterminedVertices.find((*cur_neighbor_pos)->getID()) != m_stDeterminedVertices.end())
		{
			continue;
		}

		//2.2 calculate the distance
		map<BaseVertex*, double>::const_iterator cur_pos = m_mpStartDistanceIndex.find(cur_vertex_pt);
		double distance = cur_pos != m_mpStartDistanceIndex.end() ? cur_pos->second : Graph::DISCONNECT;

		distance += is_source2sink ? m_pDirectGraph->get_edge_weight(cur_vertex_pt, *cur_neighbor_pos) :
			m_pDirectGraph->get_edge_weight(*cur_neighbor_pos, cur_vertex_pt);

		//2.3 update the distance if necessary
		cur_pos = m_mpStartDistanceIndex.find(*cur_neighbor_pos);
		if (cur_pos == m_mpStartDistanceIndex.end() || cur_pos->second > distance)
		{
			m_mpStartDistanceIndex[*cur_neighbor_pos] = distance;
			m_mpPredecessorVertex[*cur_neighbor_pos] = cur_vertex_pt;

			(*cur_neighbor_pos)->Weight(distance);

			multiset<BaseVertex*, WeightLess<BaseVertex> >::const_iterator pos = m_quCandidateVertices.begin();
			for (; pos != m_quCandidateVertices.end(); ++pos)
			{
				if ((*pos)->getID() == (*cur_neighbor_pos)->getID())
				{
					break;
				}
			}
			if (pos != m_quCandidateVertices.end())
			{
				m_quCandidateVertices.erase(pos);
			}
			m_quCandidateVertices.insert(*cur_neighbor_pos);
		}
	}
	delete neighbor_vertex_list_pt;//------------fixed neighbor_vertex_list_pt memory leak-------------------------------------------------------
}

void DijkstraShortestPathAlg::clear()
{
	m_stDeterminedVertices.clear();
	m_mpPredecessorVertex.clear();
	m_mpStartDistanceIndex.clear();
	m_quCandidateVertices.clear();
}

BasePath* DijkstraShortestPathAlg::update_cost_forward(BaseVertex* vertex)
{
	double cost = Graph::DISCONNECT;

	// 1. get the set of successors of the input vertex
	set<BaseVertex*>* adj_vertex_set = new set<BaseVertex*>();
	m_pDirectGraph->get_adjacent_vertices(vertex, *adj_vertex_set);

	// 2. make sure the input vertex exists in the index
	map<BaseVertex*, double>::iterator pos4vertexInStartDistIndex = m_mpStartDistanceIndex.find(vertex);
	if (pos4vertexInStartDistIndex == m_mpStartDistanceIndex.end())
	{
		pos4vertexInStartDistIndex =
			(m_mpStartDistanceIndex.insert(make_pair(vertex, Graph::DISCONNECT))).first;
	}

	// 3. update the distance from the root to the input vertex if necessary
	for (set<BaseVertex*>::const_iterator pos = adj_vertex_set->begin(); pos != adj_vertex_set->end(); ++pos)
	{
		// 3.1 get the distance from the root to one successor of the input vertex
		map<BaseVertex*, double>::const_iterator cur_vertex_pos = m_mpStartDistanceIndex.find(*pos);
		double distance = cur_vertex_pos == m_mpStartDistanceIndex.end() ?
			Graph::DISCONNECT : cur_vertex_pos->second;

		// 3.2 calculate the distance from the root to the input vertex
		distance += m_pDirectGraph->get_edge_weight(vertex, *pos);

		// 3.3 update the distance if necessary 
		double cost_of_vertex = pos4vertexInStartDistIndex->second;
		if (cost_of_vertex > distance)
		{
			m_mpStartDistanceIndex[vertex] = distance;
			m_mpPredecessorVertex[vertex] = cur_vertex_pos->first;
			cost = distance;
		}
	}
	delete adj_vertex_set;  //------------fixed adj_vertex_set memory leak--------------------------------------------------------------
	// 4. create the sub_path if exists
	BasePath* sub_path = NULL;
	if (cost < Graph::DISCONNECT)
	{
		vector<BaseVertex*> vertex_list;
		vertex_list.push_back(vertex);

		map<BaseVertex*, BaseVertex*>::const_iterator pos4PredVertexMap =
			m_mpPredecessorVertex.find(vertex);

		while (pos4PredVertexMap != m_mpPredecessorVertex.end())
		{
			BaseVertex* pred_vertex_pt = pos4PredVertexMap->second;
			vertex_list.push_back(pred_vertex_pt);
			pos4PredVertexMap = m_mpPredecessorVertex.find(pred_vertex_pt);
		}

		sub_path = new BasePath(vertex_list, cost);
	}
	return sub_path;
}

void DijkstraShortestPathAlg::correct_cost_backward(BaseVertex* vertex)
{
	// 1. initialize the list of vertex to be updated
	vector<BaseVertex*> vertex_pt_list;
	vertex_pt_list.push_back(vertex);

	// 2. update the cost of relevant precedents of the input vertex
	while (!vertex_pt_list.empty())
	{
		BaseVertex* cur_vertex_pt = *(vertex_pt_list.begin());
		vertex_pt_list.erase(vertex_pt_list.begin());

		double cost_of_cur_vertex = m_mpStartDistanceIndex[cur_vertex_pt];

		set<BaseVertex*> pre_vertex_set;
		m_pDirectGraph->get_precedent_vertices(cur_vertex_pt, pre_vertex_set);
		for (set<BaseVertex*>::const_iterator pos = pre_vertex_set.begin(); pos != pre_vertex_set.end(); ++pos)
		{
			map<BaseVertex*, double>::const_iterator pos4StartDistIndexMap =
				m_mpStartDistanceIndex.find(*pos);
			double cost_of_pre_vertex = m_mpStartDistanceIndex.end() == pos4StartDistIndexMap ?
				Graph::DISCONNECT : pos4StartDistIndexMap->second;

			double fresh_cost = cost_of_cur_vertex + m_pDirectGraph->get_edge_weight(*pos, cur_vertex_pt);
			if (cost_of_pre_vertex > fresh_cost)
			{
				m_mpStartDistanceIndex[*pos] = fresh_cost;
				m_mpPredecessorVertex[*pos] = cur_vertex_pt;
				vertex_pt_list.push_back(*pos);
			}
		}
	}
}

//--------------------------Graph.cpp----------------------
Graph::Graph(const string& file_name)
{
	_import_from_file(file_name);
}

Graph::Graph(const Graph& graph)
{
	m_nVertexNum = graph.m_nVertexNum;
	m_nEdgeNum = graph.m_nEdgeNum;
	m_vtVertices.assign(graph.m_vtVertices.begin(), graph.m_vtVertices.end());
	m_mpFaninVertices.insert(graph.m_mpFaninVertices.begin(), graph.m_mpFaninVertices.end());
	m_mpFanoutVertices.insert(graph.m_mpFanoutVertices.begin(), graph.m_mpFanoutVertices.end());
	m_mpEdgeCodeWeight.insert(graph.m_mpEdgeCodeWeight.begin(), graph.m_mpEdgeCodeWeight.end());
	m_mpVertexIndex.insert(graph.m_mpVertexIndex.begin(), graph.m_mpVertexIndex.end());
}

Graph::Graph(int input_m_nVertexNum, int input_m_nEdgeNum, vector <CPoint> &input_CPoint_savepoint1, vector <CPoint> &input_CPoint_savepoint2, vector <CPoint> &input_all_point_map_original)
{
	int start_vertex, end_vertex;
	int loop1, loop2;
	double edge_weight;
	int vertex_id = 0;

	m_nVertexNum = input_m_nVertexNum;

	for (loop1 = 0; loop1 < input_m_nEdgeNum; loop1++)
	{
		for (loop2 = 0; loop2 < input_m_nVertexNum; loop2++) //¬d¸ß½s¸¹
		{

			if (input_CPoint_savepoint1[loop1] == input_all_point_map_original[loop2])
				start_vertex = loop2;

			if (input_CPoint_savepoint2[loop1] == input_all_point_map_original[loop2])
				end_vertex = loop2;
		}
		edge_weight = sqrt(pow(input_CPoint_savepoint1[loop1].x - input_CPoint_savepoint2[loop1].x, 2) + pow(input_CPoint_savepoint1[loop1].y - input_CPoint_savepoint2[loop1].y, 2));


		///3.2.1 construct the vertices
		BaseVertex* start_vertex_pt = get_vertex(start_vertex);
		BaseVertex* end_vertex_pt = get_vertex(end_vertex);
		BaseVertex* start_vertex_pt2 = get_vertex(end_vertex);
		BaseVertex* end_vertex_pt2 = get_vertex(start_vertex);

		///3.2.2 add the edge weight
		//// note that the duplicate edge would overwrite the one occurring before. 
		m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_weight;
		m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt2, end_vertex_pt2)] = edge_weight;

		///3.2.3 update the fan-in or fan-out variables
		//// Fan-in
		get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);
		get_vertex_set_pt(end_vertex_pt2, m_mpFaninVertices)->insert(start_vertex_pt2);

		//// Fan-out
		get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);
		get_vertex_set_pt(start_vertex_pt2, m_mpFanoutVertices)->insert(end_vertex_pt2);

	}

	if (m_nVertexNum != m_vtVertices.size())
	{
		cerr << "The number of nodes in the graph is " << m_vtVertices.size() << " instead of " << m_nVertexNum << endl;
		exit(1);
	}

	m_nVertexNum = m_vtVertices.size();
	m_nEdgeNum = m_mpEdgeCodeWeight.size();
}

Graph::~Graph(void)
{
//	clear();
}

///////////////////////////////////////////////////////////////////////////////
///  public  _import_from_file
///  Construct the graph by importing the edges from the input file. 
///
///  @param [in]       file_name const std::string &    The input graph file
///
///  This function doesn't return a value
///
///  @remarks The format of the file is as follows:
///   1. The first line has an integer as the number of vertices of the graph
///   2. Each line afterwards contains a directed edge in the graph:
///		     starting point, ending point and the weight of the edge. 
///		 These values are separated by 'white space'.
///
///  @see <TODO: insert text here>
///
///  @author Yan Qi @date 5/29/2010
///////////////////////////////////////////////////////////////////////////////
void Graph::_import_from_file(const string& input_file_name)
{
	const char* file_name = input_file_name.c_str();

	//1. Check the validity of the file
	ifstream ifs(file_name);
	if (!ifs)
	{
		cerr << "The file " << file_name << " can not be opened!" << endl;
		exit(1);
	}

	//2. Reset the members of the class
	clear();

	//3. Start to read information from the input file. 
	/// Note the format of the data in the graph file.
	//3.1 The first line has an integer as the number of vertices of the graph
	ifs >> m_nVertexNum;

	//3.2 In the following lines, each line contains a directed edge in the graph:
	///   the id of starting point, the id of ending point, the weight of the edge. 
	///   These values are separated by 'white space'. 
	int start_vertex, end_vertex;
	double edge_weight;
	int vertex_id = 0;

	while (ifs >> start_vertex)
	{
		if (start_vertex == -1)
		{
			break;
		}
		ifs >> end_vertex;
		ifs >> edge_weight;

		///3.2.1 construct the vertices
		BaseVertex* start_vertex_pt = get_vertex(start_vertex);
		BaseVertex* end_vertex_pt = get_vertex(end_vertex);

		///3.2.2 add the edge weight
		//// note that the duplicate edge would overwrite the one occurring before. 
		m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_weight;

		///3.2.3 update the fan-in or fan-out variables
		//// Fan-in
		get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);

		//// Fan-out
		get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);

	}

	if (m_nVertexNum != m_vtVertices.size())
	{
		cerr << "The number of nodes in the graph is " << m_vtVertices.size() << " instead of " << m_nVertexNum << endl;
		exit(1);
	}

	m_nVertexNum = m_vtVertices.size();
	m_nEdgeNum = m_mpEdgeCodeWeight.size();

	ifs.close();
}

BaseVertex* Graph::get_vertex(int node_id)
{
	if (m_stRemovedVertexIds.find(node_id) != m_stRemovedVertexIds.end())
	{
		return NULL;
	}
	else
	{
		BaseVertex* vertex_pt = NULL;
		const map<int, BaseVertex*>::iterator pos = m_mpVertexIndex.find(node_id);
		if (pos == m_mpVertexIndex.end())
		{
			int vertex_id = m_vtVertices.size();
			vertex_pt = new BaseVertex();
			vertex_pt->setID(node_id);
			m_mpVertexIndex[node_id] = vertex_pt;

			m_vtVertices.push_back(vertex_pt);
		}
		else
		{
			vertex_pt = pos->second;
		}

		return vertex_pt;
	}
}

void Graph::clear()
{
	m_nEdgeNum = 0;
	m_nVertexNum = 0;

	for (map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos = m_mpFaninVertices.begin();
		pos != m_mpFaninVertices.end(); ++pos)
	{
		delete pos->second;
	}
	m_mpFaninVertices.clear();

	for (map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos = m_mpFanoutVertices.begin();
		pos != m_mpFanoutVertices.end(); ++pos)
	{
		delete pos->second;
	}
	m_mpFanoutVertices.clear();


	m_mpEdgeCodeWeight.clear();

	//clear the list of vertices objects
	for_each(m_vtVertices.begin(), m_vtVertices.end(), DeleteFunc<BaseVertex>());
	m_vtVertices.clear();
	m_mpVertexIndex.clear();

	m_stRemovedVertexIds.clear();
	m_stRemovedEdge.clear();
}

int Graph::get_edge_code(const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt) const
{
	/// Note that the computation below works only if 
	/// the result is smaller than the maximum of an integer!
	return start_vertex_pt->getID()*m_nVertexNum + end_vertex_pt->getID();
}


set<BaseVertex*>* Graph::get_vertex_set_pt(BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index)
{
	BaseVertexPt2SetMapIterator pos = vertex_container_index.find(vertex_);

	if (pos == vertex_container_index.end())
	{
		set<BaseVertex*>* vertex_set = new set<BaseVertex*>();
		pair<BaseVertexPt2SetMapIterator, bool> ins_pos =
			vertex_container_index.insert(make_pair(vertex_, vertex_set));

		pos = ins_pos.first;
	}

	return pos->second;
}


double Graph::get_edge_weight(const BaseVertex* source, const BaseVertex* sink)
{
	int source_id = source->getID();
	int sink_id = sink->getID();

	if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
		|| m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
		|| m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
	{
		return DISCONNECT;
	}
	else
	{
		return get_original_edge_weight(source, sink);
	}
}


void Graph::get_adjacent_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set)
{
	int starting_vt_id = vertex->getID();

	if (m_stRemovedVertexIds.find(starting_vt_id) == m_stRemovedVertexIds.end())
	{
		set<BaseVertex*>* vertex_pt_set = get_vertex_set_pt(vertex, m_mpFanoutVertices);
		for (set<BaseVertex*>::const_iterator pos = (*vertex_pt_set).begin();
			pos != (*vertex_pt_set).end(); ++pos)
		{
			int ending_vt_id = (*pos)->getID();
			if (m_stRemovedVertexIds.find(ending_vt_id) != m_stRemovedVertexIds.end()
				|| m_stRemovedEdge.find(make_pair(starting_vt_id, ending_vt_id)) != m_stRemovedEdge.end())
			{
				continue;
			}
			//
			vertex_set.insert(*pos);
		}
	}
}

void Graph::get_precedent_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set)
{
	if (m_stRemovedVertexIds.find(vertex->getID()) == m_stRemovedVertexIds.end())
	{
		int ending_vt_id = vertex->getID();
		set<BaseVertex*>* pre_vertex_set = get_vertex_set_pt(vertex, m_mpFaninVertices);
		for (set<BaseVertex*>::const_iterator pos = (*pre_vertex_set).begin();
			pos != (*pre_vertex_set).end(); ++pos)
		{
			int starting_vt_id = (*pos)->getID();
			if (m_stRemovedVertexIds.find(starting_vt_id) != m_stRemovedVertexIds.end()
				|| m_stRemovedEdge.find(make_pair(starting_vt_id, ending_vt_id)) != m_stRemovedEdge.end())
			{
				continue;
			}
			//
			vertex_set.insert(*pos);
		}
	}
}

double Graph::get_original_edge_weight(const BaseVertex* source, const BaseVertex* sink)
{
	map<int, double>::const_iterator pos =
		m_mpEdgeCodeWeight.find(get_edge_code(source, sink));

	if (pos != m_mpEdgeCodeWeight.end())
	{
		return pos->second;
	}
	else
	{
		return DISCONNECT;
	}
}

//-------------YenTopKShortestPathsAlg.cpp----------------------------
void YenTopKShortestPathsAlg::clear()
{
	m_nGeneratedPathNum = 0;
	//------------fixed m_mpDerivationVertexIndex, m_vResultList, m_quPathCandidates memory leak---------------------------------------------
	map<BasePath*, BaseVertex*>().swap(m_mpDerivationVertexIndex);
	vector<BasePath*>().swap(m_vResultList);
	multiset<BasePath*, WeightLess<BasePath> >().swap(m_quPathCandidates);

	//	m_mpDerivationVertexIndex.clear();
	//	m_vResultList.clear();
	//	m_quPathCandidates.clear();
}

void YenTopKShortestPathsAlg::clearAll() //------------fixed clear memory leak----------------------------
{
	m_nGeneratedPathNum = 0;
	m_pGraph->clear();
	//------------fixed m_mpDerivationVertexIndex, m_vResultList, m_quPathCandidates memory leak---------------------------------------------
	map<BasePath*, BaseVertex*>().swap(m_mpDerivationVertexIndex);
	vector<BasePath*>().swap(m_vResultList);
	multiset<BasePath*, WeightLess<BasePath> >().swap(m_quPathCandidates);

	//	m_mpDerivationVertexIndex.clear();
	//	m_vResultList.clear();
	//	m_quPathCandidates.clear();
}

void YenTopKShortestPathsAlg::_init()
{
	clear();
	if (m_pSourceVertex != NULL && m_pTargetVertex != NULL)
	{
		BasePath* pShortestPath = get_shortest_path(m_pSourceVertex, m_pTargetVertex);
		if (pShortestPath != NULL && pShortestPath->length() > 1)
		{
			m_quPathCandidates.insert(pShortestPath);
			m_mpDerivationVertexIndex[pShortestPath] = m_pSourceVertex;
		}
	}
}

BasePath* YenTopKShortestPathsAlg::get_shortest_path(BaseVertex* pSource, BaseVertex* pTarget)
{
	DijkstraShortestPathAlg dijkstra_alg(m_pGraph);
	return dijkstra_alg.get_shortest_path(pSource, pTarget);
}

bool YenTopKShortestPathsAlg::has_next()
{
	return !m_quPathCandidates.empty();
}

BasePath* YenTopKShortestPathsAlg::next()
{
	//1. Prepare for removing vertices and arcs
	BasePath* cur_path = *(m_quPathCandidates.begin());//m_quPathCandidates.top();

	//m_quPathCandidates.pop();
	m_quPathCandidates.erase(m_quPathCandidates.begin());
	m_vResultList.push_back(cur_path);

	int count = m_vResultList.size();

	BaseVertex* cur_derivation_pt = m_mpDerivationVertexIndex.find(cur_path)->second;
	vector<BaseVertex*> sub_path_of_derivation_pt;
	cur_path->SubPath(sub_path_of_derivation_pt, cur_derivation_pt);
	int sub_path_length = sub_path_of_derivation_pt.size();

	//2. Remove the vertices and arcs in the graph
	for (int i = 0; i < count - 1; ++i)
	{
		BasePath* cur_result_path = m_vResultList.at(i);
		vector<BaseVertex*> cur_result_sub_path_of_derivation_pt;

		if (!cur_result_path->SubPath(cur_result_sub_path_of_derivation_pt, cur_derivation_pt)) continue;

		if (sub_path_length != cur_result_sub_path_of_derivation_pt.size()) continue;

		bool is_equal = true;
		for (int i = 0; i < sub_path_length; ++i)
		{
			if (sub_path_of_derivation_pt.at(i) != cur_result_sub_path_of_derivation_pt.at(i))
			{
				is_equal = false;
				break;
			}
		}
		if (!is_equal) continue;

		//
		BaseVertex* cur_succ_vertex = cur_result_path->GetVertex(sub_path_length + 1);
		m_pGraph->remove_edge(make_pair(cur_derivation_pt->getID(), cur_succ_vertex->getID()));
	}

	//2.1 remove vertices and edges along the current result
	int path_length = cur_path->length();
	for (int i = 0; i < path_length - 1; ++i)
	{
		m_pGraph->remove_vertex(cur_path->GetVertex(i)->getID());
		m_pGraph->remove_edge(make_pair(
			cur_path->GetVertex(i)->getID(), cur_path->GetVertex(i + 1)->getID()));
	}

	//3. Calculate the shortest tree rooted at target vertex in the graph
	DijkstraShortestPathAlg reverse_tree(m_pGraph);
	reverse_tree.get_shortest_path_flower(m_pTargetVertex);

	//4. Recover the deleted vertices and update the cost and identify the new candidates results
	bool is_done = false;
	for (int i = path_length - 2; i >= 0 && !is_done; --i)
	{
		//4.1 Get the vertex to be recovered
		BaseVertex* cur_recover_vertex = cur_path->GetVertex(i);
		m_pGraph->recover_removed_vertex(cur_recover_vertex->getID());

		//4.2 Check if we should stop continuing in the next iteration
		if (cur_recover_vertex->getID() == cur_derivation_pt->getID())
		{
			is_done = true;
		}

		//4.3 Calculate cost using forward star form
		BasePath* sub_path = reverse_tree.update_cost_forward(cur_recover_vertex);

		//4.4 Get one candidate result if possible
		if (sub_path != NULL)
		{
			++m_nGeneratedPathNum;

			//4.4.1 Get the prefix from the concerned path
			double cost = 0;
			reverse_tree.correct_cost_backward(cur_recover_vertex);

			vector<BaseVertex*> pre_path_list;
			for (int j = 0; j < path_length; ++j)
			{
				BaseVertex* cur_vertex = cur_path->GetVertex(j);
				if (cur_vertex->getID() == cur_recover_vertex->getID())
				{
					//j = path_length;
					break;
				}
				else
				{
					cost += m_pGraph->get_original_edge_weight(
						cur_path->GetVertex(j), cur_path->GetVertex(1 + j));
					pre_path_list.push_back(cur_vertex);
				}
			}
			//
			for (int j = 0; j < sub_path->length(); ++j)
			{
				pre_path_list.push_back(sub_path->GetVertex(j));
			}

			//4.4.2 Compose a candidate
			sub_path = new Path(pre_path_list, cost + sub_path->Weight());

			//4.4.3 Put it in the candidate pool if new
			if (m_mpDerivationVertexIndex.find(sub_path) == m_mpDerivationVertexIndex.end())
			{
				m_quPathCandidates.insert(sub_path);
				m_mpDerivationVertexIndex[sub_path] = cur_recover_vertex;
			}
		}

		//4.5 Restore the edge
		BaseVertex* succ_vertex = cur_path->GetVertex(i + 1);
		m_pGraph->recover_removed_edge(make_pair(cur_recover_vertex->getID(), succ_vertex->getID()));

		//4.6 Update cost if necessary
		double cost_1 = m_pGraph->get_edge_weight(cur_recover_vertex, succ_vertex)
			+ reverse_tree.get_start_distance_at(succ_vertex);

		if (reverse_tree.get_start_distance_at(cur_recover_vertex) > cost_1)
		{
			reverse_tree.set_start_distance_at(cur_recover_vertex, cost_1);
			reverse_tree.set_predecessor_vertex(cur_recover_vertex, succ_vertex);
			reverse_tree.correct_cost_backward(cur_recover_vertex);
		}
	}

	//5. Restore everything
	m_pGraph->recover_removed_edges();
	m_pGraph->recover_removed_vertices();

	return cur_path;
}

void YenTopKShortestPathsAlg::get_shortest_paths(BaseVertex* pSource,
	BaseVertex* pTarget, int top_k, vector<BasePath*>& result_list)
{
	m_pSourceVertex = pSource;
	m_pTargetVertex = pTarget;

	_init();
	int count = 0;
	while (has_next() && count < top_k)
	{
		next();
		++count;
	}

	result_list.assign(m_vResultList.begin(), m_vResultList.end());
}