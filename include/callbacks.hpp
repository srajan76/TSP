#ifndef CALLBACKS_HPP
#define CALLBACKS_HPP

#pragma once

#include <ilcplex/ilocplex.h>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <string>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <set>


typedef const std::unordered_map<std::tuple<int, int>, int> & edgeMapType;

std::vector<IloRange> generateLazyConstraints(
    Model & model,  
    const std::vector<Edge> & edges, 
    const std::unordered_map<std::tuple<int, int>, int> & edgeMap,
    const int numVertices,
    const int source, 
    const int destination,
    std::unordered_map<std::string, IloNumArray> & variableValues) {
 std::vector<IloRange> constraints;
    lemon::ListGraph supportGraph;
    for (int i=0; i<numVertices; ++i) supportGraph.addNode();

    IloNumArray xVals = variableValues.at("x");

    for (auto i=0; i<xVals.getSize(); ++i)
        if (xVals[i] > 1E-5) {
            int from =edges[i].from();
            int to = edges[i].to();
            supportGraph.addEdge(supportGraph.nodeFromId(from), supportGraph.nodeFromId(to));
        }

    lemon::ListGraph::NodeMap<int> componentMap(supportGraph);
    int numComponents = connectedComponents(supportGraph, componentMap);

if (numComponents == 1)
        return constraints;
    
    std::vector<std::set<int>> components(numComponents);
	for (lemon::ListGraph::NodeIt n(supportGraph); n!=lemon::INVALID; ++n)
		components[componentMap[n]].insert(supportGraph.id(n));
/* for (auto & component : components) {
        auto src = component.find(source);
        auto dst = component.find(destination);
        if (src != component.end()|| dst != component.end())
        continue;
        std::set<tuple<int,int>> varID;
        for(auto i=component.begin();i!=component.end();i++){
          for(auto j=std::next(i,1);j!=component.end();j++){
        auto temp = std::make_tuple(*i,*j);
        varID.insert(temp);
          }
        }
for (auto const & i: varID)
  expr + = model.getVariables().at("x")[edgeMap.at(i)];*/
    for (auto & component : components) {
        if(component.size()<=1)
        continue;
        auto src = component.find(source);
        auto dst = component.find(destination);
        if (src != component.end()|| dst != component.end())
        continue;
        IloExpr expr(model.getEnv());
        std::set<std::tuple<int,int>> varID;
        for(auto i=component.begin();i!=component.end();i++){
          for(auto j=std::next(i,1);j!=component.end();j++){
        auto temp = std::make_tuple(*i,*j);
        varID.insert(temp);
          }
        }
        for (auto i: varID){
        expr += model.getVariables().at("x")[edgeMap.at(i)];
        }
        IloRange constr(model.getEnv(), expr, component.size()-1);
        constraints.push_back(constr);
        expr.end();
    }
    
    return constraints;
}; 

ILOLAZYCONSTRAINTCALLBACK6(addLazyCallback, 
    Model&, model, 
    std::vector<Edge>&, edges, 
    edgeMapType, edgeMap,
    int, numVertices,
    int, source, 
    int, destination) {
    
    IloEnv env = getEnv();
    std::unordered_map<std::string, IloNumArray> variableValues;
    variableValues.insert(std::make_pair("x", IloNumArray(env, model.getVariables().at("x").getSize())));
    getValues(variableValues.at("x"), model.getVariables().at("x"));

    std::vector<IloRange> constraints = generateLazyConstraints(
        model, 
        edges, 
        edgeMap,
        numVertices, 
        source, 
        destination, 
        variableValues);

    for (IloRange constr : constraints)
        add(constr);

    variableValues.at("x").end();
    constraints.clear();
};


#endif