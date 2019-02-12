#include "instance.hpp"
#include "model.hpp"
#include "edge.hpp"
#include "hpath.hpp"
#include "callbacks.hpp"
#include <algorithm>

HamiltonianPath::HamiltonianPath(const Instance & instance) :
    _instance(instance), _model(Model()), 
    _vertexCoords(), _isTarget(), _isSatellite(), 
    _vertexToTargetMap(), _vertexToSatelliteMap(), 
    _source(0), _destination(1), _edges(), 
    _edgeMap(), _path(), _pathCost(0.0) {};


void HamiltonianPath::populatePathData(std::vector<int> & targets, 
    std::vector<int> satellites, 
    int sourceTarget, 
    int destinationTarget) {

    auto targetCoords = getInstance().getTargetCoords();
    auto satelliteCoords = getInstance().getSatelliteCoords();
    std::unordered_map<int, std::tuple<double, double>> vertexCoords; 
    std::unordered_map<int, int> vertexToTargetMap;
    std::unordered_map<int, int> vertexToSatelliteMap;

    int numVertices = targets.size() + satellites.size();
    std::unordered_map<int, bool> isTarget, isSatellite;
    for (int i=0; i<numVertices; ++i) {
        isTarget.insert({i, false});
        isSatellite.insert({i, false});
    }

    int vertexIndex = 0;
    for (int target : targets) {
        if (target == sourceTarget) HamiltonianPath::setSource(vertexIndex);
        if (target == destinationTarget) HamiltonianPath::setDestination(vertexIndex);
        vertexCoords.insert({vertexIndex, targetCoords.at(target)});
        isTarget[vertexIndex] = true;
        vertexToTargetMap.insert({vertexIndex, target});
        vertexIndex++;
    }

    for (int satellite : satellites) {
        vertexCoords.insert({vertexIndex, satelliteCoords.at(satellite)});
        isSatellite[vertexIndex] = true;
        vertexToSatelliteMap.insert({vertexIndex, satellite});
        vertexIndex++;
    }

    setVertexCoords(vertexCoords);
    setIsTarget(isTarget);
    setIsSatellite(isSatellite);
    setVertexToTargetMap(vertexToTargetMap);
    setVertexToSatelliteMap(vertexToSatelliteMap);
    return;
};

void HamiltonianPath::createEdges() {

    auto vertexCoords = getVertexCoords();
    std::vector<Edge> edges;
    std::unordered_map<std::tuple<int, int>, int> edgeMap;

    for (int i=0; i<getNumVertices(); ++i) {
        for (int j=i+1; j<getNumVertices(); ++j) {
            auto ith = vertexCoords.at(i);
            auto jth = vertexCoords.at(j);
            if (i == HamiltonianPath::getSource() && j ==HamiltonianPath:: getDestination()) 
                continue;
            double cost = std::hypot(std::get<0>(ith) - std::get<0>(jth), std::get<1>(ith) - std::get<1>(jth));
            edges.push_back(Edge(i, j, cost));
            edgeMap.insert({std::make_tuple(i, j), edges.size()-1});
            edgeMap.insert({std::make_tuple(j, i), edges.size()-1});
        }
    }

    setEdges(edges);
    setEdgeMap(edgeMap);
    return;

};

void  HamiltonianPath::addVariablesandConstraints() { 

    Model model;
    
    model.getModel().setName(getInstance().getName().c_str());
    IloNumVarArray x(model.getEnv());
    model.getVariables().insert({"x", x});
    
    
    //IloNumArray cost(model.getEnv());

    for (auto const & edge : _edges) {
        std::string varname = "x_" + std::to_string(edge.from()) +  "_" + std::to_string(edge.to());
        //cost.add(edge.cost());
        IloNumVar var(model.getEnv(), 0, 1, ILOINT, varname.c_str());
        model.getVariables().at("x").add(var);
    }
   
    
    for (int i =0; i<getNumVertices(); ++i){
        IloExpr rowconstraint(model.getEnv());
        std::set<int> index ;
        for (int j=0; j < getNumVertices(); ++j) {
            auto itr = _edgeMap.find(std::make_tuple(i,j));
            if (_edgeMap.end() != itr)       
            index.insert(itr->second);
        }
            for ( auto i : index)
                rowconstraint += model.getVariables().at("x")[i];
    int EQ =2;
    if( i == HamiltonianPath::getSource()|| i == HamiltonianPath:: getDestination())
    EQ =1;
    IloRange con1(model.getEnv(),EQ,rowconstraint,EQ);
    model.getConstraints().push_back(con1);
    model.getModel().add(con1);
    }
    //model.getModel().add(IloMinimize(model.getEnv(), IloScalProd(cost, model.getVariables().at("x"))));
    setModel(model);
    return;
    
};


void HamiltonianPath::solve(){

addVariablesandConstraints();

IloNumArray cost(_model.getEnv());
    for (int i=0; i< _edges.size(); ++i) {
        auto edge = _edges[i];
        cost.add(edge.cost());
    }

    _model.getModel().add(
        IloMinimize(_model.getEnv(), 
        IloScalProd(cost, _model.getVariables().at("x"))
    ));


IloCplex cplex(_model.getEnv());
    cplex.extract(_model.getModel());
    cplex.exportModel("test.lp");
    cplex.setOut(_model.getEnv().getNullStream());
    cplex.use(addLazyCallback(_model.getEnv(), 
        _model,
        _edges, 
        _edgeMap, 
        _vertexCoords.size(), 
        _source, 
        _destination
    ));
    cplex.solve();
    //std::ofstream outfile;//added
    //outfile.open("data.txt");//added   
    std::vector<std::tuple<int, int>> solutionEdges;
    std::vector<int> omit;
    
    
    for (int i=0; i<_model.getVariables().at("x").getSize(); ++i) {
         if (cplex.getValue(_model.getVariables().at("x")[i]) > 0.9) {
            int from = _edges[i].from();
            int to = _edges[i].to();
            //outfile<<from<<" " <<to<<std::endl; //added
            solutionEdges.push_back(std::make_pair(from, to));
            _pathCost += _edges[i].cost();
          
        }
    }
//outfile.close();//added
_path.push_back(_source);
int j =_source;
while (_path.size()<= solutionEdges.size())
{ 
for (int p =0;p <solutionEdges.size();++p){
        auto src = std::find(omit.begin(),omit.end(), p);
        
        if (src != omit.end())
        continue;

if(std::get<0>(solutionEdges[p]) == j){
              _path.push_back(std::get<1>(solutionEdges[p]));
              j = std:: get<1>(solutionEdges[p]);
              omit.push_back(p);
              break;
         }
         else if( std:: get<1>(solutionEdges[p]) == j){
              _path.push_back(std::get<0>(solutionEdges[p]));
              j = std:: get<0>(solutionEdges[p]);
              omit.push_back(p);
              break;
         }
}
}   
//std::cout<<"GETPATH size:"<<getPath().size()<<std::endl;
return;    
};




